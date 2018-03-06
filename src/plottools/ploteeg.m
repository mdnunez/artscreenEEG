function ploteeg(datain,varargin)
%PLOTEEG - Plots sample*channel*trial EEG data
%
%Useage:  ploteeg(datain,varargin);
%
%INPUTS:  datain.data - sample*channel*trial EEG data
%         datain.sr - sample rate
%         datain.marker - 2*channel*trial
%
% OPTIONAL ARGUMENTS: (Passed in Name-Value Pairs)
%         chanlabels - 1*channel cell array of channel label strings
%                   Default: datain.label or 1:nchans
%
%         spacer - amplitude spacer between channels
%                  Default: standard deviation of all samples
%
%      varargin - Any "plot" inputs after the first two,
%                 i.e., plot(X,Y,varargin);
%
%See also: POWERSPEC

%% Copyright

% Copyright (C) 2018 Michael D. Nunez, <mdnunez1@uci.edu>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  01/04/18      Michael Nunez     Original code adapted from work by Beth A. Lopour
%  01/10/18      Michael Nunez              Add markers
%  02/13/18      Michael Nunez             Printing description of keys available
%  02/14/18      Michael Nunez           Add channel labels and varargin for plots

%To do:
% 1) Export markers

%% Frequency interval

nsamps = size(datain.data,1);
nchans = size(datain.data,2);
ntrials = size(datain.data,3);

[varargin,chanlabels,spacer] = parsevar(varargin,'chanlabels',[],'spacer',std(datain.data(:)));

if isempty(chanlabels)
    if isfield(datain,'label')
        chanlabels = datain.label;
    else
        chanlabels = sprintfc('%d',1:nchans); %Undocumented internal MATLAB function
    end
end


fprintf('Plotting the data...\n');
fprintf('The key commands available are left arrow and right arrow and ''e'' key...\n');

% Create figure window
screensize=get(0,'ScreenSize');
fig_h = figure('pos',[50 50 screensize(3)-100 screensize(4)-150], ...
    'DeleteFcn',@Exit_callback);
set(fig_h,'KeyPressFcn',@plotData);  % This function defines action for the keyboard

exitkey = 0;
plotepoch = 1;
epochtime = nsamps/datain.sr;
xtick = 1/datain.sr;    
plotchans = 1:nchans;

spacing = [0:(nchans-1)]*spacer;
matspacing = repmat(spacing,[nsamps 1]);

time = (epochtime*(plotepoch-1)+xtick):xtick:(epochtime*plotepoch);
cortplotx(time,squeeze(datain.data(:,plotchans,plotepoch))+matspacing(:,1:length(plotchans)),'custom', chanlabels(plotchans), varargin{:});
set(gca,'YTick',spacing(plotchans),'YTickLabel',chanlabels(plotchans));
set(gca,'YLim',[spacing(plotchans(1))-spacer*2 spacing(plotchans(end))+spacer])
xlabel('Time (sec)');

% Plot horizontal green line for any markers
if isfield(datain,'markers')
    markers = datain.markers;
    markers(markers == 0) = NaN;
    hold on;
    plot(time,squeeze(markers(:,plotchans,plotepoch))+matspacing(:,plotchans), 'g', 'linewidth', 2);
    hold off;
end

%% Nested Functions: These share a workspace with the main function********
% *************************************************************************
% *************************************************************************
    function plotData(~,eventDat)
    % This function is executed on a key press in the figure window; it
    % increments the plot time forward for a 'rightarrow', backward
    % for a 'leftarrow'
        gotoprompt = 0;
        % If the right arrow key is pressed
        if isequal(eventDat.Key,'rightarrow')
            if plotepoch <= ntrials  % make sure tStart is not too large
                plotepoch = plotepoch + 1;
            end
        % If the left arrow key is pressed
        elseif isequal(eventDat.Key,'leftarrow')
            if plotepoch >= 2  % tStart can't be < 1
                plotepoch = plotepoch - 1;
            end
        % If the "e" key is pressed
        elseif isequal(eventDat.Key,'e')
            prompt = 'Please input a vector of channels to plot: ';
            plotchans = input(prompt);
            gotoprompt = 1;
        end
        if gotoprompt,
            commandwindow; %Push focus to the command window
        end
        xtick = 1/datain.sr; 
        time = (epochtime*(plotepoch-1)+xtick):xtick:(epochtime*plotepoch);
        cortplotx(time,squeeze(datain.data(:,plotchans,plotepoch))+matspacing(:,1:length(plotchans)),'custom', chanlabels(plotchans), varargin{:});
        set(gca,'YTick',spacing(plotchans),'YTickLabel',chanlabels(plotchans));
        set(gca,'YLim',[spacing(plotchans(1))-spacer*2 spacing(plotchans(end))+spacer])
        xlabel('Time (sec)');
        % Plot horizontal green line for any markers
        if isfield(datain,'markers')
            markers = datain.markers;
            markers(markers == 0) = NaN;
            hold on;
            plot(time,squeeze(markers(:,plotchans,plotepoch))+matspacing(:,plotchans), 'g', 'linewidth', 2);
            hold off;
        end
    end %end of plotData

    function Exit_callback(~,~,~)
        % This exits the artifact screening
        exitnow=1;
        done=1;
    end % end of Exit_callback
end %end of the main function