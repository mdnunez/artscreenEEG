% function dataout = artscreen(datain,varargin)
%
% PURPOSE:
% This function is used to screen raw segmented data for bad trials, bad
% channels, and other artifacts.  Removing this bad data enables a
% better decomposotion with ICA, as the quality of the decomposition is
% linked to the homogeneity of the data. Note that only unique artifacts (such
% as head movements, bad channels, and large amplifier artifacts) should be
% removed during this screening.  Artifacts that happen frequently, such as
% blinks, muscle tension, and EKG are easily corrected with ICA.
%
% This function uses 4 primary means of identifying artifacts.  First, it
% plots the mean variance of each channel on a topographic plot of the
% head.  Bad channels can be identified by having much higher or lower
% variance than their neighbors.  In a typical eyes-open data set, the
% central channels have the lowest variance, edge channels are somewhat
% higher, while frontal channels usually have high variance due to blinks.
% Note that since it is a mean being plotted, a channel with high variance
% may not be bad on every trial, but rather may just have a large artifact
% in one or more trials that increases its mean.
%
% The second way this function identifies artifacts is by looking at the
% variance of individual channels and trials, normalized by the average
% variance of each channel.  While this tool cannot identify channels that
% are consistently bad (e.g. a channel of pure noise will have similar
% variance across trials), it is very good at finding infrequent artifacts
% such as head movements and amplifier artfacts.  Note that the behavior
% of this normalized variance method is tied to the trial length.  Very
% short trials lead to noisy variance values, while very long trials will
% make artifacts less pronounced.  Ideal lengths seem to be in the 3-10
% second range.
%
% The third way this function can help identify artifacts is by plotting the
% non-normalized epoch variances.  This can be useful for finding high
% variance artifacts that are too consistent across epochs to show up in
% the normalized view.
%
% Finally, this function also contains a raw trial viewer, which can be used
% to find artifacts.  
%
%
% USAGE:
% Trials or channels can be rejected by entering their index in the
% appropriate box then clicking the 'Reject' button. By leaving the other
% box blank, the function assumes that you want the entire trial or entire
% channel rejected. Alternatively, you can reject specific channels on
% specific trials by populating both boxes and clicking 'Reject'.
%
% Note that you can enter multiple items in the rejection boxes. Entering
% '5 16' in the channel box and hiting reject will reject channels 5 and
% 16 for all trials. Entering '3 9' in trials and '10 12' in channels will
% only reject channels 10 and 12 on trials 3 and 9.
%
% A rejection threshold can also be set, which sets the max allowable
% normed variance ratio allowed. This acts in concert with the channel
% and trial rejection boxes, so you can set thresholds for specific
% channels and/or trials. When viewing the raw (not normalized) epoch
% variances, this threshold corresponds to the unnormalized variance
% values.
%
%
% ADDITIONAL CAPABILITIES:
% Flat channels (zero variance) will automatically be labeled as artifact.
%
% The data is average referenced at the end of the artifact screening
% process, plus you have the option of average referencing during the
% screening process if you find and large artifacts that have bled into
% other channels due to recording with a hardware average reference such
% as found in the ANT EEG system.  The value that has been subtracted from
% every channel is retained in the .avgref field of the data structure.
%
% By default, any channels that are rejected will be replaced with new data
% that is interpolated from neighboring channels using a spline function.
% This can be disabled using optional inputs.
%
%
% INPUT DATA STRUCTURE:
%   datain.data - sample by channel by trial data array
%
%   datain.sr - the sampling rate of the data
%
%   datain.hm - (recommended) headmodel for the EEG or MEG system following
%               Siyi Deng's conventions. Load antwavehm.mat for example.
%               This is needed for data interpolation and plotting channel
%               variance on the scalp.
%
%   datain.artifact - (optional) if the structure contains this field,
%                  then all of the trials and channels indicated in it
%                  will be rejected when you begin the review. This can
%                  be used to reject additional bad data that you may
%                  have found during the ICA review process, without
%                  starting the artifact screening from the beginning.
%
%
% OPTIONAL ARGUMENTS: (Passed in Name-Value Pairs)
%
%   baddata - a chan by trial matrix identifying data that should start
%             as rejected. If the .artifact field is also present, this
%             baddata matrix will take precedence. Def = []
%
%   badchans - a simpler alternative to passing a full baddata matrix where
%              you can just pass a list of bad channels. Will stack with
%              the .artifact field or the baddata input matrix (i.e. both
%              will be rejected). Def = []
%
%   badtrials - a simple way to pass a bad trial list.  As with badchans,
%               this will stack with any existing .artifact field or
%               baddata input matrix. Def = []
%
%   interpflag - a flag to interpolate rejected channels from the
%                neighboring electrodes after screening.  Def = 1
%
%   noreview - a flag to skip the review process and go straight to the
%              interpolation and average referencing.  This can be useful
%              when you want to go back and delete a few bad channels or
%              trials (using the badchans and badtrials arguments) and
%              rerun the interpolation. Def = 0;
%
%
% Part 1 of artscreenEEG's basic data cleaning functions:
%    artscreen.m => icasegdata.m => icareview.m => icatochan.m (optional)
%
% Open artscreen.m for full version history/changelog
%
% Copyright (C) 2013 Cort Horton, <chorton@uci.edu>
% Copyright (C) 2016 Michael D. Nunez, <mdnunez1@uci.edu>
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

function datain = artscreen(datain,varargin)

% VERSION HISTORY:
%   1.0 - created 6/26/12 by Cort Horton
%   1.1 - Added the automatic detection of an artifact field. This will
%         aid the ability to reject additional bad data that is found in
%         the ICA review process. 7/27/12
%   1.2 - Added the ability to undo previous rejection steps.  Also sped
%         up the screening loop slightly. Fixed the font sizes to deal with
%         a wider range of monitor sizes, and can also now automatically
%         adjust when resizing the figure. 11/13/12
%   1.3 - Changes obsolete due to later changes 4/8/13
%   1.4 - Added the ability to re-average reference the data during screening.
%         This is helpful for when large artifacts have leaked into other
%         channels due to average referencing (such as that which occurs
%         at the hardware level in the ANT EEG system). 6/22/13
%   1.5 - Also added a raw data viewer similar to that found in EEGLAB and
%         changed some defaults.  7/3/13
%   1.6 - Added the ability to view raw (unnormalized) epoch variances, as
%         well as the ability to set thresholds using these vars. 7/7/13
%   1.7 - Added the ability to pass a bad trial list and the ability to
%         disable the manual review. 7/12/13
%   1.8 - Changed naming and data structure conventions to aid in working
%         with simultaneous EEG and MEG data 9/20/13
%   1.9 - Checks for cell array data (i.e. unequal trial lengths) 2/25/14
%   2.0 - Incorporated a time series mode.  No longer calls trialviwer.m
%         Moved some buttons and fixed datacursor bugs.  6/12/14
%   2.1 - Added ability to pass in a cell array data field.  This is used
%         to screen data with variable trial lengths.  This functionality
%         was previously only avaialable by using the separate
%         artscreencell.m function.  6/13/14
%   2.2 - Added a scale in microvolts to lower right corner of timeseries
%         plots.  6/24/14
%   2.3 - Added hortontopo.m capability 8/23/16 - Michael Nunez
%   2.4 - Remove automatic removal of trials if 20% of channels are bad 
%         12/14/16 - Michael Nunez
%   2.5 - No reason to keep saving headmodels. 'hm' field can now be a string
%         12/19/16 - Michael Nunez

if nargin < 1; help artscreen; return; end;

% Checks if .data field is a cell array
if iscell(datain.data)
    cellmode=1;
else
    cellmode=0;
end

% Parse inputs;
[~,baddata,badchans,badtrials,interpflag,noreview]=...
    parsevar(varargin,'baddata',[],'badchans',[],'badtrials',[],...
    'interpflag',1,'noreview',0);

% Variables set by data
if cellmode
    nchans=size(datain.data{1},2);
    ntrials=length(datain.data);
    nsamps=size(datain.data{1},1);
    tlength=nsamps/datain.sr;
else
    nsamps=size(datain.data,1);
    nchans=size(datain.data,2);
    ntrials=size(datain.data,3);
    tlength=nsamps/datain.sr;
    
    % set xlabels based on trial length
    if tlength > 5; xscale=1; else xscale=1000; end;
    xlabs=(1:nsamps)*xscale/datain.sr;
end

% If a baddata matrix was entered as an argument
if ~isempty(baddata); datain.artifact=baddata; end;

% If the artifact field exists, start with that data rejected
if isfield(datain,'artifact')
    disp('Rejecting data from .artifact field or baddata input matrix...');
else
    datain.artifact=zeros(nchans,ntrials);
end

% Reject any channels in the badchans list
if ~isempty(badchans)
    disp('Rejecting channels from the badchans list...');
    datain.artifact(badchans,:)=1;
end

% Reject any trials in the badtrials list
if ~isempty(badtrials)
    disp('Rejecting trials from the badtrials list...');
    datain.artifact(:,badtrials)=1;
end

% If no channel positions, disable the topo and interpolation
if ~isfield(datain,'hm');
    disp('No headmodel found. Disabling topo plots and interpolation.');
    nohmodel = 1;
    strhmodel = 0;
    interpflag = 0;
elseif ischar(datain.hm)
    % If headmodel is a string, add hm
    hmname = datain.hm;
    datain = addhm(datain, hmname);
    strhmodel = 1;
    nohmodel = 0;
    interpflag = 1;
else
    strhmodel =0;
    nohmodel = 0;
    interpflag = 1;
end

% Calculate variance of every channel/trial
if cellmode
    thevars=zeros(nchans,ntrials);
    for k=1:ntrials;
        thevars(:,k)=var(datain.data{k});
    end
else
    thevars=squeeze(var(datain.data));
end

% Set spacing between channels on timeseries view
chansep=5*mean(mean(sqrt(thevars)));
datascale=1;

% Sets flat channels as artifact
datain.artifact(thevars==0)=1;

% Makes the GUI
screensize=get(0,'ScreenSize');
thefig=figure('pos',[50 50 screensize(3)-100 screensize(4)-150],'toolbar','figure',...
    'menubar','none','name','ArtScreen','DeleteFcn',@Exit_callback,'ResizeFcn',@SetFontsizes);
if noreview; set(thefig,'visible','off'); end;

% Deletes unneeded buttons to limit confusion
set(0,'Showhidden','on');
tmp=get(thefig,'Children');
barbuttons=get(tmp,'Children');
delete(barbuttons([1:6 8 9 12 13 15 16]));

% GUI Elements
mainplotaxis=axes('units','norm','pos',[1.6/30 1.5/20 19.8/30 17.4/20]);

if nohmodel
    topoaxis=axes('units','norm','pos',[22.5/30 12.6/20 6.5/30 6.35/20]);
else
    topoaxis=axes('units','norm','pos',[22.45/30 11.8/20 6.45/30 6.95/20]);
end

plottrialbutton=uicontrol('units','norm','pos',[22.5/30 10/20 3.5/30 1/20],...
    'style','pushbutton','string','Plot Trial','backgroundcolor','c',...
    'callback',@Plot_callback);

plottrialedit=uicontrol('units','norm','pos',[26.5/30 10/20 2/30 1/20],...
    'style','edit');

uicontrol('units','norm','pos',[22/30 3/20 7/30 6.5/20],...
    'style','frame');

triallabel=uicontrol('units','norm','pos',[22.5/30 7.8/20 2/30 1/20],...
    'style','text','HorizontalAlignment','left',...
    'string','Trials');

trialsedit=uicontrol('units','norm','pos',[24.5/30 8/20 4/30 1/20],...
    'style','edit');

chanlabel=uicontrol('units','norm','pos',[22.5/30 6.3/20 2/30 1/20],...
    'style','text','HorizontalAlignment','left',...
    'string','Chans');

chansedit=uicontrol('units','norm','pos',[24.5/30 6.5/20 4/30 1/20],...
    'style','edit');

threshlabel=uicontrol('units','norm','pos',[22.5/30 4.8/20 3/30 1/20],...
    'style','text','HorizontalAlignment','left',...
    'string','Threshold');

threshedit=uicontrol('units','norm','pos',[26/30 5/20 2.5/30 1/20],...
    'style','edit');

rejectbutton=uicontrol('units','norm','pos',[25.5/30 3.5/20 3/30 1/20],...
    'style','pushbutton','string','Reject','BackgroundColor','r',...
    'callback',@Reject_callback);

undobutton=uicontrol('units','norm','pos',[22.5/30 3.5/20 2.5/30 1/20],...
    'style','pushbutton','string','Undo','Enable','off',...
    'callback',@Undo_callback);

finishedbutton=uicontrol('units','norm','pos',[25.4/30 .9/20 3.7/30 1.45/20],...
    'style','pushbutton','string','Finished','BackgroundColor','g',...
    'callback',@Finished_callback);

avgrefbutton=uicontrol('units','norm','pos',[22.1/30 1.1/20 2.8/30 1.05/20],...
    'style','pushbutton','string','Avg Ref',...
    'callback',@AvgRef_callback);

normvarbutton=uicontrol('units','norm','pos',[1/30 19.1/20 2/30 .75/20],...
    'style','pushbutton','string','Norm Vars','Enable','off',...
    'callback',@goNorm_callback);

rawvarbutton=uicontrol('units','norm','pos',[3.25/30 19.1/20 2/30 .75/20],...
    'style','pushbutton','string','Raw Vars','Enable','on',...
    'callback',@goRaw_callback);

tseriesbutton=uicontrol('units','norm','pos',[5.5/30 19.1/20 2/30 .75/20],...
    'style','pushbutton','string','T. Series','Enable','on',...
    'callback',@goTseries_callback);

back10button=uicontrol('units','norm','pos',[.025 .025 .04 .04],...
    'style','pushbutton','string','<<','callback',@Back10Trials_callback,'visible','off');

back1button=uicontrol('units','norm','pos',[.075 .025 .04 .04],...
    'style','pushbutton','string','<','callback',@PrevTrial_callback,'visible','off');

trialnumbox=uicontrol('units','norm','pos',[.125 .025 .06 .04],...
    'style','edit','callback',@TrialNumber_callback,'visible','off','string','1');

forward1button=uicontrol('units','norm','pos',[.195 .025 .04 .04],...
    'style','pushbutton','string','>','callback',@NextTrial_callback,'visible','off');

forward10button=uicontrol('units','norm','pos',[.245 .025 .04 .04],...
    'style','pushbutton','string','>>','callback',@Forward10Trials_callback,'visible','off');

upscalebutton=uicontrol('units','norm','pos',[.64 .025 .03 .04],...
    'style','pushbutton','string','+','callback',@UpScale_callback,'visible','off');

downscalebutton=uicontrol('units','norm','pos',[.68 .025 .03 .04],...
    'style','pushbutton','string','-','callback',@DownScale_callback,'visible','off');

timeseriestext1=uicontrol('units','norm','pos',[.31 .0375 .15 .03],'style','text','visible','off');
timeseriestext2=uicontrol('units','norm','pos',[.31 .01 .15 .03],'style','text','visible','off');

theslider=uicontrol('units','norm','pos',[.475 .025 .15 .04],...
    'style','slider','Enable','off','callback',@Slider_callback,'visible','off');

if ~cellmode && tlength > 10;
    set(theslider,'Enable','on','min',0,'max',tlength-10,'value',0,'SliderStep',[1/(tlength-10) 1/((tlength-10)/10)]);
end

% Data screening loop
finished=0;
exitnow=0;
shownorm=1;
showraw=0;
showtimeseries=0;
interpdone=0;
trialtoplot=1;

SetFontsizes;

while 1;
    if noreview; Finished_callback; else done=0; end
    
    % Get the sum of all current rejections
    artifact=sum(datain.artifact,3);
    artifact(artifact>1)=1;
    
    % Sets channels as bad if more than 66% of the trials are bad
    mostlybadchans=sum(artifact,2)>(ntrials*.66);
    datain.artifact(mostlybadchans,:,end)=1;
    
    % % Sets trials as bad if more than 20% of the channels are bad
    % mostlybadtrials=sum(artifact)>(nchans*.2);
    % datain.artifact(:,mostlybadtrials,end)=1;
    
    % Get the sum of all current rejections
    artifact=sum(datain.artifact,3);
    artifact(artifact>1)=1;
    
    % Disable undo button if no rejections have been made
    if size(datain.artifact,3)==1 && finished==0; set(undobutton,'Enable','off'); end;
    
    % Find mean variances and maxes of good trials and channels
    goodchans=find(sum(artifact,2)<ntrials)';
    goodtrials=find(sum(artifact,1)<nchans);
    chanmeanvars=zeros(nchans,1);
    normvars=zeros(size(thevars));
    for c=goodchans;
        chanmeanvars(c)=mean(thevars(c,artifact(c,:)==0));
        normvars(c,:)=thevars(c,:)/chanmeanvars(c);
    end
    normvars(artifact==1)=0;
    
    if finished==0;
        % Plot the topography of channel variances or a bar graph if no
        % headmodel is present in the data structure
        if nohmodel
            bar(topoaxis,goodchans,sqrt(chanmeanvars(goodchans)),'b');
            xlim(topoaxis,[1 nchans]);
            ylim(topoaxis,[0 1.05]*max(sqrt(chanmeanvars)));
            xlabel(topoaxis,'Channel','Fontsize',round(basefontsize*.75));
        else
            topochanfont=round(basefontsize/2);
            if nchans > 140;
                topochanfont=topochanfont-2;
            elseif nchans > 80;
                topochanfont=topochanfont-1;
            end
            hortontopo(sqrt(chanmeanvars),datain.hm,...
                'goodchans',goodchans,'plotaxes',topoaxis,...
                'chanfontsize',topochanfont,'drawcolorbar',0,'drawoutliers',1);
        end
        title(topoaxis,'Mean Channel Variance','Fontsize',round(basefontsize*.8),'FontWeight','bold');
        
        % Plot the chosen data (normed vars, raw vars, or timeseries)
        if shownorm;
            plotnormvars=normvars;
            plotnormvars(normvars==0)=NaN;
            image(plotnormvars,'Parent',mainplotaxis,'CDataMapping','scaled');
            title(mainplotaxis,'Normalized Epoch Variance','Fontsize',round(basefontsize*.8),'FontWeight','bold');
            caxis(mainplotaxis,[0 5]);
            xlabel(mainplotaxis,'Trial','Fontsize',round(basefontsize*.75),'FontWeight','bold');
            ylabel(mainplotaxis,'Channel','Fontsize',round(basefontsize*.75),'FontWeight','bold');
            colorbar('peer',mainplotaxis,'Fontsize',round(basefontsize*.6));
            
        elseif showraw;
            plotnormvars=sqrt(thevars);
            plotnormvars(normvars==0)=NaN;
            image(plotnormvars,'Parent',mainplotaxis,'CDataMapping','scaled');
            title(mainplotaxis,'Raw Epoch Variance','Fontsize',round(basefontsize*.8),'FontWeight','bold');
            xlabel(mainplotaxis,'Trial','Fontsize',round(basefontsize*.75),'FontWeight','bold');
            ylabel(mainplotaxis,'Channel','Fontsize',round(basefontsize*.75),'FontWeight','bold');
            colorbar('peer',mainplotaxis,'Fontsize',round(basefontsize*.6));
            
        elseif showtimeseries
            if cellmode
                plotdata=datain.data{trialtoplot};
                nsamps=size(plotdata,1);
                tlength=nsamps/datain.sr;
                if tlength > 10;
                    set(theslider,'Enable','on');
                else
                    set(theslider,'Enable','off');
                end
                
                % set xlabels based on trial length
                if tlength > 5; xscale=1; else xscale=1000; end;
                xlabs=(1:nsamps)*xscale/datain.sr;
                
            else
                plotdata=datain.data(:,:,trialtoplot);
            end
            sumart=sum(datain.artifact,3);
            plotdata(:,sumart(:,trialtoplot)>0)=0;
            
            plotdata=datascale*plotdata-chansep*ones(nsamps,1)*(1:nchans);
            
            cortplotx(mainplotaxis,xlabs,plotdata);
            hold(mainplotaxis,'on');
            plot(mainplotaxis,[1 1]*prctile(xlabs,90),[-1*(nchans+3)*chansep -1*(nchans+3)*chansep+100*datascale],'k','linewidth',2);
            plot(mainplotaxis,[.9925 1.0075]*prctile(xlabs,90),[-1*(nchans+3)*chansep -1*(nchans+3)*chansep],'k','linewidth',2);
            plot(mainplotaxis,[.9925 1.0075]*prctile(xlabs,90),[-1*(nchans+3)*chansep+100*datascale -1*(nchans+3)*chansep+100*datascale],'k','linewidth',2);
            text(prctile(xlabs,91),-1*(nchans+2.5)*chansep,'100\muV','fontsize',16,'parent',mainplotaxis);
            hold(mainplotaxis,'off');
            ylim(mainplotaxis,[-1*(nchans+5)*chansep chansep*5]);
            title(mainplotaxis,'Trial Time Series','Fontsize',round(basefontsize*.8),'FontWeight','bold');
            if tlength>10;
                xlim(mainplotaxis,[0 10]);
            else
                xlim(mainplotaxis,[0 max(xlabs)]);
            end
            set(mainplotaxis,'YTickLabel',[]);
            if xscale > 1;
                xlabel(mainplotaxis,'Time (ms)','fontsize',round(basefontsize*.7));
            else
                xlabel(mainplotaxis,'Time (s)','fontsize',round(basefontsize*.7));
            end
        end

        % setup modified data cursor
        if shownorm || showraw;
            dch=datacursormode(thefig);
            datacursormode on;
            set(dch,'UpdateFcn',@dc_updatefcn);
        else
            datacursormode off;
        end
    end
    
    % wait for user action
    while done==0;
        pause(.05);
    end
    
    if finished
        break
    elseif exitnow
        return
    end
end

% Collapse all rejections
datain.artifact=sum(datain.artifact,3);
datain.artifact(datain.artifact>1)=1;

% Zero out all rejected data
if cellmode
    for t=1:ntrials;
        datain.data{t}(:,datain.artifact(:,t)==1)=0;
    end
else
    datain.data(:,datain.artifact==1)=0;
end

% Data Interpolation
if interpflag
    disp('Interpolating rejected channels...');
    chanpos=datain.hm.Electrode.CoordOnSphere;
    interpchans=setdiff(1:size(chanpos,1),datain.hm.Electrode.NoInterp);
    
    for t=goodtrials
        mat=splineinterp(.1,chanpos(interpchans(datain.artifact(interpchans,t)==0),:),chanpos(interpchans,:));
        if cellmode
            datain.data{t}(:,interpchans)=datain.data{t}(:,interpchans(datain.artifact(interpchans,t)==0))*mat';
        else
            datain.data(:,interpchans,t)=datain.data(:,interpchans(datain.artifact(interpchans,t)==0),t)*mat';
        end
    end
    interpdone=1;
end

% If headmodel was a string, convert 'hm' field back to a string
if strhmodel
    datain.hm = hmname;
end

% Average reference final output
doAvgRef;


%% Nested Functions: These share a workspace with the main function********
% *************************************************************************
% *************************************************************************
    function doAvgRef
        % This function average references the data. Artifact channels are
        % not included in the average.
        
        disp('Average referencing the data...');
        
        if ~interpdone
            % Zero out all rejected data first
            datain.artifact=sum(datain.artifact,3);
            datain.artifact(datain.artifact>1)=1;
            
            if cellmode
                for t=1:ntrials;
                    datain.data{t}(:,datain.artifact(:,t)==1)=0;
                end
            else
                datain.data(:,datain.artifact==1)=0;
            end
        end
        
        % Concatenate all data
        [catdata,tsamples]=segtocat(datain.data);
        
        % Find flat chans/trials
        zerodata=catdata==0;
        
        % Calculate average reference of good data
        catref=sum(catdata,2)./sum(catdata~=0,2);
        catdata=catdata-catref*ones(1,nchans);
        catdata(zerodata)=0;
        
        if cellmode
            for t=1:ntrials;
                datain.data{t}=catdata(sum(tsamples(1:(t-1)))+1:sum(tsamples(1:t)),:);
            end
        else
            datain.data=cattoseg(catdata,nsamps);
        end
        
    end % end of doAvgRef

    function AvgRef_callback(~,~,~)
        % Performs average referencing of the data
        doAvgRef;
        if cellmode
            thevars=zeros(nchans,ntrials);
            for t=1:ntrials;
                thevars(:,t)=var(datain.data{t});
            end
        else
            thevars=squeeze(var(datain.data));
        end
        done=1;
    end % end of AvgRef_callback

    function goNorm_callback(~,~,~)
        % Changes view to normalized epoch variances
        shownorm=1;
        showraw=0;
        showtimeseries=0;
        set(normvarbutton,'Enable','off');
        set(rawvarbutton,'Enable','on');
        set(tseriesbutton,'Enable','on');
        SetTimeseriesInvisible;
        done=1;
    end % end of goNorm_callback

    function goRaw_callback(~,~,~)
        % Changes view to raw epoch variances
        shownorm=0;
        showraw=1;
        showtimeseries=0;
        set(normvarbutton,'Enable','on');
        set(rawvarbutton,'Enable','off');
        set(tseriesbutton,'Enable','on');
        SetTimeseriesInvisible;
        done=1;
    end % end of goRaw_callback

    function goTseries_callback(~,~,~)
        % Changes view to trial timeseries
        shownorm=0;
        showraw=0;
        showtimeseries=1;
        set(normvarbutton,'Enable','on');
        set(rawvarbutton,'Enable','on');
        set(tseriesbutton,'Enable','off');
        SetTimeseriesVisible;
        SetFontsizes;
        done=1;
    end % end of goRaw_callback

    function Plot_callback(~,~,~)
        % This Callback function is used to plot individual trials. The spacing
        % between channels is determined by the data.
        
        newtrial=str2double(get(plottrialedit,'string'));
        if ismember(newtrial,1:ntrials);
            trialtoplot=newtrial;
        else
            set(plottrialedit,'string',num2str(trialtoplot));
        end
        
        showtimeseries=1;
        showraw=0;
        shownorm=0;
        
        set(normvarbutton,'Enable','on');
        set(rawvarbutton,'Enable','on');
        set(tseriesbutton,'Enable','off');
        
        set(plottrialedit,'string',[]);
        set(trialnumbox,'string',num2str(trialtoplot));
        
        SetTimeseriesVisible;
        SetFontsizes;
        
        done=1;
    end % end of Plot_callback

    function Reject_callback(~,~,~)
        % This callback function is used to reject data and advance the loop to the
        % next iteration.  It also clears the edit boxes. If channel or trial are
        % left blank, it assumes you want all trials or channels.
        
        rejtrials=str2num(get(trialsedit,'string'));
        rejchans=str2num(get(chansedit,'string'));
        rejthresh=str2num(get(threshedit,'string'));
        
        if isempty(rejtrials) && isempty(rejchans) && isempty(rejthresh)
            done=1;
            return;
        end
        
        if isempty(rejtrials); rejtrials=1:ntrials; end;
        if isempty(rejchans); rejchans=1:nchans; end;
        if isempty(rejthresh); rejthresh=0; end;
        
        % Checks for bad inputs
        if sum(ismember(rejtrials,1:ntrials))<length(rejtrials);
            set(trialsedit,'string',[]);
            done=1;
            return;
        elseif sum(ismember(rejchans,1:nchans))<length(rejchans);
            set(chansedit,'string',[]);
            done=1;
            return;
        elseif rejthresh < 0;
            set(threshedit,'string',[]);
            done=1;
            return;
        end
        
        datain.artifact=cat(3,datain.artifact,zeros(nchans,ntrials));
        
        if shownorm;
            for trial=rejtrials;
                for chan=rejchans;
                    if normvars(chan,trial)>rejthresh;
                        datain.artifact(chan,trial,end)=1;
                    end
                end
            end
        else
            for trial=rejtrials;
                for chan=rejchans;
                    if sqrt(thevars(chan,trial))>rejthresh;
                        datain.artifact(chan,trial,end)=1;
                    end
                end
            end
        end
        
        set(trialsedit,'string',[]);
        set(chansedit,'string',[]);
        set(threshedit,'string',[]);
        set(undobutton,'Enable','on');
        done=1;
    end % end of Reject_callback

    function Finished_callback(~,~,~)
        % This callback function breaks the data review loop and closes the figure
        set(thefig,'DeleteFcn','');
        close(thefig);
        finished=1;
        done=1;
        drawnow;
    end % end of Finished_callback

    function Exit_callback(~,~,~)
        % This exits the artifact screening and returns an empty structure
        datain=[];
        exitnow=1;
        done=1;
    end % end of Exit_callback

    function Undo_callback(~,~,~)
        % This undoes the last rejection choice
        datain.artifact=datain.artifact(:,:,1:end-1);
        done=1;
    end % end of Undo_callback

    function txt = dc_updatefcn(~,event_obj)
        pos=get(event_obj,'Position');
        whichplot=get(event_obj,'Target');
        if get(whichplot,'parent')==mainplotaxis && shownorm;
            txt={['Trial: ' num2str(pos(1))],['Chan: ' num2str(pos(2))],...
                ['N.Var: ' num2str(plotnormvars(pos(2),pos(1)),3)]};
        elseif get(whichplot,'parent')==mainplotaxis && showraw;
            txt={['Trial: ' num2str(pos(1))],['Chan: ' num2str(pos(2))],...
                ['Var: ' num2str(plotnormvars(pos(2),pos(1)),3)]};    
        elseif get(whichplot,'parent')==topoaxis && nohmodel;
            txt={['Chan: ' num2str(pos(1))],['Var: ' num2str(pos(2),3)]};
        else
            txt=[];
        end
    end % end of dc_updatefcn

    function SetFontsizes(~,~,~)
        % Sets the font sizes for buttons and titles. Allows for better scaling
        % for various different-sized monitors
        guipos=get(thefig,'pos');
        basefontsize=round(guipos(3)/70);
        if guipos(3)/guipos(4) > 3;
            basefontsize=ceil(basefontsize*.6);
        elseif guipos(3)/guipos(4) > 2;
            basefontsize=ceil(basefontsize*.75);
        end
        set(mainplotaxis,'fontsize',round(basefontsize*.6));
        set(topoaxis,'fontsize',round(basefontsize*.6));
        set(plottrialedit,'fontsize',basefontsize);
        set(trialsedit,'fontsize',basefontsize);
        set(chansedit,'fontsize',basefontsize);
        set(threshedit,'fontsize',basefontsize);
        set(plottrialbutton,'fontsize',basefontsize);
        set(triallabel,'fontsize',basefontsize);
        set(chanlabel,'fontsize',basefontsize);
        set(threshlabel,'fontsize',basefontsize);
        set(rejectbutton,'fontsize',basefontsize);
        set(finishedbutton,'fontsize',round(basefontsize*1.25));
        set(undobutton,'fontsize',basefontsize);
        set(avgrefbutton,'fontsize',basefontsize);
        set(normvarbutton,'fontsize',round(basefontsize*.6));
        set(rawvarbutton,'fontsize',round(basefontsize*.6));
        set(tseriesbutton,'fontsize',round(basefontsize*.6));
        
        if showtimeseries
            set(back10button,'fontsize',round(basefontsize*.6));
            set(back1button,'fontsize',round(basefontsize*.6));
            set(trialnumbox,'fontsize',round(basefontsize*.6));
            set(forward1button,'fontsize',round(basefontsize*.6));
            set(forward10button,'fontsize',round(basefontsize*.6));
            set(upscalebutton,'fontsize',round(basefontsize*.6));
            set(downscalebutton,'fontsize',round(basefontsize*.6));
            
            tmp=get(thefig,'Color');
            set(timeseriestext1,'fontsize',basefontsize*.5,'BackgroundColor',tmp,'string',['Channels: ' num2str(nchans) '   Trials: ' num2str(ntrials)]);
            set(timeseriestext2,'fontsize',basefontsize*.5,'BackgroundColor',tmp,'string',['Trial Length: ' num2str(tlength) 's   SR: ' num2str(datain.sr) 'Hz']);
        end
        done=1;
    end % end of SetFontsizes

    function SetTimeseriesInvisible(~,~,~)
        set(mainplotaxis,'units','norm','pos',[1.6/30 1.5/20 18.72/30 17.4/20]);
        set(back10button,'visible','off');
        set(back1button,'visible','off');
        set(trialnumbox,'visible','off');
        set(forward1button,'visible','off');
        set(forward10button,'visible','off');
        set(upscalebutton,'visible','off');
        set(downscalebutton,'visible','off');
        set(timeseriestext1,'visible','off');
        set(timeseriestext2,'visible','off');
        set(theslider,'visible','off');
    end % end of SetTimeseriesInvisible

    function SetTimeseriesVisible(~,~,~)
        set(mainplotaxis,'units','norm','pos',[1.6/30 2.9/20 19.8/30 16.0/20]);
        set(back10button,'visible','on');
        set(back1button,'visible','on');
        set(trialnumbox,'visible','on');
        set(forward1button,'visible','on');
        set(forward10button,'visible','on');
        set(upscalebutton,'visible','on');
        set(downscalebutton,'visible','on');
        set(timeseriestext1,'visible','on');
        set(timeseriestext2,'visible','on');
        set(theslider,'visible','on');
    end % end of SetTimeseriesVisible

    function Back10Trials_callback(~,~)
        trialtoplot=max([1 trialtoplot-10]);
        set(trialnumbox,'string',num2str(trialtoplot));
        done=1;
    end % end of nested function

    function PrevTrial_callback(~,~)
        trialtoplot=max([1 trialtoplot-1]);
        set(trialnumbox,'string',num2str(trialtoplot));
        done=1;
    end % end of nested function

    function TrialNumber_callback(~,~)
        newtrial=str2double(get(trialnumbox,'string'));
        if ismember(newtrial,1:ntrials);
            trialtoplot=newtrial;
        else
            set(trialnumbox,'string',num2str(trialtoplot));
        end
        done=1;
    end % end of nested function

    function NextTrial_callback(~,~)
        trialtoplot=min([trialtoplot+1 ntrials]);
        set(trialnumbox,'string',num2str(trialtoplot));
        done=1;
    end % end of nested function

    function Forward10Trials_callback(~,~)
        trialtoplot=min([ntrials trialtoplot+10]);
        set(trialnumbox,'string',num2str(trialtoplot));
        done=1;
    end % end of nested function

    function UpScale_callback(~,~)
        datascale=datascale/.75;
        done=1;
    end % end of nested function

    function DownScale_callback(~,~)
        datascale=datascale*.75;
        done=1;
    end % end of nested function

    function Slider_callback(hobj,~)
        tmp=get(hobj,'value');
        xlim(mainplotaxis,[tmp tmp+10]);
    end % end of nested function

% End of nested functions *******************************************************
%********************************************************************************
%********************************************************************************

end % end of main function
