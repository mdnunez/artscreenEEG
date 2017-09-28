function [xfreqs, outpower, fourier] = powerspec(datain,varargin)
%POWERSPEC - Calculates & plots power spectrum of sample*channel*trial EEG data
%
%Useage:  [xfreqs, outpower, fourier] = powerspec(datain,varargin);
%
%INPUTS:  datain.data - sample*channel*trial EEG data
%         datain.sr - sample rate
%
% OPTIONAL ARGUMENTS: (Passed in Name-Value Pairs)
%
%   baddata - a chan by trial matrix identifying data that should start
%             as rejected. If the .artifact field is also present, this
%             baddata matrix will take precedence. Def = []
%    freqs - lower and upper boundary of frequencies to plot,
%            default: [minimum_freq 50] (Hz)
%
%       dB - Units in standardized dB instead of standardized 
%            power, default: 0 (false)
%
%      noplot - Suppresses plot, default: 0 (false)
%
%      varargin - Any "plot" inputs after the first two,
%                 i.e., plot(X,Y,varargin);
%
%Outputs:  xfreqs - Frequencies plotted (x-axis)
%          outpower - Power values plotted (y-axis)
%          fourier - Fouier coefficients
%
%See also: EEGPLI

%% Copyright

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

%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  09/08/16       Michael Nunez              Adapted from 
%  09/16/17       Michael Nunez         Don't plot artifact trials
%  09/25/17       Michael Nunez         Fix for no .artifact field

%% Frequency interval

[plotvarargin,freqs,dB,noplot] = parsevar(varargin,'baddata',[],'freqs',[0 50],'dB',0,'noplot',0);

nchans=size(datain.data,2);
ntrials=size(datain.data,3);

% If a baddata matrix was entered as an argument
if ~isempty(baddata),
    datain.artifact=baddata; 
end

if ~isfield(datain,'artifact'),
    datain.artifact=zeros(nchans,ntrials);
end

%Recalculate if maximum larger than Nyquist frequency
nyquist = (2/5)*datain.sr;
if freqs(2) > nyquist
    fprintf('User defined maximum frequency %0.3f is larger than the Nyquist frequency %0.3f! \n',freqs(2),nyquist);
    fprintf('Using Nyquist frequnecy as maximum \n');
    freqs(2) = nyquist;
end
%Recalculate if minimum is smaller than Nyquist sampling rate
nsr = datain.sr/size(datain.data,1); %Nyquist sampling rate (Hz)
if freqs(1) < nsr && freqs(1) ~= 0
    fprintf('User defined minimum frequency %0.3f is smaller than Nyquist sampling rate %0.3f! \n',freqs(1),nsr);
    fprintf('Using Nyquist sampling rate as minimum frequency \n');
    freqs(1) = nsr;
end
if freqs(1) <= 1
    freqs(1) = 0;
    %freqs(1) = nsr;
end


%% Power
fourier = fft(datain.data)/size(datain.data,1);
plotfreqs = 0:nsr:freqs(2);
[~,minindex] = min(abs(freqs(1)-plotfreqs));
maxindex = length(plotfreqs);

%Find good trials
disp('Rejecting bad trials from .artifact field or baddata input matrix...');
goodtrials = ~all(datain.artifact,1);


power = mean(abs(fourier(:,:,goodtrials)).^2,3)*(2/(nsr)); %Power in standardized units (\muV^2/Hz)
if dB
    power = 10*log10(power);
end

% Set the power of bad channels to NaN
disp('Rejecting bad channels from .artifact field or baddata input matrix...');
power(:,all(datain.artifact,2)) = NaN;

if ~noplot
    cortplotx(plotfreqs(minindex:maxindex),power(minindex:maxindex,:),plotvarargin{:});
    xlabel('Frequency (Hz)');
    if dB
        ylabel('Standardized Log Power (10*log_{10}(\muV^2/Hz); dB)');
    else
        ylabel('Standardized Power (\muV^2/Hz)');
    end
    title('EEG Power Spectrum');
end

xfreqs = plotfreqs(minindex:maxindex);
outpower = power(minindex:maxindex,:);
fourier = fourier(minindex:maxindex,:,:);
