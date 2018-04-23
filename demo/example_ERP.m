%EXAMPLE_ERP - Script that calculates a simple visual evoked potential
%              or ERP (event-related potential)

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
%  04/23/18       Michael Nunez                 Original code

%% Initial
sub = 'subject1';

%% Code
%Download data from figshare.com
if ~(exist('subject1_cleaned.mat') == 2)
	error('Please run example_steps.m first!')
end

fprintf('\n');
fprintf('%%Note that the stimulus onset occurred at sample 1251 in each trial\n');
fprintf('\n');

%The following data has been bandpass filtered from 1 to 100 Hz with a notch at 60 Hz
%Some channels have been premarked as artifact

%Filter data from 1 to 5 Hz (data has previously been filtered from 1 to 50 Hz)
% Lowpass the data
Fpass = 5;          % Passband Frequency
Fstop = 10;          % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 10;          % Stopband Attenuation (dB)
match = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, eeg.sr);
Hd = design(h, 'butter', 'MatchExactly', match);

% Carry out the filtering
fprintf('%%Lowpass filtering the data at 5 Hz...\n');
filtered =filtfilthd(Hd,eeg.data);

% Baseline the data
fprintf('%%Baseline the data 100 ms before stimulus onset...\n');
[baselined, baselines] = eegBaseline(filtered, 1151:1250);

% Calculate the ERP 
fprintf('%%Calculate the ERP by averaging ''good'' baselined trials (for only ''good'' channels)\n');
[goodchans, goodtrials] = goodbad(eeg);
erp = mean(baselined(1151:1350, goodchans, goodtrials));

fprintf('%%Plot the visual ERP\n');
figure;
plotx(-99:1000,erp);
ylims = get(gca,'YLim');
lhandle = line([1 1],get(gca,'YLim'));
set(gca,'YLim',ylims);
set(lhandle,'LineStyle','--','Color','r','LineWidth',2);
xlabel('Time after stimulus onset (ms)');
ylabel('Amplitude (Microvolts)');
title('ERP for each non-artifact electrode');



