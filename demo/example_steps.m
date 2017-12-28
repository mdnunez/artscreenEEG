%EXAMPLE_STEPS - Script that follows the steps to clean data
%
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
%  8/23/16        Michael Nunez                 Original code
%  8/25/16        Michael Nunez       figshare link, updated feedback
%  12/1/16        Michael Nunez     Download warning only if .mat is non-local
%  1/17/17        Michael Nunez  Changing name of structure to avoid confusion
%  8/15/17        Michael Nunez         Use of icaadjust().
%  12/27/17       Michael Nunez  Open icareview() after icaadjust(), fix last save
%  12/28/17       Michael Nunez         Laplacian

%% Initial
sub = 'subject1';

%% Code
%Download data from figshare.com
if ~(exist('subject1.mat') == 2)
	%% See warning to advance
	warning('The following script will download and create data in the current working directory!');
	fprintf('Press any key to advance (CTRL+C to exit)...\n');
	pause;
	fprintf('%%Downloading example data (~240MB) in the current working directory from figshare.com...\n');
	urlwrite(...
		'https://ndownloader.figshare.com/files/5865048?private_link=e2405744c7a1346ea17b',...
		'subject1.mat');
	fprintf('%%Loading downloaded data into MATLAB workspace...\n');
else
	fprintf('%%subject1.mat file found locally. Loading this data...\n');
end

fprintf('\n');
fprintf('%%Note that all of the base artscreenEEG functions take MATLAB structure inputs with fields ''data'' and ''sr'' at least.\n');
fprintf('%%Field ''hm'' is a recommended head model field. ''hm'' itself is a structure with fields ''CoordOnSphere'', ''SphereExtendAngle'', and ''Scaling2D''.\n');
fprintf('\n');

%The following data has been bandpass filtered from 1 to 100 Hz with a notch at 60 Hz
%Some channels have been premarked as artifact
eeg = load(sprintf('%s.mat',sub)); %Load all variables in structure 'eeg'
fprintf('%%Adding head model...\n');
fprintf('%%Alternate spherical head model .mat files are contained in src/hmodels.\n');
eeg = addhm(eeg,'eginn128');

fprintf('%%Loading downloaded data into artscreen()...\n');
fprintf('%%Try a 300 abs variance cutoff...\n');
%Use 300 abs variance cutoff, reject poor channels not in parietal and occipital locations
eeg = artscreen(eeg); 

%Independent component analysis on data
fprintf('%%Running ICA on artscreened data...\n');
ica = icasegdata(eeg,'ncomps',60,'nkeep',60,'fftfreq',100);

%Save ICA data
fprintf('%%Saving ICA data...\n');
save(sprintf('%s_ica.mat',sub),'-struct','ica');

% % Manual review of ICA components
% fprintf('%%Review ICA components...\n');
% fprintf('%%Suggested component evaluations are included. Note that IC evaluation is inherently subjective...\n');
% ica.compevals = [0,0,0,2,2,2,0,2,2,0,0,2,1,2,2,0,0,0,1,1,0,0,0,2,0,0,0,0,1,zeros(1,31)];
% ica = icareview(ica);

% Automatic review of ICA components
fprintf('%%Automatic review of ICA components...\n');
ica = icaadjust(ica);
fprintf('%%Review ICA components...\n');
ica = icareview(ica);

%Save ICA data with component evaluations
fprintf('%%Overwriting ICA data with evaluated component information...\n');
save(sprintf('%s_ica.mat',sub),'-struct','ica');

%ICs to Chan data: Keep components which are marked as "unsure"
fprintf('%%Keeping all components marked as "good" or "unsure"...\n');
eeg = icatochan(ica,1); 

%Filter data
% Lowpass the data
Fpass = 50;          % Passband Frequency
Fstop = 60;          % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 10;          % Stopband Attenuation (dB)
match = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, 1000);
Hd = design(h, 'butter', 'MatchExactly', match);

% Carry out the filtering
fprintf('%%Lowpass filtering the data at 50 Hz...\n');
eeg.data=filtfilthd(Hd,eeg.data);

%Remove any probable bad trials observed during "icareview"
fprintf('%%A second pass of artscreen() is recommended...\n');
eeg = artscreen(eeg);

%Calculate Spherical Laplacian
fprintf('%%Calculate spherical Laplacian using sphlapdata()...\n');
eeg = sphlapdata(eeg);

%Save cleaned data
fprintf('%%Saving the cleaned data...\n');
save(sprintf('%s_cleaned.mat',sub),'-struct','eeg');

