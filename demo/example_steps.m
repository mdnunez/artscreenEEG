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

%% Initial
sub = 'subject1';

%% See warning to advance
warning('The following script will download and create data in the current working directory!');
fprintf('Press any key to advance (CTRL+C to exit)...\n');
pause;

%% Code
%Download data from figshare.com
if ~(exist('subject1.mat') == 2)
	fprintf('%%Downloading example data (~240MB) in the current working directory from figshare.com...\n');
	urlwrite(...
		'https://ndownloader.figshare.com/files/5865048?private_link=e2405744c7a1346ea17b',...
		'subject1.mat');
	fprintf('%%Loading downloaded data into MATLAB workspace...\n');
else
	fprintf('%%subject1.mat file found locally. Loading this data...\n');
end

%The following data has been bandpass filtered from 1 to 100 Hz with a notch at 60 Hz
%Some channels have been premarked as artifact
data = load(sprintf('%s.mat',sub)); %Can this be loaded from figshare.com?
fprintf('%%Adding headmodel...\n');
data = addhm(data,'eginn128');

fprintf('%%Loading downloaded data into artscreen()...\n');
fprintf('%%Try a 300 abs variance cutoff...\n');
%Use 300 abs variance cutoff, reject poor channels not in parietal and occipital locations
data = artscreen(data); 

%Independent component analysis on data
fprintf('%%Running ICA on artscreened data...\n');
ica = icasegdata(data,'ncomps',60,'nkeep',60,'fftfreq',100);

%Save ICA data
fprintf('%%Saving ICA data...\n');
save(sprintf('%s_ica.mat',sub),'-struct','ica');

%Review ICA components
fprintf('%%Review ICA components...\n');
fprintf('%%Suggested component evaluations are included. Note that IC evaluation is inherently subjective...\n');
ica.compevals = [0,0,0,2,2,2,0,2,2,0,0,2,1,2,2,0,0,0,1,1,0,0,0,2,0,0,0,0,1,zeros(1,31)];
ica = icareview(ica);

%Save ICA data with component evaluations
fprintf('%%Overwriting ICA data with evaluated component information...\n');
save(sprintf('%s_ica.mat',sub),'-struct','ica');

%ICs to Chan data: Keep components which are marked as "unsure"
fprintf('%%Keeping all components marked as "good" or "unsure"...\n');
data = icatochan(ica,1); 

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
data.data=filtfilthd(Hd,data.data);

%Remove any probable bad trials observed during "icareview"
fprintf('%%A second pass of artscreen() is recommended...\n');
data = artscreen(data);

%Save cleaned data
fprintf('%%Saving the cleaned data...\n');
save(sprintf('%s_cleaned.mat',sub),'-struct','data');

