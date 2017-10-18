function data=filtereeg(data,sr,passband,stopband,stopdB)
% function fdatastruct=filtereeg(datastruct,[],[passband],[stopband],[stopdB])
% function fdata=filtereeg(dataarray,sr,[passband],[stopband],[stopdB])
%
% A basic set of Butterworth IIR filters for EEG data that will not affect 
% phase. By default is uses a passband from 1 to 50 Hz with 1 dB ripple and 
% a stopband of .25 and 60 Hz with 10dB attenuation. Additionally employs a 
% 60 Hz notch filter to remove line noise.
%
% Input can be either a structure with fields .sr and (.eeg or .data), or
% you can pass the data and sampling rate as two arguments.  If you give it
% a structure, it returns a structure.  If you give two arguments, it
% returns just the filtered data array.
%
% Required Inputs:
%   data - structure with fields .data and .sr 
%
%            ****** OR ****** 
%
%   data - time by chan (by trial) data array
%   sr - sampling rate in Hz
%  
% Optional Inputs:
%   passband - 2 element vector indicating passband. Def [1 50]
%   stopband - 2 element vector indicating cutoff freqs. Def [.25 60]
%   stopdB - amount of attenuation at the stopband freqs. Can be either a 
%            single number or a 3 element vector specifying attenuation at 
%            the low cutoff, high cutoff, and notch. Def 10
%
%
% Copyright (C) 2013 Cort Horton, <chorton@uci.edu>
% Copyright (C) 2017 Michael D. Nunez, <mdnunez1@uci.edu>
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
%
% To do:
% Use: parsevar()
if nargin < 5; stopdB = 10; end;
if nargin < 3 || isempty(passband) || isempty(stopband); 
    passband=[1 50];
    stopband=[.25 60];
end
if length(stopdB)==1; stopdB=ones(1,3) * stopdB; end;
stopdB=stopdB/2; % filtfilt doubles attenuation, so we divide by two

% detects if input is a structure
convertback = 0;
if isstruct(data)
    convertback = 1;
    tmp=data;
    sr=data.sr;
    if isfield(data,'eeg')
        data=data.eeg;
        datafieldname='eeg';
    elseif isfield(data,'data')
        data=data.data;
        datafieldname='data';
    else
        error('Unrecognized data field');
    end
end

% Detrend the data
data=ndetrend(data);

% Highpass the data
Fpass = passband(1);  % Passband Frequency
Fstop = stopband(1);  % Stopband Frequency
Apass = 1;            % Passband Ripple (dB)
Astop = stopdB(1);    % Stopband Attenuation (dB)
match = 'passband';   % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, sr);
Hd = design(h, 'butter', 'MatchExactly', match);

% Carry out the filtering
data=filtfilthd(Hd,data);

% Notch at 60 Hz
Fpass1 = 59;          % First Passband Frequency
Fstop1 = 59.9;        % First Stopband Frequency
Fstop2 = 60.1;        % Second Stopband Frequency
Fpass2 = 61;          % Second Passband Frequency
Apass1 = 1;           % First Passband Ripple (dB)
Astop  = stopdB(3);   % Stopband Attenuation (dB)
Apass2 = 1;           % Second Passband Ripple (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                      Apass2, sr);
Hd = design(h, 'butter', 'MatchExactly', match);

% Carry out the filtering
data=filtfilthd(Hd,data);

% Lowpass the data
Fpass = passband(2);  % Passband Frequency
Fstop = stopband(2);  % Stopband Frequency
Apass = 1;            % Passband Ripple (dB)
Astop = stopdB(2);    % Stopband Attenuation (dB)
match = 'passband';   % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, sr);
Hd = design(h, 'butter', 'MatchExactly', match);

% Carry out the filtering
data=filtfilthd(Hd,data);

% if structure input, return output to structure form
if convertback;
    tmp.(datafieldname)=data;
    data=tmp;
end
