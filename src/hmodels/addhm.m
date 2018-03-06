% function datain = addhm(datain,hm)
%
% Loads the named headmodel and adds it to the data structure.
%
%   Currently available headmodels: 
%       acticap32 - 10-20 system 32 electrode cap
%       extnd1020 - Extended 10-20 system 64 electrode cap
%       egihc64 - newer 64 channel EGI hydrocell nets?
%       egihc128 - newer 128 channel EGI hydrocell nets
%       egihc256 - newer 256 channel EGI hydrocell nets
%       egihc256red - newer 256 channel EGI hydrocell minus useless face electrodes
%       eginn128 - newer 128 channel EGI hydrocell net with a few electrodes switched (specific usage for Human Neuroscience Lab)
%       antwave - ANT waveguard caps
%       ucsdeeg - eeg from the UCSD simultaneous EEG/MEG center
%       ucsdmeg - meg from the UCSD simultaneous EEG/MEG center
%       nscan128 - neuroscan 128 channel quikcaps like those in G107
%       biosemi - biosemi data from UCSD collaberators
%       
% Ex: segdata = addhm(segdata,'egi128');
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

function datain = addhm(datain,hm)
tmp=load([hm 'hm']);
datain.hm=tmp.(upper(hm));