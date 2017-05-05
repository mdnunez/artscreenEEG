function [newdata baselines] = eegBaseline(data,wind)
%EEGBASELINE - Shifts each trial of EEG data by the mean of data(wind,:,:) 
%              This function is useful for ERP studies
%
%Useage:  [newdata baselines] = eegBaseline(data,wind);
%
%Inputs:  data - sample*channel*trial or sample*trial*channel EEG data
%         wind - Window in which to calculate basline
%
%Outputs:  newdata = baselined data
%          baselines = Means used to shift data at each channel and trial
%
% Copyright (C) 2015 Michael D. Nunez, <mdnunez1@uci.edu>
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
%  10/21/15       Michael Nunez                Original code

%% Code

baselines = squeeze(mean(data(wind,:,:),1));
opbases = permute(repmat(baselines,[1 1 size(data,1)]),[3 1 2]);
newdata = data - opbases;