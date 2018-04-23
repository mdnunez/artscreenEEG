function [newdata lockindex badtrials] = epochsubset(data,newindex,varargin)
%EPOCHSUBSET - Reepochs each epoch of EEG data by timelocking each epoch
%              "i" to "newindex(i)" with maximum window size available
%
%Useage:  [newdata lockindex badtrials] = epochsubset(data,newindex,lockindex);
%    OR:  [newdata lockindex badtrials] = epochsubset(data,newindex);
%
%Inputs:  data - sample*channel*trial EEG data
%         newindex - Vector of length "trial"
%
%Optional inputs
%          lockindex - Sample in which newdata is timelocked
%                      Default: nanmin(newindex)       
%
%Outputs:  newdata - Re-timelocked EEG data
%          lockindex - Sample in which newdata is timelocked
%          badtrials - Index of trials where newindex contained nans
%

%% Copyright
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


%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  04/24/17       Michael Nunez                Original code
%  05/01/17       Michael Nunez             Return bad trial index
%  05/15/17       Michael Nunez             Lockindex optional input

%% Code

if ~isvector(newindex)
    error('newindex must be a vector!');
end

if size(data,3) ~= length(newindex)
    error('size(data,3) must be equal to length(newindex)');
end

if ~isempty(varargin) & isnumeric(varargin{3})
	lockindex = round(varargin{3});
else
	lockindex = round(nanmin(newindex));
end
fprintf('Using lockindex %d \n',lockindex);

windsize = (size(data,1) - nanmax(newindex)) + lockindex;

newdata = zeros([windsize, size(data,2), size(data,3)]);

for t=1:length(newindex)
    if ~isnan(newindex(t))
        begin_index = newindex(t)-lockindex + 1;
        end_index = windsize + begin_index -1;
        newdata(:,:,t) = data(begin_index:end_index,:,t);
    end
end

badtrials = find(isnan(newindex));