function [goodchans,goodtrials,badchans,badtrials] = goodbad(datain)
%GOODBAD - Finds suggested goodchans and goodtrials from datain.artifact
%
%Useage:  [goodchans,goodtrials,badchans,badtrials] = goodbad(datain);
%
%Inputs:  datain - data structure with 'artifact' nchan*ntrial field
%         
%
%Outputs:  goodchans - good channel index
%          goodtrials - good trial index
%          badchans - bad channel index
%          badtrials - bad trial index
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
%  05/15/17       Michael Nunez                Original code

%% Code

if ~isfield(datain,'artifact')
    help goodbad;
end

ntrials = size(datain.artifact,2);
nchans = size(datain.artifact,1);

badchans = find(sum(datain.artifact,2) == ntrials);
badtrials = find(sum(datain.artifact,1) == nchans);
goodchans = setdiff(1:nchans,badchans);
goodtrials = setdiff(1:ntrials,badtrials);