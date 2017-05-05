% function dataout = lapdata(datain,varargin)
%
% This function performs a spherical Laplacian
%
% Required Input:
%    Data structure containing the fields .data and .sr 
%
% Optional Arguments: 
%
%    badchans - channels that are not zeroed out but should be considered
%               bad and excluded from the ICA.
%                    
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

function datain = lapdata(datain,varargin)

% VERSION HISTORY:
%   1.0 - created by Michael Nunez  4/24/17

% To do:
% 1) Display spherical Laplacian?
% 2) Record the outer ring in head models to be ignored after the spherical-spline-Laplacian filter


if nargin < 1; help lapdata; return; end;

% Parse inputs;
[~,badchans]=...
    parsevar(varargin,'badchans',[]);


% Performs a Laplacian filter on the data
fprintf('Calculating a spherical spline Laplacian filter...\n');
chanpos=datain.hm.Electrode.CoordOnSphere;
%Need to remove the outer ring for Laplacian channels
if isfield(datain.hm.Electrode,'NoInterp')
    lapchans=setdiff(1:size(chanpos,1),datain.hm.Electrode.NoInterp);
    if isfield(datain.hm.Electrode,'NoLap')
        lapchans=setdiff(lapchans,datain.hm.Electrode.NoLap);
    end
else
    lapchans = 1:size(chanpos,1);
end


for t=goodtrials
    mat=splinenlap(.1,chanpos(lapchans(datain.artifact(lapchans,t)==0),:),chanpos(lapchans,:));
    if cellmode
        datain.data{t}(:,lapchans)=datain.data{t}(:,lapchans(datain.artifact(lapchans,t)==0))*mat';
    else
        datain.data(:,lapchans,t)=datain.data(:,lapchans(datain.artifact(lapchans,t)==0),t)*mat';
    end
end