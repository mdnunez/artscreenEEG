% function dataout = sphlapdata(datain,varargin)
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

function datain = sphlapdata(datain,varargin)

% VERSION HISTORY:
%   1.0 - created by Michael Nunez  5/05/17

% To do:
% 1) Apply spherical Laplacian on good channels for each trial or overall?
% 2) Display spherical Laplacian?
% 3) Record the outer ring in head models to be ignored after the spherical-spline-Laplacian filter


if nargin < 1; help sphlapdata; return; end;

% Parse inputs;
[~,badchans]=...
    parsevar(varargin,'badchans',[]);

% Determined from the data
nsamps=size(datain.data,1);
nchans=size(datain.data,2);
ntrials=size(datain.data,3);
tlength=nsamps/datain.sr;

disp('Assuming all non-zero and non-NaN data are good...');
thevars=squeeze(var(datain.data));
artifact=thevars==0 | isnan(thevars);

if isfield(datain,'artifact')
     disp('Rejecting data from .artifact field matrix...');
     artifact = datain.artifact | artifact;
end

goodtrials=setdiff(1:ntrials,find(sum(artifact)==nchans));
goodchans=setdiff(1:nchans,find(sum(artifact,2)==ntrials));
goodchans=setdiff(goodchans,badchans);
ngoodchans=length(goodchans);
ngoodtrials=length(goodtrials);

% Checks if .data field is a cell array
if iscell(datain.data)
    cellmode=1;
else
    cellmode=0;
end

% Performs a Laplacian filter on the data
fprintf('Calculating a spherical spline Laplacian filter...\n');
chanpos=datain.hm.Electrode.CoordOnSphere;
%Need to remove the outer ring for Laplacian channels
if isfield(datain.hm.Electrode,'NoInterp')
    poschans=setdiff(1:size(chanpos,1),datain.hm.Electrode.NoInterp);
    lapchans=setdiff(poschans,badchans);
    if isfield(datain.hm.Electrode,'NoLap')
        lapchans=setdiff(lapchans,datain.hm.Electrode.NoLap);
    end
else
    poschans = 1:size(chanpos,1);
    lapchans=setdiff(poschans,badchans);
end

datain.sphlap = zeros(size(datain.data));
for t=goodtrials
    mat=splinenlap(.1,chanpos(lapchans(datain.artifact(lapchans,t)==0),:),chanpos(lapchans,:)); %Check this behavior
    if cellmode
        datain.sphlap{t}(:,lapchans)=datain.data{t}(:,lapchans(datain.artifact(lapchans,t)==0))*mat';
    else
        datain.sphlap(:,lapchans,t)=datain.data(:,lapchans(datain.artifact(lapchans,t)==0),t)*mat';
    end
end