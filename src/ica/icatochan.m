% function chandata = icatochan(icadata,evalthresh)
%
% Converts data from ICA space into channel space.  Input is a segmented
% data structure with fields (.ica, .data, or .att) and .mix.  
%
% Second optional argument is the minimum compeval to include
% when projecting the ICA components into channel space.  Def: 2
%   Good Only:2
%   Good + Unsure:1
%   All: 0
%
% Copyright (C) 2013 Cort Horton, <chorton@uci.edu>
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
function datain = icatochan(datain,evalthresh)

% VERSION HISTORY:
%   1.0 - Created by Cort Horton
%   1.1 - Updates to make compatble with old naming conventions 4/8/13
%   1.2 - Default option is to keep both Good and Unsure 12/27/17

if nargin <2; evalthresh=1; end;

% Zero out components that are not going to be projected back into channels
datain.mix(find(datain.compevals<evalthresh),:)=0;

% For typical use
if isfield(datain,'ica');
    nsamps=size(datain.ica,1);
    ntrials=size(datain.ica,3);
    nchans=size(datain.mix,2);
    
    datain.data=zeros(nsamps,nchans,ntrials);
    
    for t=1:ntrials
        datain.data(:,:,t)=datain.ica(:,:,t)*datain.mix;
    end
    datain=rmfield(datain,'ica');

% For legacy support
elseif isfield(datain,'data');
    nsamps=size(datain.data,1);
    ntrials=size(datain.data,3);
    nchans=size(datain.mix,2);
    
    temp=zeros(nsamps,nchans,ntrials);
    
    for t=1:ntrials
        temp(:,:,t)=datain.data(:,:,t)*datain.mix;
    end
    datain.data=temp;
    
% For use with dichotic listening cross-correlations
elseif isfield(datain,'att');
    nsamps=size(datain.att,1);
    ntrials=size(datain.att,3);
    nchans=size(datain.mix,2);
    
    temp1=zeros(nsamps,nchans,ntrials);
    temp2=zeros(nsamps,nchans,ntrials);
    temp3=zeros(nsamps,nchans,ntrials);
    temp4=zeros(nsamps,nchans,ntrials);
    temp5=zeros(nsamps,nchans,ntrials);
    
    for t=1:ntrials;
        temp1(:,:,t)=datain.att(:,:,t)*datain.mix;
        temp2(:,:,t)=datain.unatt(:,:,t)*datain.mix;
        temp3(:,:,t)=datain.cont(:,:,t)*datain.mix;
        temp4(:,:,t)=datain.left(:,:,t)*datain.mix;
        temp5(:,:,t)=datain.right(:,:,t)*datain.mix;
    end
    datain.att=temp1;
    datain.unatt=temp2;
    datain.cont=temp3;
    datain.left=temp4;
    datain.right=temp5;
end

% Remove unneeded fields
fieldstoremove={'sep' 'mix' 'cpvars' 'wspecs' 'compevals'};
for k=1:length(fieldstoremove)
    if isfield(datain,fieldstoremove{k});
        datain=rmfield(datain,fieldstoremove{k});
    end
end