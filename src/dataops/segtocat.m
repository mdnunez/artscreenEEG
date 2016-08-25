function [dataout,tsamples] = segtocat(datain,goodtrials)
%function [dataout,[tsamples]] = segtocat(datain,[goodtrials])
%
% Concatenates data from 3D arrays (time by chan by trial) or 1 by ntrials 
% cell arrays (each cell containing a time by chan matrix) into a single 
% 2D array (time by chan).  
%
% If optional input goodtrials is used, only data from those trials is 
% returned.  Optional output tsamples is only needed for inequal trial
% lengths - used to properly resegment the data.  
%
% Copyright (C) 2014 Cort Horton, <chorton@uci.edu>
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% Version History
%   1.0 Created by Cort Horton 2/09
%   1.1 Now operates on 3D arays for cell arrays 6/13/14

if iscell(datain)
    % values read from data
    nchans=size(datain{1},2);
    ntrials=length(datain);
    
    % Specify good trials
    if nargin < 2;
        goodtrials=1:ntrials;
    end
    ngoodtrials=length(goodtrials);
    
    % Find number of samples in each trial
    tsamples=zeros(1,ngoodtrials);
    k=1;
    for t=goodtrials
        tsamples(k)=size(datain{t},1);
        k=k+1;
    end
    
    % Put all data into one matrix
    dataout=zeros(sum(tsamples),nchans);
    k=1;
    for t=goodtrials
        dataout(sum(tsamples(1:(k-1)))+1:sum(tsamples(1:k)),:) = datain{t};
        k=k+1;
    end
    
else
    tsamples=size(datain,1);
    nchans=size(datain,2);
    ntrials=size(datain,3);
    
    % Specify good trials
    if nargin < 2;
        goodtrials=1:ntrials;
    end
    ngoodtrials=length(goodtrials);
    
    dataout=zeros(tsamples*ngoodtrials,nchans);
    k=1;
    for t=goodtrials
        dataout(tsamples*(k-1)+1:tsamples*k,:)=datain(:,:,t);
        k=k+1;
    end
end