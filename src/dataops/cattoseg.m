function dataout = cattoseg(datain,tsamps,mat)
% function dataout = cattoseg(datain,tsamps,[mat])
%
% Takes a 2D array (i.e. sample by channel) and reshapes it into a 3D
% array (i.e. sample by channel by trial). tsamps indicates the size
% of the first dimension (i.e. samples per trial).
%
% Optionally, if you also pass in a matrix, it will multiply the data 
% through that matrix.  This can be used to apply spatial filters, etc.
%
% Copyright (C) 2009 Cort Horton, <chorton@uci.edu>
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
totalsamps=size(datain,1);
if rem(totalsamps,tsamps)~=0
	disp('Warning: Data length is not an integer multiple of tsamps.');
	disp('         Samples beyond the last multiple will be dropped.');
end

ntrials=floor(totalsamps/tsamps);
nchans=size(datain,2);

if nargin <3; mat=eye(nchans); end;

dataout=zeros(tsamps,nchans,ntrials);

for t=1:ntrials
    dataout(:,:,t)=datain(tsamps*(t-1)+1:tsamps*t,:)*(mat*mat');
end
