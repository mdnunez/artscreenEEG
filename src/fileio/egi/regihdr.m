function [rate version n nevent N scale event header] = regihdr(fid)
% [rate version n nevent N scale event header] = regihdr(fid)
%
% reads in header versions 2 or 4 of NSFragger output files 
%    rate: sampling rate
% version: version of file
%       n: total number of channels
%  nevent: total number of events
%       N: total number of samples
%   scale: uV per unit
%   event: event codes 
%  header: header without event codes

% Adapted from rd_fragger_hdr
% Release 1.0 2/9/97 R.S. 
% 

% Copyright (C) 1997 Ramesh Srinivasan
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

if nargin < 1 || fid < 3, error('No valid file id.'); end
version = fread(fid,1,'int32');
if all(version ~= [2 4]), error('Incompatible version.'); end
year    = fread(fid,1,'int16');
month   = fread(fid,1,'int16');
day     = fread(fid,1,'int16');
hour    = fread(fid,1,'int16');
min     = fread(fid,1,'int16');
sec     = fread(fid,1,'int16');
msec    = fread(fid,1,'int32');
rate    = fread(fid,1,'int16');
n       = fread(fid,1,'int16');
gain    = fread(fid,1,'int16');
bits    = fread(fid,1,'int16');
range   = fread(fid,1,'int16');
N       = fread(fid,1,'int32');
nevent  = fread(fid,1,'int16');
event   = char(ones(nevent,4));
for k = 1:nevent, event(k,:) = fread(fid,[1 4],'char'); end
header  = [version year month day hour min sec msec rate n gain bits range N nevent];
scale   = range/(2^bits);