function [g z] = regigz(file,varargin)
% [gain zero] = regigz(file)
%
% reads gain and zero file corresponding to rawfname

% Copyright (C) 2006 William Winter
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


env = environment(varargin{:});
if ~strcmp(env,{'egi' 'hnl'}),
    warning(['Environment "' env '" is not egi-based.']);
    env = 'hnl';
end
f = strfind(lower(file),'.raw');
file(f:end) = [];                               % drop extension

fid = fopen([file '.GAIN']);
if fid == -1, fid = fopen([env '.GAIN']); end   % try file
if fid == -1, error('invalid environment'); end % try default
flook(fid,'Gain event %s');
g = fscanf(fid,'%*s : %e',inf);
fclose(fid);

fid = fopen([file '.ZERO']);
if fid == -1, fid = fopen([env '.ZERO']); end   % try file
if fid == -1, error('invalid environment'); end % try default
flook(fid,'Zero event %s');
z = fscanf(fid,'%*s : %e',inf);
fclose(fid);
z(129) = 0;