function [out index] = regievt(fid,rate)
% out = regievt(filename)
%
% REGIEVT reads EGI event files and returns the salient content
%
% out: structure array containing event labels and event times

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

fid = fopen(fid);
if fid < 3, error('inappropriate file'); end
for k = 1:3, fgetl(fid); end
label = textscan(fid,'%s %*s %*s %*s _%n:%n:%n %*s','Delimiter','\t');
fclose(fid);

index = (label{2}*60 + label{3})*60 + label{4};
if nargin > 1, index = floor(index*rate); end
[label,~,in] = unique(label{1});

% see parsetechmarkup below for multi-part files
for k = 1:length(label), out.(label{k}) = index(in == k); end

% [names,~,nin] = unique(Tech_Markup(1,:));
% [session,~,sin] = unique(cat(2,Tech_Markup{3,:}));
% 
% index = cat(1,Tech_Markup{4,:});
% 
% for s = 1:length(session)
%     for n = 1:length(names)
%         out(s).(names{n}) = index(nin == n & sin == s);
%     end
% end
