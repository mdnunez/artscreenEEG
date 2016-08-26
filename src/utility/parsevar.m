function [inProcessed,varargout] = parsevar(in,varargin)

%PARSEVAR parses varargin 'PropertyName', 'PropertyValue' pairs.
%   [processedVarargin,out1,out2,...,outN] = PARSEVAR(originalVarargin,...
%       parameter1,value1,p2,v2,...,pN,vN);
%   
%   PARSEVAR takes the original varargin, extract property name/value
%   pairs to seperate variables, and set default values to unsupplied
%   pairs if nessary. The processed varargin is then passed out. This could
%   be useful if only part of the varargin is to be censored and the rest
%   is to be fed to another function which also use property pairs.
% 
%   Note: all parameter1, ... ,N properties are ripped from varargin after
%   parseing, which means processedVarargin doesn't contain these pairs. If
%   reconstruction is necessary, use [processedVarargin,parameter1,out1,...];
%   
%   Example:
%       varargin = {'nation','MEX','gend',0,'ethn',...
%           'latino','birth',[1963 3 11],'name','LUIS'};
%       [varargin,nation,sex,race] = parsevar(varargin,...
%       'nation','USA','gender','1','ethnic','hispanic')
%
%   See also: VARARGIN, VARARGOUT.
%
% Copyright (C) 2007 Siyi Deng, <siyideng@live.com>
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
%
% Written by Siyi Deng; 05-03-2007;
% 05-26-07;

pair = varargin;
nPair = numel(pair);
if isempty(in), in = {sprintf('\n'),sprintf('\n')}; end
if ~iscell(in), in = {in}; end
if mod(nPair,2), error('Parameter and Value not in pair.'); end
if (nargout ~= nPair/2+1) & (nargout ~= 1), error('Argument number mismatch.'); end
nIn = numel(in);
inField = {in{1:2:end}};
inValue = {in{2:2:end}};
indexToDump = [];
for k = 1:max(nPair/2,1)
    %loc = strmatch(pair{2*k-1},inField);
    loc = [];
    tmp = pair{2*k-1};
    for j = 1:nIn/2
        tmp(end+1:length(inField{j})) = ' ';
        if strcmpi(tmp(1:length(inField{j})),inField{j})
                loc = [loc;j]; 
        end
    end        
    if length(loc) > 1, loc = strmatch(pair{2*k-1},inField,'exact'); end
    if isempty(loc)
        varargout{k} = pair{2*k};
    else
        varargout{k} = inValue{loc};
        indexToDump = [indexToDump,loc];
    end
end
inField(indexToDump) = [];
inValue(indexToDump) = [];
inProcessed = [inField;inValue];
inProcessed = inProcessed(:);
if strcmp(in{1},sprintf('\n')), inProcessed = {}; end

end % PARSEVAR;