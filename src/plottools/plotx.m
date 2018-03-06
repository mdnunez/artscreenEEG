function out = plotx(varargin)
%PLOTX shows number of lineseries when clicked.
%   H = PLOTX(...) accepts exactly same input arguments as PLOT, and
%   returns the handles of lineseries as vector H; click on any lineseries
%   object and the index will be shown in upper left corner.

% Copyright (C) 2008 Siyi Deng, <siyideng@live.com>
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
% Written by Siyi Deng; 08-22-2008;

h = plot(varargin{:});
hTxt = text('unit','norm','position',[.05 .9],...
    'string',' ','edgecolor','none','linewidth',2);
for k = 1:numel(h)
    set(h(k),'HitTest','on','DisplayName',num2str(k),...
        'ButtonDownFcn',{@lineseriescall,k,hTxt});
end    
if nargout > 0, out = h; end
end % PLOTX;

function lineseriescall(src,evt,k,h)
set(h,'edgecolor',get(src,'color'),'string',num2str(k));
end % LINESERIESCALL;