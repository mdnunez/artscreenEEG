function out = cortplotx(varargin)
% CORTPLOTX shows number of lineseries when clicked.
%   H = CORTPLOTX(...) accepts exactly same input arguments as PLOT, and
%   returns the handles of lineseries as vector H; click on any lineseries
%   object and the index will be shown in upper left corner.
%
% Copyright (C) 2008 Siyi Deng, <siyideng@live.com>
% Copyright (C) 2013 Cort Horton, <chorton@uci.edu>
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
% Written by Siyi Deng; 08-22-2008;
% Modified by Cort Horton to make the numbers easier to read. Now has a
% larger font, thicker border, and a white background 7/6/13

h = plot(varargin{:});
hTxt = text('unit','norm','position',[.05 .9],...
    'string',' ','linewidth',3,'parent',get(h(1),'parent'),'fontsize',16,'margin',4);
for k = 1:numel(h)
    set(h(k),'HitTest','on','DisplayName',num2str(k),...
        'ButtonDownFcn',{@lineseriescall,k,hTxt});
end

if nargout > 0, out = h; end

end % PLOTX;

% Nested Functions*************************************************************
	function lineseriescall(src,evt,k,h)
		set(h,'edgecolor',get(src,'color'),'backgroundcolor',[1 1 1],'string',num2str(k));
	end % end of function LINESERIESCALL;
