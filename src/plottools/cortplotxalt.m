function out = cortplotxalt(varargin)
% CORTPLOTXALT shows number of lineseries when clicked.
%   H = CORTPLOTXALT(...) accepts exactly same input arguments as PLOT, and
%   returns the handles of lineseries as vector H; click on any lineseries
%   object and the index will be shown in upper left corner.
%
% Optional Inputs: 
%     'custom': Custom labels, cell array of strings the same length as the second dimension
%
% Copyright (C) 2008 Siyi Deng, <siyideng@live.com>
% Copyright (C) 2013 Cort Horton, <chorton@uci.edu>
% Copyrigth (C) 2018 Michael D. Nunez, <mdnunez1@uci.edu>
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
% Modified by Cort Horton to make the numbers easier to read. Now has a
% larger font, thicker border, and a white background 7/6/13
% Modified by Michael D. Nunez to include custom labels 2/13/18

[varargin,custom]=parsevar(varargin,'custom','emptyfill');

h = plot(varargin{:});
hTxt = text('unit','norm','position',[.05 .9],...
    'string',' ','linewidth',3,'parent',get(h(1),'parent'),'fontsize',16,'margin',4);
if strcmp(custom,'emptyfill'),
	for k = 1:numel(h)
	    set(h(k),'HitTest','on','DisplayName',num2str(k),...
	        'ButtonDownFcn',{@lineseriescall,k,hTxt});
	end
else
	for k = 1:numel(h)
	    set(h(k),'HitTest','on','DisplayName',custom{k},...
	        'ButtonDownFcn',{@lineseriescustom,k,hTxt,custom});
	end
end

if nargout > 0, out = h; end

end % PLOTX;

% Nested Functions*************************************************************
	function lineseriescall(src,evt,k,h)
		set(h,'edgecolor',get(src,'color'),'backgroundcolor',[1 1 1],'string',num2str(k));
	end % end of function LINESERIESCALL;

	function lineseriescustom(src,evt,k,h,custom)
		set(h,'edgecolor',get(src,'color'),'backgroundcolor',[1 1 1],'string',custom{k});
	end % end of function LINESERIESCUSTOM;