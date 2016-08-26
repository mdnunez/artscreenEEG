% function H = hortontopo(data,hm,varargin)
%
% Plots a 2-dimensional topography of given data using the supplied or
% named headmodel. If the data is complex, this calls complextopo.m
% instead. There are many options (described below) for how you want the
% topo plot to be drawn.
%
% Output (optional):
%   H - structure with handles and interpolated data
%
% Required Input:
%   data - the data vector to plot (nchans x 1)
%   hm - a headmodel with 3D channel positions following Siyi's structure
%        Example headmodels can be found in egihm.mat and antwavehm.mat.
%        Alternatively, you can pass a string with the name of an existing
%        system that has a headmodel, such as 'egi128' or 'antwave'
%
% Optional Inputs: (default values)
%   plotaxes - axes in which to create plot (current axes)
%   badchans - use this or goodchans, not both ([])
%   goodchans - use this or badchans, not both ([])
%
% Color Options
%   drawcolorbar - toggles drawing a colorbar (0)
%   cmap - the colormap to use ('jet')
%   clim - length 2 vector that sets the caxis limits (auto)
%   weights - toggles using the cmap hotncold and sets the clim equal to
%             [-1 1]*max(abs(data)). Useful for plotting channel weights
%             of ICA components or topos of condition differences. (0)
%   whiteback - sets the color of the masked area to white.  This looks
%               bad against the default grey figure color, but without
%               this set to 1 your saved figures will have a grey square (0)
%
% Electrode Options
%   drawelectrodes - toggles drawing electrodes on head (1)
%   drawoutliers - toggles drawing electrodes outsize the circle (0)
%   channumbers - toggles plotting channel numbers instead of dots (1)
%   chanfontsize - fontsize for channel numbers (10)
%   markbadchans - draws bad electrodes dark (1)
%   goodelectcolor - triplet color for the good electrodes ([.9,.9,.9])
%   badelectcolor - triplet color for the bad electrodes ([.2,.2,.2])
%   colorthese - colors a subset of electrodes ([])
%   otherelectcolor - triplet color for the subset of electrodes ([.5 .5 0])
%
% Contour Options
%   drawcontours - toggles drawing contour lines (0)
%   ncontours - number of countour lines (6)
%
% Advanced Options
%   nPixel - resolution of the plot (256)
%   electrodeWidth - lower this for very high density recordings (.01)
%
% Copyright (C) 2006 Bill Winter, <wwinter@hs.uci.edu>
% Copyright (C) 2007 Siyi Deng, <siyideng@live.com>
% Copyright (C) 2013 Cort Horton, <chorton@uci.edu>
% Copyright (C) 2016 Michael D. Nunez, <mdnunez1@uci.edu>
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

function HH = hortontopo(data,hm,varargin)
% Version History:
%   1.0 - Originally adapted from Siyi Deng's topohead2d.m   All further
%          revisions by Cort Horton
%   1.1 - Fixed errors with plotaxes, head models, alpha mapping
%   1.2 - Added features: weights flag, bad channel marking
%   1.3 - No longer hard-coded to a recording environment. Headmodel is
%          now expected as an input
%   1.4 - Can now plot channel numbers instead of dots
%   1.5 - Allowed faux transparency with non-opengl rendereres by setting
%          the lowest colormap value to the figure background color. Also
%          added the ability to name a headmodel to use.  7/12/13
%   1.6 - Changed name of function to reduce confusion with Bill's topo
%         plotting functions. Also now calls complextopo when it detects
%         complex data. 7/14/13
%   1.7 - Include color electrode indications - Michael Nunez on 7/01/14
%   1.8 - Added dependencies - Michae Nunez on 8/25/16

% Transposes data into column vector if needed
if size(data,1)==1; data=data'; end;

% Detects complex input and calls complextopo.m instead
if ~isreal(data);
    warning('''complextopo'' does not exist in the artscreenEEG repository as of 8/23/16!\n')
    H = complextopo(data,hm,varargin{:});
else
    % Parse inputs;
    [~,nPixel,cmap,badchans,goodchans,drawcontours,ncontours,drawelectrodes,...
        drawoutliers,channumbers,chanfontsize,markbadchans,goodelectcolor,...
        badelectcolor,drawcolorbar,plotaxes,weights,electrodeWidth,clim,...
        colorthese,otherelectcolor,colorthese2,otherelectcolor2,...
        colorthese3,otherelectcolor3, colorthese4,otherelectcolor4,...
        colorthese5,otherelectcolor5,whiteback] = ...
        parsevar(varargin,'nPixel',256,'cmap',[],'badchans',[],'goodchans',[],...
        'drawcontours',0,'ncontours',6,'drawelectrodes',1,'drawoutliers',0,...
        'channumbers',1,'chanfontsize',10,...
        'markbadchans',1,'goodelectcolor',[1 1 1],'badelectcolor',[.3 .3 .3],...
        'drawcolorbar',0,'plotaxes',gca,'weights',0,'electrodeWidth',.01,...
        'clim',[],'colorthese',[],'otherelectcolor',[0 0 1],...
        'colorthese2',[],'otherelectcolor2',[0 .9 .4],...
        'colorthese3',[],'otherelectcolor3',[1 .3 0],...
        'colorthese4',[],'otherelectcolor4',[0 .7 .7],...
        'colorthese5',[],'otherelectcolor5',[.3 .6 1],...
        'whiteback',0);
  
    
    % Default colormaps
    if weights && isempty(cmap);
        cmap='hotncold';
    elseif isempty(cmap);
        cmap='jet';
    end
    
    if isempty(colorthese)
        markthese = 0;
    else
        markthese = 1;
    end
    if isempty(colorthese2)
        markthese2 = 0;
    else
        markthese2 = 1;
    end
    if isempty(colorthese3)
        markthese3 = 0;
    else
        markthese3 = 1;
    end
    if isempty(colorthese4)
        markthese4 = 0;
    else
        markthese4 = 1;
    end
    if isempty(colorthese5)
        markthese5 = 0;
    else
        markthese5 = 1;
    end
    
    % loads a named headmodel
    if ischar(hm);
        tmp=load([hm 'hm']);
        fn=fieldnames(tmp);
        hm=getfield(tmp,fn{1});
    end
    
%     % Remove channels without coordinates (ie EOG channels) from the data
%     if length(data) > size(hm.Electrode.CoordOnSphere,1);
%         data=data(1:size(hm.Electrode.CoordOnSphere,1));
%     end
    
    % Resolve good and bad channel lists
    nchans=length(data);
    if isempty(goodchans) && isempty(badchans);
        goodchans=1:nchans;
    elseif isempty(goodchans) && ~isempty(badchans);
        goodchans=setdiff(1:nchans,badchans);
    elseif ~isempty(goodchans) && isempty(badchans);
        goodchans=intersect(1:nchans,goodchans);
        badchans=setdiff(1:nchans,goodchans);
    elseif ~isempty(goodchans) && ~isempty(badchans);
        error('Use a goodchan list or a badchan list, not both');
    end
    badchans(badchans > nchans) = [];
    
    % Other preprocessing
    es = hm.Electrode.CoordOnSphere;
    H = struct;
    if size(data,1)==1;
        data=data';
    end
    
    % Remove bad channels;
    data(badchans) = [];
    es(badchans,:) = [];
    
    % Data interpolation from 3D to 2D
    [xi,yi,zi] = plane2ball(nPixel,hm.Electrode.SphereExtendAngle,'AED');
    vi = interpolate3d(es,data,[xi(:),yi(:),zi(:)],electrodeWidth);
    
    % Restrict image to a disc; get the range of data;
    diameter= nPixel;
    vi(logical(1-discmask(nPixel,diameter-2))) = NaN;
    vi = reshape(vi,nPixel,nPixel);
    H.Data=vi;
    
    % Get limits of the data
    dataRange = [min(min(vi)) max(max(vi))];
    
    % Prepare the figure and axes;
    H.Axes=plotaxes;
    cla(H.Axes);
    axis(H.Axes,'equal','off');
    hold(H.Axes,'on');
    
    % Construct the colormap;
    if ischar(cmap);
        c = eval([cmap '(256);']);
    else
        c=cmap;
    end
    
    % Forces opengl renderer for transparency
    set(get(plotaxes,'parent'),'renderer','opengl');
    
%     % Needed to simulate tranparency without opengl
%     figrenderer=get(get(plotaxes,'parent'),'renderer');
%     if ~strcmp(figrenderer,'OpenGL')
%         if whiteback
%             c(256,:)=[1 1 1];
%         else
%             c(256,:)=get(get(plotaxes,'parent'),'color');
%         end
%     end
    
    % Set the colormap
    H.Colormap = c;
    colormap(H.Axes,H.Colormap);
    
    % Draw the interpolated data
    H.Image = image(vi,'Parent',H.Axes,'CDataMapping','scaled');
    
    % Set the color axis
    if isempty(clim)
        if weights
            clim=[-1.05 1.05]*max(abs(dataRange));
        else
            clim=dataRange;
        end
    end
    caxis(H.Axes,clim);
    
    % Draw the colorbar;
    if drawcolorbar
        H.Colorbar = colorbar;
    end
    
    % Draw contour lines;
    if drawcontours
        H.ContourStep = linspace(dataRange(1),...
            dataRange(2),ncontours);
        [~,H.ContourLow] = ...
            contour(H.Data,H.ContourStep(1:ceil(ncontours/2)),...
            'LineStyle','--','color',[0.2,0.2,0.2]);
        [~,H.ContourHigh] = ...
            contour(H.Data,...
            H.ContourStep(ceil(ncontours/2+1):ncontours),...
            'LineStyle','-','color',[0.2,0.2,0.2]);
    end
    
    % Draw the ear and nose sketch;
    [x,y] = pol2cart([linspace(-pi/4,pi/4,20);...
        linspace(0.75*pi,1.25*pi,20)],diameter/6);
    y = y+nPixel/2;
    x = x+nPixel/2+sign(x)*6*diameter/16+1/2;
    H.Ears = plot(H.Axes,x',y','k-','linewidth',3);
    x = linspace(-nPixel/10,nPixel/10,20);
    y = nPixel*cos(x*5*pi/nPixel)/12+nPixel-1-ceil((nPixel-diameter)/2);
    x = x+nPixel/2;
    H.Nose = plot(H.Axes,x,y,'k-','linewidth',3);
    
    [x,y] = pol2cart(linspace(-pi/2,1.5*pi,128),diameter/2-1);
    H.Ring = plot(H.Axes,x+nPixel/2+0.5,y+nPixel/2+0.5,'k-','linewidth',3);
    
    % Transparency effects
    H.Mask = isnan(vi);
%     if strcmp(figrenderer,'OpenGL')
        set(H.Image,'AlphaData',double(~H.Mask));
%     else
%         tmp=get(H.Image,'CData');
%         tmp(H.Mask)=clim(2);
%         set(H.Image,'CData',tmp);
% %         if ~weights
% %             tmp=get(H.Image,'CData');
% %             tmp(~H.Mask)=tmp(~H.Mask)+(clim(2)-clim(1))/30;
% %             set(H.Image,'CData',tmp);
% %         end
%     end
    
    % Draw electrodes;
    if drawelectrodes
        [x,y] = ball2plane(hm.Electrode.CoordOnSphere(:,1),...
            hm.Electrode.CoordOnSphere(:,2),hm.Electrode.CoordOnSphere(:,3));
        
        if ~drawoutliers
            outlier = ((x.^2+y.^2) > 1.05*hm.Electrode.Scaling2D^2/4);
            x(outlier) = NaN;
            y(outlier) = NaN;
        end
        
        x = x*nPixel/hm.Electrode.Scaling2D+nPixel/2;
        y = y*nPixel/hm.Electrode.Scaling2D+nPixel/2;
        
        if channumbers
            for c = 1:nchans, thestrings{c} = num2str(c); end;
            if markbadchans
                H.gctxt = text(x(setdiff(goodchans,colorthese)),y(setdiff(goodchans,colorthese)),...
                    ones(size(setdiff(goodchans,colorthese))),thestrings(setdiff(goodchans,colorthese)),...
                    'Parent',H.Axes,'HorizontalAlignment','center');
                set(H.gctxt,'fontsize',chanfontsize,'fontweight','bold','color',goodelectcolor);
                H.bctxt = text(x(badchans),y(badchans),ones(size(badchans)),thestrings(badchans),...
                    'Parent',H.Axes,'HorizontalAlignment','center');
                set(H.bctxt,'fontsize',chanfontsize,'fontweight','bold','color',badelectcolor);
            else
                H.txt = text(x,y,ones(size(x)),thestrings);
                set(H.txt,'fontsize',chanfontsize,'fontweight','bold','color',goodelectcolor,'Parent',H.Axes,...
                    'HorizontalAlignment','center');
            end
            if markthese
                H.gctxt = text(x(colorthese),y(colorthese),ones(size(colorthese)),thestrings(colorthese),'Parent',H.Axes,...
                    'HorizontalAlignment','center');
                set(H.gctxt,'fontsize',chanfontsize,'fontweight','bold','color',otherelectcolor);
            end
            if markthese2
                H.gctxt = text(x(colorthese2),y(colorthese2),ones(size(colorthese2)),thestrings(colorthese2),'Parent',H.Axes,...
                    'HorizontalAlignment','center');
                set(H.gctxt,'fontsize',chanfontsize,'fontweight','bold','color',otherelectcolor2);
            end
            if markthese3
                H.gctxt = text(x(colorthese3),y(colorthese3),ones(size(colorthese3)),thestrings(colorthese3),'Parent',H.Axes,...
                    'HorizontalAlignment','center');
                set(H.gctxt,'fontsize',chanfontsize,'fontweight','bold','color',otherelectcolor3);
            end
            if markthese4
                H.gctxt = text(x(colorthese4),y(colorthese4),ones(size(colorthese4)),thestrings(colorthese4),'Parent',H.Axes,...
                    'HorizontalAlignment','center');
                set(H.gctxt,'fontsize',chanfontsize,'fontweight','bold','color',otherelectcolor4);
            end
            if markthese5
                H.gctxt = text(x(colorthese5),y(colorthese5),ones(size(colorthese5)),thestrings(colorthese5),'Parent',H.Axes,...
                    'HorizontalAlignment','center');
                set(H.gctxt,'fontsize',chanfontsize,'fontweight','bold','color',otherelectcolor5);
            end
        else
            if markbadchans
                H.GoodElectrode = plot(H.Axes,x(setdiff(goodchans,colorthese)),y(setdiff(goodchans,colorthese)),...
                    'Marker','o','MarkerEdgeColor',goodelectcolor,...
                    'MarkerFaceColor',goodelectcolor,...
                    'MarkerSize',4,'LineStyle','none');
                H.BadElectrode = plot(H.Axes,x(badchans),y(badchans),'Marker','o',...
                    'MarkerEdgeColor',badelectcolor,...
                    'MarkerFaceColor',badelectcolor,...
                    'MarkerSize',4,'LineStyle','none');
            else
                H.Electrode = plot(H.Axes,x,y,'Marker','o',...
                    'MarkerEdgeColor',goodelectcolor,...
                    'MarkerFaceColor',goodelectcolor,...
                    'MarkerSize',4,'LineStyle','none');
            end
            if markthese
                H.GoodElectrode = plot(H.Axes,x(colorthese),y(colorthese),'Marker','o',...
                    'MarkerEdgeColor',otherelectcolor,...
                    'MarkerFaceColor',otherelectcolor,...
                    'MarkerSize',4,'LineStyle','none');
            end
            if markthese2
                H.GoodElectrode = plot(H.Axes,x(colorthese2),y(colorthese2),'Marker','o',...
                    'MarkerEdgeColor',otherelectcolor2,...
                    'MarkerFaceColor',otherelectcolor2,...
                    'MarkerSize',4,'LineStyle','none');
            end
            if markthese3
                H.GoodElectrode = plot(H.Axes,x(colorthese3),y(colorthese3),'Marker','o',...
                    'MarkerEdgeColor',otherelectcolor3,...
                    'MarkerFaceColor',otherelectcolor3,...
                    'MarkerSize',4,'LineStyle','none');
            end
            if markthese4
                H.GoodElectrode = plot(H.Axes,x(colorthese4),y(colorthese4),'Marker','o',...
                    'MarkerEdgeColor',otherelectcolor4,...
                    'MarkerFaceColor',otherelectcolor4,...
                    'MarkerSize',4,'LineStyle','none');
            end
            if markthese5
                H.GoodElectrode = plot(H.Axes,x(colorthese5),y(colorthese5),'Marker','o',...
                    'MarkerEdgeColor',otherelectcolor5,...
                    'MarkerFaceColor',otherelectcolor5,...
                    'MarkerSize',4,'LineStyle','none');
            end
        end
    end
    
end

if nargout>0; HH=H; end

% Nested Functions*************************************************************
    function [x,y,z] = plane2ball(n,angle,method)

        %PLANE2BALL creates coordinates on a 3-D unit ball.
        %   [X,Y,Z] = PLANE2BALL(N) creates coordinates on a unit sphere
        %   which, when projected to a plane, give N x N matrices;
        %   
        %   [X,Y,Z] = PLANE2BALL(N,ANGLE) specifes the angle the spherical
        %   coordinate extend; default ANGLE = pi;
        %
        %   [X,Y,Z] = PLANE2BALL(N,ANGLE,METHOD) uses METHOD to perform the inverse
        %   projection; default METHOD = 'AED'; METHOD could be:
        %       'AED': Azimuthal Equidistant Projection from center;
        %       'STG': Stereographic Projection from remote pole;
        %       'ORT': Orthographic Projection from infinity; note that angles
        %              grather than PI are treated as 2*PI-angle;
        %
        %   Example:
        %       [x,y,z] = plane2ball(20);
        %       plot3(x,y,z,'b.');
        %       axis equal vis3d off;

        % Written by Siyi; 04-22-2007;

        if nargin < 2 | isempty(angle), angle = pi; end
        if nargin < 3, method = 'AED'; end

        if strcmpi(method,'AED')
            [x,y] = meshgrid(linspace(-angle/2,angle/2,n));
            [t,r] = cart2pol(x,y);
            [x,y,z] = sph2cart(t,pi/2-r,1);
        elseif strcmpi(method,'STG')
            [x,y] = meshgrid(linspace(-tan(angle/4),tan(angle/4),n));
            [t,r] = cart2pol(x,y);
            [x,y,z] = sph2cart(t,pi/2-2*atan(r),1);
        elseif strcmpi(method,'ORT')
            [x,y] = meshgrid(linspace(-sin(angle/2),sin(angle/2),n));
            [t,r] = cart2pol(x,y);
            [x,y,z] = sph2cart(t,acos(r),1);
        else
            error('Unknown projection method.');
        end

    end % end of function PLANE2BALL;

    function [x1,y1] = ball2plane(x,y,z)

        %BALL2PLANE projects the coordinates from a ball to a plane.
        %   [X1,Y1] = BALL2PLANE(X,Y,Z) uses the Azimuthal Equal-Distant projection
        %   method to project X,Y,Z onto X1,Y1;
        %
        %   See also PLANE2BALL.

        % Written by Siyi; 04-22-2007;

        [t,p,r] = cart2sph(x,y,z);
        p1 = min(p);
        p2 = min(p(p~=p1));
        angle = pi-p1-p2;
        [x1,y1] = pol2cart(t,(pi/2-p)*max(r)*2/angle);

    end % end of function BALL2PLANE;


    function vi = interpolate3d(xyzD,vd,xyzI,w)
        %INTERPOLATE3D 3-D spline interpolation;
        %   VI = INTERPOLATE3D([X,Y,Z],V,[XI,YI,ZI],W) uses min-region radius W, to
        %   interpolate values V on points defined in euclidean space [X, Y, Z]
        %   at interpolation points [XI, YI, ZI];
        %
        %   See also INTERP_3D, MATEQS, MATMAKE, K_AND_E, OFCP.

        % Written by Siyi Deng; 04-25-2007;
        [xd,yd,zd,xi,yi,zi] = deal(xyzD(:,1),xyzD(:,2),xyzD(:,3),...
            xyzI(:,1),xyzI(:,2),xyzI(:,3));
        if nargin < 4, w = (max(yi(:))-min(yi(:)))/20; end

        ww = w.^2;
        % ed = zeros(length(xd),10);
        % ed = [1+rand(size(xd))*(1e-9),... % In order to avoid rank deficiency of ed;
        %     xd,yd,xd.^2,xd.*yd,yd.^2,zd,zd.*xd,zd.*yd,zd.^2];
        ed = [ones(size(xd)),... % In order to avoid rank deficiency of ed;
            xd,yd,xd.^2,xd.*yd,yd.^2,zd,zd.*xd,zd.*yd,zd.^2];
        kd = eucdist([xd,yd,zd]); % disteu is bill's funcion;

        % Is this necessary?
        %kd(kd < w/2) = w/2;

        kd = kd.^2 + ww;
        kd = log(kd).*(kd.^2);

        % kinv = inv(kd);
        % qd = inv(ed'*(kinv*ed))*(ed'*(kinv*vd));
        % pd = kinv*vd-(kinv*(ed*qd));
        kinv = pinv(kd);
        qd = pinv(ed'*(kinv*ed))*(ed'*(kinv*vd));
        pd = kinv*vd-(kinv*(ed*qd));

        % qd = ed\vd;
        % pd = kd\(vd-ed*qd);

        vi = zeros(size(xi));

        % Based on Ramesh's algorism;
        for id = 1:numel(xd)
            kid = ww+d2(xd(id),yd(id),zd(id),xi,yi,zi);
            vi = vi+pd(id)*(log(kid).*(kid.^2));
        end

        % Osculating function serves to smooth the interpolated result;
        vi = vi+qd(1)+qd(2).*xi+qd(3).*yi+qd(4).*(xi.^2)+...
            qd(5).*xi.*yi+qd(6).*(yi.^2)+qd(7).*zi+...
            qd(8).*zi.*xi+qd(9).*zi.*yi+qd(10).*(zi.^2);

    end % end of function INTERPOLATE3D;


    function out = d2(x0,y0,z0,x,y,z)
        % Squared distance to a certain point;
        out = (x-x0).^2+(y-y0).^2+(z-z0).^2;
    end % end of function D2;


    function mat = eucdist(x)
        %EUCDIST finds the euclidean distance between points.
        % Slightly modified DISTEU, written by Bill Winter, December 2005;
        [x1,y1] = meshgrid(sum(x.^2,2));
        mat = real(sqrt(x1+y1-2*x*x'));
        mat(mat<length(mat)*eps(norm(mat))) = 0;
        mat(1:length(x)+1:numel(x)) = 0;
    end % end of function EUCDIST;

    function out = discmask(n,d)

        %DISCMASK creates a N x N square mask inscribed by a disc with diameter D.
        %   M = DISCMASK(N, D)
        %   M is a N x N matrix of 0, inscribed by a disc with value 1;
        %   N should be equal to or larger than D; by default N equals D.
        %
        %   Example:
        %       imagesc(discmask(131)); axis image;

        % Written by Siyi Deng; 04-03-2007;

        if nargin < 2, d = n; end;

        % Generate the coordinate vectors x, y;
        [x,y] = meshgrid(linspace(-n/2,n/2,n));

        out = (x.^2+y.^2 < d*d/4);

    end % DISCMASK;

    function h = hotncold(m)
        % h = hotncold(m)
        %
        % HOTNCOLD returns a colormap, cycling through white, cyan, blue, black,
        % red, yellow, white.
        %
        % m: length of colormap.  default is length of current colormap
        % h: mx3 rgb colormap
        %
        % See also: COLORMAP
        % Written by Bill Winter, March 2006

        if nargin < 1, m = size(get(gcf,'colormap'),1); end
        h = symmap(hot(m),1); 
    end % end of hotncold

    % Written by Bill Winter, March 2006
    function h = symmap(h,flip)
        % h = symmap(map,flip)
        %
        % SYMMAP operates on a colormap 'map' and creates a symmetric map of equal
        % length.  With 'flip' true, the red/blue components are reversed in the
        % lower map, as in hotncold.
        %
        %  map: colormap.  default is current colormap
        % flip: whether to flip the red/blue component in the lower half of the map
        %
        %    h: mx3 rgb colormap
        %
        % See also: COLORMAP, HOTNCOLD
        if nargin < 1, h = colormap; end
        if nargin < 2, flip = 0; end
        h = [h(2*floor(end/2):-2:1,:);h(1:2:end,:)];
        if flip, h(1:floor(end/2),[1 3]) = h(1:floor(end/2),[3 1]); end
    end

end % end of main function