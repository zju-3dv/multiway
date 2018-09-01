function [tout,rout] = roseBar(theta_bins, hist1, hist2)
% by Minsu Cho
hold on;
rmax = max( max(max(hist1),max(hist2)), eps);
drawCircle(0,0, rmax, 'k-');
cst = cos(theta_bins);
snt = sin(theta_bins);
for k=1:numel(theta_bins)
    plot( [0 rmax * cst(k)], [0 rmax * snt(k)], 'k:', 'linewidth', 1);
end
axis ij square;

rt = 1.1*rmax;
for i = 1:numel(theta_bins)
%    text(rt*cst(i),rt*snt(i),int2str(i*30),...
%         'horizontalalignment','center',...
%         'handlevisibility','off','parent',cax);
    loc = sprintf('%.0f',theta_bins(i)*180/pi);
    %text(ceil(rt*cst(i)),ceil(rt*snt(i)),loc);
end

% first bar (wide)
for i = 1:numel(theta_bins)
    if i == 1
        x = hist1(i)*[ 0 cst(end) cst(i)];
        y = hist1(i)*[ 0 snt(end) snt(i)];
    else
        x = hist1(i)*[ 0 cst(i-1) cst(i)];
        y = hist1(i)*[ 0 snt(i-1) snt(i)];
    end
    fill(x,y,'b')
end

% second bar (narrow)
d_theta = (theta_bins(3) - theta_bins(2))/4;
for i = 1:numel(theta_bins)
    if i == 1
        x = hist2(i)*[ 0 cos(theta_bins(end)+d_theta) cos(theta_bins(i)-d_theta)];
        y = hist2(i)*[ 0 sin(theta_bins(end)+d_theta) sin(theta_bins(i)-d_theta)];
    else
        x = hist2(i)*[ 0 cos(theta_bins(i-1)+d_theta) cos(theta_bins(i)-d_theta)];
        y = hist2(i)*[ 0 sin(theta_bins(i-1)+d_theta) sin(theta_bins(i)-d_theta)];
    end
    fill(x,y,'r')
end

xlim([-rmax rmax]);
ylim([-rmax rmax]);
        
hold off;
return;
%     
% [cax,args,nargs] = axescheck(varargin{:});
% error(nargchk(1,2,nargs,'struct'));
% 
% theta = args{1};
% if nargs > 1, 
%   x = args{2}; 
% end
% 
% if ischar(theta)
%   error(id('NonNumericInput'),'Input arguments must be numeric.');
% end
% %theta = rem(rem(theta,2*pi)+2*pi,2*pi); % Make sure 0 <= theta <= 2*pi
% if nargs==1,
%   x = (0:19)*pi/10+pi/20;
% 
% elseif nargs==2,
%   if ischar(x)
%     error(id('NonNumericInput'),'Input arguments must be numeric.');
%   end
%   if length(x)==1,
%     x = (0:x-1)*2*pi/x + pi/x;
%   else
%     x = sort(rem(x(:)',2*pi));
%   end
% 
% end
% if ischar(x) || ischar(theta)
%   error(id('NonNumericInput'),'Input arguments must be numeric.');
% end
% 
% % Determine bin edges and get histogram
% edges = sort(rem([(x(2:end)+x(1:end-1))/2 (x(end)+x(1)+2*pi)/2],2*pi));
% edges = [edges edges(1)+2*pi];
% % nn = histc(rem(theta+2*pi-edges(1),2*pi),edges-edges(1));
% % nn(end-1) = nn(end-1)+nn(end);
% % nn(end) = [];
% nn = theta(:);
% 
% % Form radius values for histogram triangle
% if min(size(nn))==1, % Vector
%   nn = nn(:); 
% end
% [m,n] = size(nn);
% mm = 4*m;
% r = zeros(mm,n);
% r(2:4:mm,:) = nn;
% r(3:4:mm,:) = nn;
% 
% % Form theta values for histogram triangle from triangle centers (xx)
% zz = edges;
% 
% t = zeros(mm,1);
% t(2:4:mm) = zz(1:m);
% t(3:4:mm) = zz(2:m+1);
% 
% if nargout<2
%   if ~isempty(cax)
%     h = polar(cax,t,r);
%   else
%     h = polar(t,r);
%   end
%   
%   % Register handles with m-code generator
%   %if ~isempty(h)
%   %   mcoderegister('Handles',h,'Target',h(1),'Name','rose');
%   %end
%   
%   if nargout==1, tout = h; end
%   return
% end
% 
% if min(size(nn))==1,
%   tout = t'; rout = r';
% else
%   tout = t; rout = r;
% end
% 
% function str=id(str)
% str = ['MATLAB:rose:' str];
