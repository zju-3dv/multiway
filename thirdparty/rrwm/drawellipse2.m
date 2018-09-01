%%%%%%%DISPLAY ELLIPSE
function drawellipse2( xyabc ,scaling,col,width)

x = xyabc(1); y = xyabc(2); a = xyabc(3); b = xyabc(4); c = xyabc(5);

if( ~exist('width') )
    width = 1;
end

[v e]=eig([a b;b c]);

l1=1/sqrt(e(1));
l2=1/sqrt(e(4));

alpha=atan2(v(4),v(3));
t = 0:pi/50:2*pi;
yt=scaling*(l2*sin(t));
xt=scaling*(l1*cos(t));

p=[xt;yt];
R=[cos(alpha) sin(alpha);-sin(alpha) cos(alpha)];
pt=R*p;
plot(pt(2,:)+x,pt(1,:)+y,'Color',col,'LineWidth',width);
%set(gca,'Position',[0 0 1 1]);