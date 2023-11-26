function [U] = vortex( xe1, ye1, xe2, ye2, x, y) % it seems ok
l=sqrt((xe2-xe1)^2+(ye2-ye1)^2);
xl2 = l; yl2=0;
xl=(x-xe1)*(xe2-xe1)/l+(y-ye1)*(ye2-ye1)/l; yl=-(x-xe1)*(ye2-ye1)/l+(y-ye1)*(xe2-xe1)/l;
r1=sqrt(xl^2+yl^2); r2=sqrt((xl-xl2)^2+(yl-yl2)^2);
theta1=atan2(yl,xl); theta2=atan2(yl-yl2,xl-xl2);
if (abs(theta1)<10^(-12) && abs(theta2)>3); theta1=0; theta2=pi; end
if (abs(theta2)<10^(-12) && abs(theta1)>3); theta2=0; theta1=-pi; end
ul=1/(2*pi)*(theta2-theta1); vl=-1/(2*pi)*log(r1/r2);
u=ul*(xe2-xe1)/l-vl*(ye2-ye1)/l; v=ul*(ye2-ye1)/l+vl*(xe2-xe1)/l;
U=[u,v];
end