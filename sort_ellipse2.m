function [state] = sort_ellipse2(p, ellipse_param)
%sort_ellipse.
%This function determines if one x,y point is inside a defined ellipse.
%Inputs:
%p: it is a vector containing x and y coordinates of a point p.
%ellipse_param: it is a vector containing the ellipse parameters: center
%coordinates, ellipse radii and main ellipse axe angle.
%Output:
%state: its 1 if the point is inside the ellipse and 0 if it is not.

%Value extraction
%Point coordinates:
x = p(1);
y = p(2);
%Coordinates of the center of the ellipse:
xc = ellipse_param(1);
yc = ellipse_param(2);
%Ellipse radii:
A = ellipse_param(3);
B = ellipse_param(4);

%Main ellipse axe angle:
theta = ellipse_param(5);

xt = x - xc;
yt = y - yc;

xr = xt*cos(theta) + yt*sin(theta);
yr = -xt*sin(theta) + yt*cos(theta);

if xr^2/A^2 + yr^2/B^2 <= 1
    state = 1;
else
    state = 0;
end