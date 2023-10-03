%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 1 - Problem 4 - Flow Field: Rotating Cylinder
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section sets up coordinates

xmin = -2;
xmax = 2;
dx = 0.1;

% X changes with i, Y changes with j
[Y, X] = meshgrid(xmin:dx:xmax);
ngrid = max(size(X));
x = linspace(xmin, xmax, ngrid);
y = linspace(xmin, xmax, ngrid);

% Freestream flow from left to right
vinf = 1;

% Negative gamma is clockwise rotation
Gamma = 0;

% Doublet strength dictates radius
C = 0.5;

% Find radius to each cartesian point
for i = 1:ngrid
    for j = 1:ngrid
        r(i,j) = sqrt(x(i)^2+y(j)^2);
        cos(i,j) = x(i)/r(i,j);
        sin(i,j) = y(j)/r(i,j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates velocity vectors and pressure contours for a doublet in uniform flow:
% Cylinder

for i = 1:ngrid
    for j = 1:ngrid
        vdubx(i,j) = -2*C*(cos(i,j)^2)/(r(i,j)^2);
        vduby(i,j) = -2*C*(cos(i,j)*sin(i,j))/(r(i,j)^2);
    end
end

for i = round(ngrid/2)
    for j = round(ngrid/2)
        vdubx(i,j) = 0;
        vduby(i,j) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates velocity vectors and pressure contours for a vortex

for i = 1:ngrid
    for j = 1:ngrid
        vvortx(i,j) = -Gamma/(2*pi*r(i,j))*(y(j)/r(i,j));
        vvorty(i,j) = Gamma/(2*pi*r(i,j))*(x(i)/r(i,j));
    end
end

vvortx(round(ngrid/2),round(ngrid/2)) = 0;
vvorty(round(ngrid/2),round(ngrid/2)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine velocity components for cylinder and vortex

for i = 1:ngrid
    for j = 1:ngrid
        vtotx(i,j) = vvortx(i,j)+vdubx(i,j)+vinf;
        vtoty(i,j) = vvorty(i,j)+vduby(i,j);
        lim = 2;

        if vtotx(i,j)>lim
            vtotx(i,j) = lim;
        end

        if vtotx(i,j)<-lim
            vtotx(i,j) = -lim;
        end

        if vtoty(i,j)>lim
            vtoty(i,j) = lim;
        end

        if vtoty(i,j)<-lim
            vtoty(i,j) = -lim;
        end
    end
end

for i = 1:ngrid
    for j = 1:ngrid
        vtot(i,j) = sqrt(vtotx(i,j)^2+vtoty(i,j)^2);
        cptot(i,j) = 1-((vtot(i,j)^2)/(vinf^2));
        lim = 25;

        if cptot(i,j)>lim
            cptot(i,j) = lim;
        end

        if cptot(i,j)<-lim
            cptot(i,j) = -lim;
        end
    end
end

cptot = cptot';
vtot = vtot';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp = min(min(vtot));
[row, col] = find(vtot<1.2*temp);

if size(row) == 2
    theta1 = atan(y(row(2))/x(col(2)));
    theta2 = (pi)+atan(y(row(1))/x(col(1)));
    theta1d = (180*theta1)/pi
    theta2d = (180*theta2)/pi
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)

subplot(2, 2, 1)
quiver(X, Y, vtotx, vtoty, 3)
axis square
axis([xmin, xmax, xmin, xmax])
colorbar
title('Flow Over Cylinder: Velocity Vectors')

subplot(2, 2, 2)
contour(x, y, vtot, 100)
axis square
axis([xmin, xmax, xmin, xmax])
colorbar
title('Flow Over Cylinder: Velocity Contours')

subplot(2, 2, 3)
contour(x,y,cptot,20)
axis square
axis([xmin, xmax, xmin, xmax])
colorbar
title('Flow Over Cylinder: Pressure Contours')

subplot(2, 2, 4)
contour(x, y, cptot, 20)
axis square
axis([xmin, xmax, xmin, xmax])
hold on
quiver(X, Y, vtotx, vtoty, 2)
colorbar
title('Flow Over Cylinder: Pressure Contours and Velocity Vectors')
hold off
