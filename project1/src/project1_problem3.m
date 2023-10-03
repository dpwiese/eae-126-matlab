%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 1 - Problem 3 - Flow Field: Rankine Body, Kelvin Oval, Cylinder
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section sets up coordinates

xmin = -1;
xmax = 1;
dx = 0.05;

% X changes with i, Y changes with j
[Y, X] = meshgrid(xmin:dx:xmax);
ngrid = max(size(X));
x = linspace(xmin, xmax, ngrid);
y = linspace(xmin, xmax, ngrid);

rho = 1;
Gamma = 50;
Q = 50;
vinf = 10;
a = 0.23;

for i = 1:ngrid
    for j = 1:ngrid
        % Find radius to each cartesian point
        r(i,j) = sqrt(x(i)^2+y(j)^2);
        rsorc(i,j) = sqrt(((a+x(i))^2)+(y(j)^2));
        rsink(i,j) = sqrt(((x(i)-a)^2)+(y(j)^2));
        rcw(i,j) = sqrt(((a+x(i))^2)+(y(j)^2));
        rccw(i,j) = sqrt(((x(i)-a)^2)+(y(j)^2));
        uinf(i,j) = vinf;
        cos(i,j) = x(i)/r(i,j);
        sin(i,j) = y(j)/r(i,j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates velocity vectors and pressure contours for a doublet with uniform flow
% around it: Rankine Body

% Calculate source and sink velocity magnitudes
for i = 1:ngrid
    for j = 1:ngrid
        vsorc(i,j) = Q/(2*pi*rho*rsorc(i,j));
        vsink(i,j) = -Q/(2*pi*rho*rsink(i,j));
    end
end

% Calculate source and sink velocity components
for i = 1:ngrid
    for j = 1:ngrid
        vsorcx(i,j) = vsorc(i,j)*((x(i)+a)/rsorc(i,j));
        vsorcy(i,j) = vsorc(i,j)*(y(j)/rsorc(i,j));
        vsinkx(i,j) = vsink(i,j)*((x(i)-a)/rsink(i,j));
        vsinky(i,j) = vsink(i,j)*(y(j)/rsink(i,j));
    end
end

% Add source and sink velocity components
for i = 1:ngrid
    for j = 1:ngrid
        vrankx(i,j) = vsorcx(i,j)+vsinkx(i,j)+uinf(i,j);
        vranky(i,j) = vsorcy(i,j)+vsinky(i,j);
        lim = 50;

        if vrankx(i,j)>lim
            vrankx(i,j) = lim;
        end

        if vrankx(i,j)<-lim
            vrankx(i,j) = -lim;
        end

        if vranky(i,j)>lim
            vranky(i,j) = lim;
        end

        if vranky(i,j)<-lim
            vranky(i,j) = -lim;
        end
    end
end

% Calculate rankine body cp and velocity magnitudes
for i = 1:ngrid
    for j = 1:ngrid
        vrank(i,j) = sqrt((vrankx(i,j)^2)+(vranky(i,j)^2));
        prank(i,j) = 1-((vrank(i,j)^2)/(vinf^2));
        lim = 10;

        if prank(i,j)>lim
            prank(i,j) = lim;
        end

        if prank(i,j)<-lim
            prank(i,j) = -lim;
        end
    end
end

prank = prank';
vrank = vrank';

figure(1)
subplot(2,2,1)
quiver(X,Y,vrankx,vranky,1)
axis square
axis([xmin,xmax,xmin,xmax])
colorbar
title('Rankine Body: Velocity Vectors')
subplot(2,2,2)
contour(y,x,vrank,20)
axis square
colorbar
title('Rankine Body: Velocity Contours')
subplot(2,2,3)
contour(y,x,prank,20)
axis square
colorbar
title('Rankine Body: Pressure Contours')
subplot(2,2,4)
contour(y,x,prank,20)
axis square
hold on
quiver(X,Y,vrankx,vranky,1)
colorbar
title('Rankine Body: Pressure Contours and Velocity Vectors')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates velocity vectors and pressure contours for a two vortices with uniform flow
% around it: Kelvin Oval

% Calculate cw and ccw vortex velocity magnitudes
for i = 1:ngrid
    for j = 1:ngrid
        vcw(i,j) = Gamma/(2*pi*rcw(i,j));
        vccw(i,j) = Gamma/(2*pi*rccw(i,j));
    end
end

% Calculate cw and ccw vortex velocity components
for i = 1:ngrid
    for j = 1:ngrid
        vcwx(i,j) = vcw(i,j)*(y(j)/rcw(i,j));
        vcwy(i,j) = -vcw(i,j)*((x(i)+a)/rcw(i,j));
        vccwx(i,j) = -vccw(i,j)*(y(j)/rccw(i,j));
        vccwy(i,j) = vccw(i,j)*((x(i)-a)/rccw(i,j));
    end
end

% Add cw and ccw vortex components
for i = 1:ngrid
    for j = 1:ngrid
        vkelvx(i,j) = vccwx(i,j)+vcwx(i,j);
        vkelvy(i,j) = vccwy(i,j)+vcwy(i,j)+uinf(i,j);
        lim = 50;

        if vkelvx(i,j)>lim
            vkelvx(i,j) = lim;
        end

        if vkelvx(i,j)<-lim
            vkelvx(i,j) = -lim;
        end

        if vkelvy(i,j)>lim
            vkelvy(i,j) = lim;
        end

        if vkelvy(i,j)<-lim
            vkelvy(i,j) = -lim;
        end
    end
end

% Calculate doublet velocity magnitudes and pressures
for i = 1:ngrid
    for j = 1:ngrid
        vkelv(i,j) = sqrt((vkelvx(i,j)^2)+(vkelvy(i,j)^2));
        pkelv(i,j) = 1-((vkelv(i,j)^2)/(vinf^2));
        lim = 10;

        if pkelv(i,j)>lim
            pkelv(i,j) = lim;
        end

        if pkelv(i,j)<-lim
            pkelv(i,j) = -lim;
        end
    end
end

pkelv = pkelv';
vkelv = vkelv';

figure(2)
subplot(2,2,1)
quiver(X,Y,vkelvx,vkelvy,1)
axis square
axis([xmin,xmax,xmin,xmax])
colorbar
title('Kelvin Oval: Velocity Vectors')

subplot(2,2,2)
contour(x,y,vkelv,20)
axis square
axis([xmin,xmax,xmin,xmax])
colorbar
title('Kelvin Oval: Velocity Contours')

subplot(2,2,3)
contour(x,y,pkelv,20)
axis square
axis([xmin,xmax,xmin,xmax])
colorbar
title('Kelvin Oval: Pressure Contours')

subplot(2,2,4)
contour(x,y,pkelv,20)
axis square
axis([xmin,xmax,xmin,xmax])
hold on
quiver(X,Y,vkelvx,vkelvy,1)
colorbar
title('Kelvin Oval: Pressure Contours and Velocity Vectors')
hold off
