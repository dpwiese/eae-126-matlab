%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 1 - Problem 2 - Flow Field: Source, Sink, Doublet
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
[Y, X]  =  meshgrid(xmin:dx:xmax);
ngrid = max(size(X));
x = linspace(xmin, xmax, ngrid);
y = linspace(xmin, xmax, ngrid);

rho = 1;
Q = 4;
Gamma = 4;
vinf = 1;
a = 0.33;
C = 0.2;

for i = 1:ngrid
    for j = 1:ngrid
        % Find radius to each cartesian point
        r(i,j) = sqrt(x(i)^2+y(j)^2);
        cos(i,j) = x(i)/r(i,j);
        sin(i,j) = y(j)/r(i,j);
    end
end

sin(round(ngrid/2),round(ngrid/2)) = 0;
cos(round(ngrid/2),round(ngrid/2)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates velocity vectors and pressure contours for a source

for i = 1:ngrid
    for j = 1:ngrid
        vsorc(i,j) = Q/(2*pi*rho*r(i,j));
        cpsorc(i,j) = 1-((vsorc(i,j)^2)/(vinf^2));
        lim = 10;

        if cpsorc(i,j)>lim
            cpsorc(i,j) = lim;
        end

        if cpsorc(i,j)<-lim
            cpsorc(i,j) = -lim;
        end

        if vsorc(i,j)>lim
            vsorc(i,j) = lim;
        end

        if vsorc(i,j)<-lim
            vsorc(i,j) = -lim;
        end
    end
end

for i = 1:ngrid
    for j = 1:ngrid
        vsorcx(i,j) = vsorc(i,j)*(x(i)/r(i,j));
        vsorcy(i,j) = vsorc(i,j)*(y(j)/r(i,j));
    end
end

vsorcx(round(ngrid/2),round(ngrid/2)) = 0;
vsorcy(round(ngrid/2),round(ngrid/2)) = 0;

figure(1)
subplot(2,2,1)
quiver(X,Y,vsorcx,vsorcy,3)
axis square
axis([-1,1,-1,1])
colorbar
title('Source: Velocity Vectors')

subplot(2,2,2)
contour(x,y,vsorc,40)
axis square
colorbar
title('Source: Velocity Contours')

subplot(2,2,3)
contour(x,y,cpsorc,40)
axis square
colorbar
title('Source: Pressure Contours')

subplot(2,2,4)
contour(x,y,cpsorc,40)
hold on
quiver(X,Y,vsorcx,vsorcy,3)
axis square
colorbar
title('Source: Pressure Contours and Velocity Vectors')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates velocity vectors and pressure contours for a vortex

for i = 1:ngrid
    for j = 1:ngrid
        vvort(i,j) = Gamma/(2*pi*r(i,j));
        cpvort(i,j) = 1-((vvort(i,j)^2)/(vinf^2));
        lim = 10;

        if cpvort(i,j)>lim
            cpvort(i,j) = lim;
        end

        if cpvort(i,j)<-lim
            cpvort(i,j) = -lim;
        end

        if vvort(i,j)>lim
            vvort(i,j) = lim;
        end

        if vvort(i,j)<-lim
            vvort(i,j) = -lim;
        end
    end
end

for i = 1:ngrid
    for j = 1:ngrid
        vvortx(i,j) = -vvort(i,j)*(y(j)/r(i,j));
        vvorty(i,j) = vvort(i,j)*(x(i)/r(i,j));
    end
end

vvortx(round(ngrid/2),round(ngrid/2)) = 0;
vvorty(round(ngrid/2),round(ngrid/2)) = 0;

figure(2)
subplot(2,2,1)
quiver(X,Y,vvortx,vvorty,3)
axis square
axis([-1,1,-1,1])
colorbar
title('Vortex: Velocity Vectors')

subplot(2,2,2)
contour(x,y,vvort,40)
axis square
colorbar
title('Vortex: Velocity Contours')
subplot(2,2,3)
contour(x,y,cpvort,40)
axis square
colorbar
title('Vortex: Pressure Contours')

subplot(2,2,4)
contour(x,y,cpvort,40)
hold on
quiver(X,Y,vvortx,vvorty,3)
axis square
colorbar
title('Vortex: Pressure Contours and Velocity Vectors')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates velocity vectors and pressure contours for a doublet

for i = 1:ngrid
    for j = 1:ngrid
        vdubx(i,j) = -2*C*(cos(i,j)^2)/(r(i,j)^2);
        vduby(i,j) = -2*C*(cos(i,j)*sin(i,j))/(r(i,j)^2);
        lim = 5;

        if vdubx(i,j)>lim
            vdubx(i,j) = lim;
        end

        if vdubx(i,j)<-lim
            vdubx(i,j) = -lim;
        end

        if vduby(i,j)>lim
            vduby(i,j) = lim;
        end

        if vduby(i,j)<-lim
            vduby(i,j) = -lim;
        end
    end
end

for i = round(ngrid/2)
    for j = round(ngrid/2)
        vdubx(i,j) = vdubx(i+1,j);
        vduby(i,j) = vduby(i+1,j);
    end
end

% Calculate velocity magnitude and pressure at each point
for i = 1:ngrid
    for j = 1:ngrid
        vdub(i,j) = sqrt((vdubx(i,j)^2)+(vduby(i,j)^2));
        cpdub(i,j) = 1-((vdub(i,j)^2)/(vinf^2));
        lim = 10;

        if cpdub(i,j)>lim
            cpdub(i,j) = lim;
        end

        if cpdub(i,j)<-lim
            cpdub(i,j) = -lim;
        end
    end
end

cpdub = cpdub';
vdub = vdub';

figure(3)
subplot(2,2,1)
quiver(X,Y,vdubx,vduby,3)
axis square
axis([-1,1,-1,1])
colorbar
title('Doublet: Velocity Vectors')

subplot(2,2,2)
contour(x,y,vdub,40)
axis square
colorbar
title('Doublet: Velocity Contours')

subplot(2,2,3)
contour(x,y,cpdub,40)
axis square
colorbar
title('Doublet: Pressure Contours')

subplot(2,2,4)
contour(x,y,cpdub,40)
hold on
quiver(X,Y,vdubx,vduby,3)
colorbar
axis square
title('Doublet: Pressure Contours and Velocity Vectors')
hold off
