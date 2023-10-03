%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 8/9 - Transonic Flow - SUPERSONIC
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uinf = 1;
rhoinf = 1;
gamma = 1.4;

% ny must be even
nx = 100;
ny = 100;

xmin = -10;
xmax = 10;
ymin = 0;
ymax = 20;

dx = (xmax-xmin)/(nx-1);
dy = (ymax-ymin)/(ny-1);

x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);

nLE = round(2*nx/5);
nTE = round(3*nx/5)+1;

% tau - thickness ratio of parabolic profiles and NACA
% sh - vertical shift for parabolic profiles points on grid
% coe - 'A' coefficient for parabolic profiles points on grid
% chord - chordlength of airfoil points on grid
tau = 0.1;
sh = tau*x(nTE);
coe = sh/x(nLE)^2;
chord = (x(nTE)-x(nLE));

% For diamond airfoil
diamond_slope = tau;

itermax = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nLE-1
    yB(i) = 0;
end

% Biconvex
% for i = nLE:nTE
%     yB(i) = -coe*x(i)^2+sh;
% end

% Diamond Airfoil (first half)
for i = nLE:nx/2
    yB(i) = diamond_slope*x(i)+tau*chord/2;
end

% Diamond Airfoil (second half)
for i = nx/2+1:nTE
    yB(i) = -diamond_slope*x(i)+tau*chord/2;
end

% NACA 00xx
% for i = nLE:nTE
%     yB(i) = 10*tau*chord*(0.2969*sqrt((x(i)+chord/2)/chord)-0.1260*((x(i)+chord/2)/chord)-0.3537*((x(i)+chord/2)/chord)^2+0.2843*((x(i)+chord/2)/chord)^3-0.1015*((x(i)+chord/2)/chord)^4);
% end

for i = nTE+1:nx
    yB(i) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nLE-1
    d2yBdx2(i) = 0;
end

for i = nLE:nTE
    d2yBdx2(i) = (yB(i+1)-2*yB(i)+yB(i-1))/dx^2;
end

for i = nTE+1:nx
    d2yBdx2(i) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Minf = 1.08;

u = ones(nx,ny);

u(1,:) = 0;
u(2,:) = 0;
u(nx,:) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solve for U

for i = 3:nx-1
    iter = 1;
    while iter<itermax
        for j = 1;
            fm2 = (1-Minf^2)*u(i-2,j)-((gamma+1)*Minf^2*u(i-2,j)^2)/2;
            fm1 = (1-Minf^2)*u(i-1,j)-((gamma+1)*Minf^2*u(i-1,j)^2)/2;
            fm0 = (1-Minf^2)*u(i,j)-((gamma+1)*Minf^2*u(i,j)^2)/2;

            TKm2 = 1-Minf^2-(gamma+1)*Minf^2*u(i-2,j);
            TKm1 = 1-Minf^2-(gamma+1)*Minf^2*u(i-1,j);
            TKm0 = 1-Minf^2-(gamma+1)*Minf^2*u(i,j);

            Au(j) = 0;
            Bu(j) = -1/dy^2+TKm0/dx^2;
            Cu(j) = 1/dy^2;
            Du(j) = (-1/dx^2)*(-2*TKm1*u(i-1,j)+TKm2*u(i-2,j)) - ...
                    (fm0-2*fm1+fm2)/dx^2 + ...
                    (TKm0*u(i,j)-2*TKm1*u(i-1,j)+TKm2*u(i-2,j))/dx^2 + ...
                    d2yBdx2(i)/dy;
        end
        for j = 2:ny-1
            fm2 = (1-Minf^2)*u(i-2,j)-((gamma+1)*Minf^2*u(i-2,j)^2)/2;
            fm1 = (1-Minf^2)*u(i-1,j)-((gamma+1)*Minf^2*u(i-1,j)^2)/2;
            fm0 = (1-Minf^2)*u(i,j)-((gamma+1)*Minf^2*u(i,j)^2)/2;

            TKm2 = 1-Minf^2-(gamma+1)*Minf^2*u(i-2,j);
            TKm1 = 1-Minf^2-(gamma+1)*Minf^2*u(i-1,j);
            TKm0 = 1-Minf^2-(gamma+1)*Minf^2*u(i,j);

            Au(j) = 1/dy^2;
            Bu(j) = -2/dy^2+TKm0/dx^2;
            Cu(j) = 1/dy^2;
            Du(j) = (-1/dx^2)*(-2*TKm1*u(i-1,j)+TKm2*u(i-2,j)) - ...
                    (fm0-2*fm1+fm2)/dx^2 + ...
                    (TKm0*u(i,j)-2*TKm1*u(i-1,j)+TKm2*u(i-2,j))/dx^2;
        end
        for j = ny
            Au(j) = 0;
            Bu(j) = 1;
            Cu(j) = 0;
            Du(j) = 0;
        end
        u(i,:) = tridiagscalar(Au,Bu,Cu,Du);
        iter = iter+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
contour(x,y,u',200)
hold on
plot(x,yB,'-k')
fill(x,yB,'k')
title('u Countours Around Airfoil')
xlabel('u-direction')
ylabel('v-direction')
colorbar
grid on
