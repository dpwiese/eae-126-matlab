%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 8/9 - BOUNDARY LAYER
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Re = 1000;

% ny must be even
nx = 20;
ny = 20;

xmin = -4;
xmax = 4;
ymin = 0;
ymax = 10/sqrt(Re);

dx = (xmax - xmin) / (nx - 1);
dy = (ymax - ymin) / (ny - 1);

x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);

for i = 1:nx
    for j = 1:ny
        X(i,j) = x(i);
        Y(i,j) = y(j);
    end
end

nLE = round(2*nx/5);
nTE = round(3*nx/5)+1;

% tau - thickness ratio of parabolic profiles and NACA
% sh - vertical shift for parabolic profiles points on grid
% coe - 'A' coefficient for parabolic profiles points on grid
% chord - chordlength of airfoil points on grid
tau = 0.1;
sh = tau*x(nTE);
coe = sh/x(nLE)^2;
chord = x(nTE)-x(nLE);

% For diamond airfoil
diamond_slope = tau;

itermax = 20;
iter2max = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case 1: flat plate
for i = 1:nx
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

u = ones(nx,ny);
v = ones(nx,ny);
w = ones(nx,ny);

v(1,:) = 0;
w(1,:) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for U
iter2 = 1;

while iter2<iter2max

    % Points before plate
    for i = 2:nLE-1

        iter = 1;

        while iter<itermax

            for j = 1;
                Au(j) = 0;
                Bu(j) = -2/dy^2-2/dx^2;
                Cu(j) = 2/dy^2;
                Du(j) = (-1/dx^2)*(u(i+1,j)+u(i-1,j));

                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = 0;
            end

            for j = 2:ny-1
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2-2/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = (-1/dx^2)*(u(i+1,j)+u(i-1,j)) - (1/(2*dy))*(w(i,j+1)-w(i,j-1));

                Av(j) = 1/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 1/dy^2;
                Dv(j) = (-1/dx^2)*(v(i+1,j)+v(i-1,j)) + (1/(2*dx))*(w(i+1,j)-w(i-1,j));
            end

            for j = ny
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = 1;

                Av(j) = 2/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 0;
                Dv(j) = (-1/dx^2)*(v(i+1,j)+v(i-1,j)) + (1/(2*dx))*(w(i+1,j)-w(i-1,j));
            end

            u(i,:) = tridiagscalar(Au, Bu, Cu, Du);
            v(i,:) = tridiagscalar(Av, Bv, Cv, Dv);

            for j = 1;
                Aw(j) = 0;
                Bw(j) = 1;
                Cw(j) = 0;
                Dw(j) = 0;
            end

            for j = 2:ny-1
                Aw(j) = -1/(Re*dy^2)-v(i,j)/(2*dy);
                Bw(j) = u(i,j)/dx+2/(Re*dx^2)+2/(Re*dy^2);
                Cw(j) = v(i,j)/(2*dy)-1/(Re*dy^2);
                Dw(j) = (w(i-1,j)*u(i,j))/dx + (w(i+1,j)+w(i-1,j))/(Re*dx^2);
            end

            for j = ny
                Aw(j) = 0;
                Bw(j) = 1;
                Cw(j) = 0;
                Dw(j) = 0;
            end
            w(i,:) = tridiagscalar(Aw, Bw, Cw, Dw);
            iter = iter+1;
        end
    end

    % Points across plate
    for i = nLE:nTE

        iter = 1;

        while iter<itermax

            for j = 1;
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = 0;

                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = 0;
            end

            for j = 2:ny-1
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2-2/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = (-1/dx^2)*(u(i+1,j)+u(i-1,j)) - (1/(2*dy))*(w(i,j+1)-w(i,j-1));

                Av(j) = 1/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 1/dy^2;
                Dv(j) = (-1/dx^2)*(v(i+1,j)+v(i-1,j)) + (1/(2*dx))*(w(i+1,j)-w(i-1,j));
            end

            for j = ny
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = 1;

                Av(j) = 2/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 0;
                Dv(j) = (-1/dx^2)*(v(i+1,j)+v(i-1,j)) + (1/(2*dx))*(w(i+1,j)-w(i-1,j));
            end

            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);
            v(i,:) = tridiagscalar(Av,Bv,Cv,Dv);

            % Boundary condition for 'w' across flate plate
            for j = 1;
                Aw(j) = 0;
                Bw(j) = 1;
                Cw(j) = 0;
                Dw(j) = -u(i,j+1)/dy;
            end

            for j = 2:ny-1
                Aw(j) = -1/(Re*dy^2)-v(i,j)/(2*dy);
                Bw(j) = u(i,j)/dx+2/(Re*dx^2)+2/(Re*dy^2);
                Cw(j) = v(i,j)/(2*dy)-1/(Re*dy^2);
                Dw(j) = (w(i-1,j)*u(i,j))/dx + (w(i+1,j)+w(i-1,j))/(Re*dx^2);
            end

            for j = ny
                Aw(j) = 0;
                Bw(j) = 1;
                Cw(j) = 0;
                Dw(j) = 0;
            end

            w(i,:) = tridiagscalar(Aw,Bw,Cw,Dw);
            iter = iter+1;
        end
    end

    % Trailing edge to end of grid
    for i = nTE+1:nx-1

        iter = 1;

        while iter<itermax

            for j = 1;
                Au(j) = 0;
                Bu(j) = -2/dy^2-2/dx^2;
                Cu(j) = 2/dy^2;
                Du(j) = (-1/dx^2)*(u(i+1,j)+u(i-1,j));

                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = 0;
            end

            for j = 2:ny-1
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2-2/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = (-1/dx^2)*(u(i+1,j)+u(i-1,j)) - (1/(2*dy))*(w(i,j+1)-w(i,j-1));

                Av(j) = 1/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 1/dy^2;
                Dv(j) = (-1/dx^2)*(v(i+1,j)+v(i-1,j)) + (1/(2*dx))*(w(i+1,j)-w(i-1,j));
            end

            for j = ny
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = 1;

                Av(j) = 2/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 0;
                Dv(j) = (-1/dx^2)*(v(i+1,j)+v(i-1,j)) + (1/(2*dx))*(w(i+1,j)-w(i-1,j));
            end

            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);
            v(i,:) = tridiagscalar(Av,Bv,Cv,Dv);

            for j = 1;
                Aw(j) = 0;
                Bw(j) = 1;
                Cw(j) = 0;
                Dw(j) = 0;
            end

            for j = 2:ny-1
                Aw(j) = -1/(Re*dy^2)-v(i,j)/(2*dy);
                Bw(j) = u(i,j)/dx+2/(Re*dx^2)+2/(Re*dy^2);
                Cw(j) = v(i,j)/(2*dy)-1/(Re*dy^2);
                Dw(j) = (w(i-1,j)*u(i,j))/dx+(w(i+1,j)+w(i-1,j))/(Re*dx^2);
            end

            for j = ny
                Aw(j) = 0;
                Bw(j) = 1;
                Cw(j) = 0;
                Dw(j) = 0;
            end

            w(i,:) = tridiagscalar(Aw,Bw,Cw,Dw);
            iter = iter+1;
        end
    end

    % End of grid
    for i = nx

        iter = 1;

        while iter<itermax

            for j = 1;
                Au(j) = 0;
                Bu(j) = -2/dy^2-2/dx^2;
                Cu(j) = 2/dy^2;
                Du(j) = (-1/dx^2)*(2*u(i-1,j));

                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = 0;
            end

            for j = 2:ny-1
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2-2/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = (-1/dx^2)*(2*u(i-1,j)) - (1/(2*dy))*(w(i,j+1)-w(i,j-1));

                Av(j) = 1/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 1/dy^2;
                Dv(j) = (-1/dx^2)*(2*v(i-1,j));
            end

            for j = ny
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = 1;

                Av(j) = 2/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 0;
                Dv(j) = (-1/dx^2)*(2*v(i-1,j));
            end

            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);
            v(i,:) = tridiagscalar(Av,Bv,Cv,Dv);

            for j = 1;
                Aw(j) = 0;
                Bw(j) = 1;
                Cw(j) = 0;
                Dw(j) = 0;
            end

            for j = 2:ny-1
                Aw(j) = -1/(Re*dy^2)-v(i,j)/(2*dy);
                Bw(j) = u(i,j)/dx+2/(Re*dx^2)+2/(Re*dy^2);
                Cw(j) = v(i,j)/(2*dy)-1/(Re*dy^2);
                Dw(j) = (w(i-1,j)*u(i,j))/dx+(2*w(i-1,j))/(Re*dx^2);
            end

            for j = ny
                Aw(j) = 0;
                Bw(j) = 1;
                Cw(j) = 0;
                Dw(j) = 0;
            end

            w(i,:) = tridiagscalar(Aw,Bw,Cw,Dw);
            iter = iter+1;
        end
    end
    iter2 = iter2+1
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
contour(x,y,u',200)
hold on
plot(x,yB,'-k')
fill(x,yB,'k')
title('u Countours Around Airfoil')
xlabel('u-direction')
ylabel('v-direction')
colorbar
grid on

figure(2)
contour(x,y,v',200)
hold on
plot(x,yB,'-k')
fill(x,yB,'k')
title('v Countours Around Airfoil')
xlabel('u-direction')
ylabel('v-direction')
colorbar
grid on

figure(3)
contour(x,y,w',200)
hold on
plot(x,yB,'-k')
fill(x,yB,'k')
title('w Countours Around Airfoil')
xlabel('u-direction')
ylabel('v-direction')
colorbar
grid on
