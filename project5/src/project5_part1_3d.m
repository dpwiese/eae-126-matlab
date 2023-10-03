%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 5 - Part 1 - Problem 1
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uinf = 1;
aoa_deg = 0;
Ma = 0.0;

% ny must be even number
nx = 40;
ny = 40;
nz = 40;

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
zmin = -10;
zmax = 10;

dx = (xmax-xmin)/(nx-1);
dy = (ymax-ymin)/(ny-1);
dz = (zmax-zmin)/(nz-1);

x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);
z = linspace(zmin, zmax, nz);

nLE = round(2*nx/5);
nTE = round(3*nx/5)+1;
nLEFT = round(2*nz/5);
nRGHT = round(3*nz/5)+1;

maxiter = 100;
omega = 1.95;
resmax = 10^-6;

% tau - thickness ratio of parabolic profiles and NACA
% sh - vertical shift for parabolic profiles
% coe - 'A' coefficient for parabolic profiles
% chord - chordlength of airfoil points on grid
tau = 0.1;
sh = tau*(x(nTE)+dx/2);
coe2 = sh/x(nLE)^2;
chord = (x(nTE)-x(nLE));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Joukowski Stuff
ntheta = 1000;
theta = linspace(0,2*pi,ntheta);
epsilon = 0.1;
mu = 0.1;

% Adjust size of 'b' so chord length is correct
b = (chord/4)*0.99;
a = sqrt(mu^2+(b-epsilon)^2);

for i = 1:ntheta
    xj(i) = epsilon+a*cos(theta(i));
    yj(i) = mu+a*sin(theta(i));
    X(i) = xj(i)*(1+(b^2)/(xj(i)^2+yj(i)^2));
    Y(i) = yj(i)*(1-(b^2)/(xj(i)^2+yj(i)^2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nLE-1
    yBT(i) = 0;
    yBB(i) = 0;
end

% % Flat Plate
% for i = nLE:nTE
%     yBT(i) = 0;
%     yBB(i) = 0;
% end

% % Cambered Plate points on grid
% for i = nLE:nTE
%     yBT(i) = -coe2*x(i)^2+sh;
%     yBB(i) = -coe2*x(i)^2+sh;
% end

% Biconvex points on grid
for i = nLE:nTE
    yBT(i) = -coe2*x(i)^2+sh;
    yBB(i) = coe2*x(i)^2-sh;
end

% % NACA 00xx
% for i = nLE:nTE
%     yBT(i) = 10*tau*chord*(0.2969*sqrt((x(i)+chord/2)/chord)-0.1260*((x(i)+chord/2)/chord)-0.3537*((x(i)+chord/2)/chord)^2+0.2843*((x(i)+chord/2)/chord)^3-0.1015*((x(i)+chord/2)/chord)^4);
%     yBB(i) = -yBT(i);
% end

% % Joukowski
% for i = nLE:nTE
%     [r] = find(X>x(i-1) & X<x(i+1)); %find location in X that I want for each x
%     rlen(i) = length(r);
%     yBB(i) = Y(r(round(0.15*rlen(i))));
%     yBT(i) = Y(r(round(0.85*rlen(i))));
%     % yBB(i) = Y(r(1));
%     % yBT(i) = Y(r(rlen(i)));
% end

for i = nTE+1:nx
    yBT(i) = 0;
    yBB(i) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1
    d2yBBdx2(i) = 0;
end

for i = 2:nx-1
    d2yBBdx2(i) = (yBB(i+1)-2*yBB(i)+yBB(i-1))/dx^2;
end

for i = nx
    d2yBBdx2(i) = 0;
end

for i = 1
    d2yBTdx2(i) = 0;
end

for i = 2:nx-1
    d2yBTdx2(i) = (yBT(i+1)-2*yBT(i)+yBT(i-1))/dx^2;
end

for i = nx
    d2yBTdx2(i) = 0;
end

for i = 1:nLE-1
    alpha(i) = 0;
end

for i = nLE:nTE
    alpha(i) = aoa_deg*pi/180;
end

for i = nTE+1:nx
    alpha(i) = 0;
end

Gamma = 0;

% for i = 1:nx
%     for j = 1:ny
%         Beta(i,j) = sqrt(1-Ma);
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = zeros(nx,ny,nz);
uold = zeros(nx,ny,nz);

% Planes before left tip of wing
for k = 2:nLEFT-1
    res = 1;
    iter = 1;
    while (iter<maxiter && res>resmax)
        for i = 2:nx-1
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end
            for j = 2:ny-1
                A(j) = 1/dy^2;
                B(j) = -2/dy^2-2/dx^2-2/dz^2;
                C(j) = 1/dy^2;
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
            end
            for j = ny
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end
            u(i,:,k) = tridiagscalar(A,B,C,D);
            u(i,:,k) = uold(i,:,k)+omega*(u(i,:,k)-uold(i,:,k));
            uold(i,:,k) = u(i,:,k);
        end
        for i = 2:nx-1
            for j = 2:ny-1
                resid(i,j,k) = abs((u(i+1,j,k)-2*u(i,j,k)+u(i-1,j,k))/dx^2+(u(i,j+1,k)-2*u(i,j,k)+u(i,j-1,k))/dy^2+(u(i,j,k+1)-2*u(i,j,k)+u(i,j,k-1))/dz^2);
            end
        end
        for i = nLE-1:nTE+1
            for j = ny/2-1:ny/2+2
                resid(i,j,k) = 0;
            end
        end
            res = max(max(max(resid)));
            resplot(iter) = res;
            iter = iter+1;
    end
end

% Planes across wing
for k = nLEFT:nRGHT
    res = 1;
    iter = 1;
    while (iter<maxiter && res>resmax)

        % Points upstream of airfoil
        for i = 2:nLE-1
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end
            for j = 2:ny-1
                A(j) = 1/dy^2;
                B(j) = -2/dy^2-2/dx^2-2/dz^2;
                C(j) = 1/dy^2;
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
            end
            for j = ny
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end
            u(i,:,k) = tridiagscalar(A,B,C,D);
            u(i,:,k) = uold(i,:,k)+omega*(u(i,:,k)-uold(i,:,k));
            uold(i,:,k) = u(i,:,k);
        end

        % Points where airfoil is
        for i = nLE:nTE-1
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end

            % Region below body
            for j = 2:ny/2-1;
                A(j) = 1/dy^2;
                B(j) = -2/dy^2-2/dx^2-2/dz^2;
                C(j) = 1/dy^2;
                D(j) = -(uold(i-1,j)+uold(i+1,j))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
            end

            % Just before body
            for j = ny/2;
                A(j) = 1/dy^2;
                B(j) = -1/dy^2-2/dx^2-2/dz^2;
                C(j) = 0;
                D(j) = -(uold(i-1,j)+uold(i+1,j))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2-(yBB(i+1)-2*yBB(i)+yBB(i-1))/(dx^2*dy);
            end

            % Just after body
            for j = ny/2+1;
                A(j) = 0;
                B(j) = -1/dy^2-2/dx^2-2/dz^2;
                C(j) = 1/dy^2;
                D(j) = -(uold(i-1,j)+uold(i+1,j))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2+(yBB(i+1)-2*yBB(i)+yBB(i-1))/(dx^2*dy);
            end

            % Region above body
            for j = ny/2+2:ny-1
                A(j) = 1/dy^2;
                B(j) = -2/dy^2-2/dx^2-2/dz^2;
                C(j) = 1/dy^2;
                D(j) = -(uold(i-1,j)+uold(i+1,j))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
            end
            for j = ny
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end
            u(i,:,k) = tridiagscalar(A,B,C,D);
            u(i,:,k) = uold(i,:,k)+omega*(u(i,:,k)-uold(i,:,k));
            uold(i,:,k) = u(i,:,k);
        end

        % Trailing Edge
        for i = nTE
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end

            % Region below body
            for j = 2:ny/2-1;
                A(j) = 1/dy^2;
                B(j) = -2/dy^2-2/dx^2-2/dz^2;
                C(j) = 1/dy^2;
                D(j) = -(uold(i-1,j)+uold(i+1,j))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
            end

            % Just before body
            for j = ny/2;
                A(j) = 1/dy^2;
                B(j) = -1/dy^2-2/dx^2-2/dz^2;
                C(j) = 0;
                D(j) = -(uold(i-1,j)+uold(i+1,j))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2-(yBB(i+1)-2*yBB(i)+yBB(i-1))/(dx^2*dy);
            end

            % Just after body
            for j = ny/2+1;
                A(j) = 0;
                B(j) = -1/dy^2-2/dx^2-2/dz^2;
                C(j) = 1/dy^2;
                D(j) = -(uold(i-1,j)+uold(i+1,j))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2+(yBB(i+1)-2*yBB(i)+yBB(i-1))/(dx^2*dy);
            end

            % Region above body
            for j = ny/2+2:ny-1
                A(j) = 1/dy^2;
                B(j) = -2/dy^2-2/dx^2-2/dz^2;
                C(j) = 1/dy^2;
                D(j) = -(uold(i-1,j)+uold(i+1,j))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
            end
            for j = ny
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end
            u(i,:,k) = tridiagscalar(A,B,C,D);
            u(i,:,k) = uold(i,:,k)+omega*(u(i,:,k)-uold(i,:,k));
            uold(i,:,k) = u(i,:,k);
        end

        % Points downstream of airfoil
        for i = nTE+1:nx-1
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end
            for j = 2:ny-1
                A(j) = 1/dy^2;
                B(j) = -2/dy^2-2/dx^2-2/dz^2;
                C(j) = 1/dy^2;
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
            end
            for j = ny
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end
            u(i,:,k) = tridiagscalar(A,B,C,D);
            u(i,:,k) = uold(i,:,k)+omega*(u(i,:,k)-uold(i,:,k));
            uold(i,:,k) = u(i,:,k);
        end

        for i = 2:nx-1
            for j = 2:ny-1
                resid(i,j,k) = abs((u(i+1,j,k)-2*u(i,j,k)+u(i-1,j,k))/dx^2+(u(i,j+1,k)-2*u(i,j,k)+u(i,j-1,k))/dy^2+(u(i,j,k+1)-2*u(i,j,k)+u(i,j,k-1))/dz^2);
            end
        end

        for i = nLE-1:nTE+1
            for j = ny/2-1:ny/2+2
                resid(i,j,k) = 0;
            end
        end

        res = max(max(max(resid)));
        resplot(iter) = res;
        iter  =  iter + 1;
    end
end

% Planes after right tip of wing
for k = nRGHT+1:nz-1
    res = 1;
    iter = 1;
    while (iter<maxiter && res>resmax)
        for i = 2:nx-1
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end
            for j = 2:ny-1
                A(j) = 1/dy^2;
                B(j) = -2/dy^2-2/dx^2-2/dz^2;
                C(j) = 1/dy^2;
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
            end
            for j = ny
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end
            u(i,:,k) = tridiagscalar(A,B,C,D);
            u(i,:,k) = uold(i,:,k)+omega*(u(i,:,k)-uold(i,:,k));
            uold(i,:,k) = u(i,:,k);
        end
        for i = 2:nx-1
            for j = 2:ny-1
                resid(i,j,k) = abs((u(i+1,j,k)-2*u(i,j,k)+u(i-1,j,k))/dx^2+(u(i,j+1,k)-2*u(i,j,k)+u(i,j-1,k))/dy^2+(u(i,j,k+1)-2*u(i,j,k)+u(i,j,k-1))/dz^2);
            end
        end
        for i = nLE-1:nTE+1
            for j = ny/2-1:ny/2+2
                resid(i,j,k) = 0;
            end
        end
            res = max(max(max(resid)));
            resplot(iter) = res;
            iter = iter+1;
    end
end

% Residual

% % Calculate Gamma
%     uBot = 0;
%     uTop = 0;
% for i = nLE:nTE
%     for j = ny/2; %just before body
%         uBot = uBot+u(i,j);
%     end
%     for j = ny/2+1; %just after body
%         uTop = uTop+u(i,j);
%     end
% end
%     Gamma = dx*(uBot-uTop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = nLE:nTE
%     xcp(i-nLE+1) = x(i);
%     cpbot(i-nLE+1) = -2*u(i,ny/2)/uinf;
%     cptop(i-nLE+1) = -2*u(i,ny/2+1)/uinf;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(x,yBB,'--b')
title('Airfoil Profile')
hold on
plot(x,yBT,'-b')
daspect([1 1 1])
axis([xmin xmax -4 4])
legend('Bottom of Profile','Top of Profile')
grid on
hold off

figure(2)
semilogy(resplot)
title('Residual versus Iteration Count')
xlabel('Number of Iterations')
ylabel('Residual')
grid on

figure(3)
contour(x,y,u(:,:,20)')
title('Pressure Contours Around Airfoil')
xlabel('u-direction')
ylabel('v-direction')
colorbar

% for i = 1:nz
%     figure(i+2)
%     contour(x,y,u(:,:,i)')
%     title('Pressure Contours Around Airfoil')
%     xlabel('u-direction')
%     ylabel('v-direction')
%     colorbar
% end

% contour(x,y,u(:,:,round(nz/2))')
% title('Pressure Contours Around Airfoil')
% xlabel('u-direction')
% ylabel('v-direction')
% colorbar

% figure(4)
% plot(xcp,-cptop,'-b')
% title('-C_P Across Airfoil')
% xlabel('Location Along Chord')
% ylabel('-C_P')
% hold on
% plot(xcp,-cpbot,'--b')
% legend('C_P Top','C_P Bottom')
% grid on
% hold off

% figure(5)
% plot(x,alpha)
% title('Angle of Attack on Airfoil')
% grid on

% figure(6)
% plot(x,u(:,ny/2),'--b')
% hold on
% plot(x,u(:,ny/2+1),'-b')
% title('Perturbation Velocity Above and Below Airfoil')
% legend('Perturbation Velocity on Bottom','Perturbation Velocity on Top')
% grid on
% hold off

% figure(7)
% plot(x,d2yBBdx2,'--b')
% hold on
% plot(x,d2yBTdx2,'-b')
% grid on
% title('Second Derivative of Airfoil')
% legend('Second Derivative on Bottom','Second Derivative on Top')
% hold off
