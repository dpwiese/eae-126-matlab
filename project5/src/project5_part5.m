%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 5 - Part 1 - Problem 5: Oblique/Crescent Wings
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The wing has a sweep of 45 degrees, and constant chord

uinf = 1;
aoa_deg = 5;
Ma = 0.0;

% nx and ny must be even
nx = 100;
ny = 20;
nz = nx;

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
zmin = -10;
zmax = 10;

dx = (xmax-xmin)/(nx-1);
dy = (ymax-ymin)/(ny-1);
dz = (zmax-zmin)/(nz-1);

x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);
z = linspace(zmin,zmax,nz);

span = 10;
chord = 4;
area = span*chord;
AR = span^2/area;

for k = 1:nz/2
    [rspan] = find(z(k+1)>-(span/2) & z(k)<-span/2);
    rlens(k) = length(rspan);
end
    nLEFT = find(rlens);
    nRGHT = nz-nLEFT+1;

for i = 1:nx/2
    [rchord] = find(x(i+1)>-(chord/2) & x(i)<-chord/2);
    rlenc(i) = length(rchord);
end

nLE(nz/2) = find(rlenc)+1;
nLE(nz/2+1) = nLE(nz/2)-1;
nCHRD = round(chord/dx);
nTE(nz/2) = nLE(nz/2)+nCHRD;
nTE(nz/2+1) = nLE(nz/2+1)+nCHRD;

maxiter = 300;
omega = 1.95;
resmax = 10^-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:nLEFT-1
    nLE(k) = 1;
    nTE(k) = 1;
end

% Set up sweep on right side
for k = nz/2+2:nRGHT
    nLE(k) = nLE(nz/2+1)-k+nz/2+1;
    nTE(k) = nLE(k)+nCHRD;
end

% Set up sweep on left side
for k = nz/2-1:-1:nLEFT
    nLE(k) = nLE(nz/2)-k+nz/2;
    nTE(k) = nLE(k)+nCHRD;
end

for k = nRGHT+1:nz
    nLE(k) = 1;
    nTE(k) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:nz
    xLE(k) = x(nLE(k));
    xTE(k) = x(nTE(k));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flat wing
for k = 1:nz
    for i = 1:nx
        yBT(i,k) = 0;
        yBB(i,k) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:nz
    for i = 1:nLE(k)-1
        alpha(i,k) = 0;
    end
    for i = nLE(k):nTE(k)
        alpha(i,k) = aoa_deg*pi/180;
    end
    for i = nTE(k)+1:nx
        alpha(i,k) = 0;
    end
end
Gamma = 0;

for k = 1:nz
    for i = 1:nx
        for j = 1:ny
            Beta(i,j,k) = sqrt(1-Ma);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = zeros(nx,ny,nz);
uold = zeros(nx,ny,nz);
res = 1;
iter = 1;
while (iter<maxiter && res>resmax)
    % Plaens before left tip of wing
    for k = 2:nLEFT-1
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
    end

    % Planes across wing
    for k = nLEFT:nRGHT
        % Points upstream of airfoil
        for i = 2:nLE(k)-1
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
        for i = nLE(k):nTE(k)
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
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
            end

            % Just before body
            for j = ny/2;
                A(j) = 1/dy^2;
                B(j) = -1/dy^2-2/dx^2-2/dz^2;
                C(j) = 0;
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2-(yBB(i+1,k)-2*yBB(i,k)+yBB(i-1,k))/(dx^2*dy);
            end

            % Just after body
            for j = ny/2+1;
                A(j) = 0;
                B(j) = -1/dy^2-2/dx^2-2/dz^2;
                C(j) = 1/dy^2;
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2-(uold(i,j,k-1)+uold(i,j,k+1))/dz^2+(yBT(i+1,k)-2*yBT(i,k)+yBT(i-1,k))/(dx^2*dy);
            end

            % Region above body
            for j = ny/2+2:ny-1
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

        % Points downstream of airfoil
        for i = nTE(k)+1:nx-1
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
    end

    % Planes after right tip of wing
    for k = nRGHT+1:nz-1
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
    end

    % Residual
    for k = 2:nz-1
        for i = 2:nx-1
            for j = 2:ny-1
                resid(i,j,k) = abs((u(i+1,j,k)-2*u(i,j,k)+u(i-1,j,k))/dx^2+(u(i,j+1,k)-2*u(i,j,k)+u(i,j-1,k))/dy^2+(u(i,j,k+1)-2*u(i,j,k)+u(i,j,k-1))/dz^2);
            end
        end
    end

    for k = nLEFT:nRGHT
        for i = nLE(k):nTE(k)
            for j = ny/2:ny/2+1
                resid(i,j,k) = 0;
            end
        end
    end

    % % Calculate Gamma
    % uBot = 0;
    % uTop = 0;
    % for i = nLE:nTE
    %     % Just before body
    %     for j = ny/2;
    %         uBot = uBot+u(i,j);
    %     end
    %     % Just after body
    %     for j = ny/2+1;
    %         uTop = uTop+u(i,j);
    %     end
    % end
    % Gamma = dx*(uBot-uTop);

    res = max(max(max(resid)))
    resplot(iter) = res;
    iter = iter+1
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:nz
    for i = 1:nx
        cpbot(i,k) = -2*u(i,ny/2,k)/uinf;
        cptop(i,k) = -2*u(i,ny/2+1,k)/uinf;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(xLE,z)
title('Wing Planform')
hold on
plot(xTE,z)
grid on
hold off
axis([-10 10 -span/2 span/2])
daspect([1 1 1])
xlabel('x-axis')
ylabel('z-axis')

figure(2)
semilogy(resplot)
title('Residual versus Iteration Count')
xlabel('Number of Iterations')
ylabel('Residual')
grid on

figure(3)
surf(x,z,-cptop')
title('-C_P Across Top of Wing')

figure(4)
plot(x,alpha(:,nz/2))
title('Angle of Attack on Airfoil')
grid on
