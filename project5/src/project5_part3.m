%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 5 - THICKNESS - Problem 3: Rectangular Wing
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uinf  =  1;
Ma  =  0.5;

% Grid - ny must be even number
nx  =  10;
ny  =  10;
nz  =  10;

xmin  =  -10;
xmax  =  10;
ymin  =  -10;
ymax  =  10;
zmin  =  -10;
zmax  =  10;

dx  =  (xmax-xmin)/(nx-1);
dy  =  (ymax-ymin)/(ny-1);
dz  =  (zmax-zmin)/(nz-1);

x  =  linspace(xmin, xmax, nx);
y  =  linspace(ymin, ymax, ny);
z  =  linspace(zmin, zmax, nz);

nLEFT  =  round(2*nz/5);
nRGHT  =  round(3*nz/5)+1;

for k = 1:nLEFT-1
    nLE(k) = 1;
    nTE(k) = 1;
end

for k = nLEFT:nRGHT
    nLE(k) = round(2*nx/5);
    nTE(k) = round(3*nx/5)+1;
end

for k = nRGHT+1:nz
    nLE(k) = 1;
    nTE(k) = 1;
end

maxiter = 300;
omega = 1.95;
resmax = 10^-6;

% tau - thickness ratio of parabolic profiles and NACA
% sh - vertical shift for parabolic profiles
% coe - 'A' coefficient for parabolic profiles
% chord - chordlength of airfoil points on grid
tau = 0.06;
sh = tau*(x(nTE(nz/2))+dx/2);
coe2 = sh/x(nLE(nz/2))^2;
chord = (x(nTE(nz/2))-x(nLE(nz/2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:nLEFT-1
    for i = 1:nx
        yBT(i,k) = 0;
        yBB(i,k) = 0;
    end
end

% Case 1: Biconvex Rectangular Wing
for k = nLEFT:nRGHT
    for i = 1:nLE(k)-1
        yBT(i,k) = 0;
        yBB(i,k) = 0;
    end
    for i = nLE(k):nTE(k)
        yBT(i,k) = -coe2*x(i)^2+sh;
        yBB(i,k) = coe2*x(i)^2-sh;
    end
    for i = nTE(k)+1:nx
        yBT(i,k) = 0;
        yBB(i,k) = 0;
    end
end

% % Case 2: NACA Rectangular Wing
% for k = nLEFT:nRGHT
%     for i = 1:nLE(k)-1
%         yBT(i,k) = 0;
%         yBB(i,k) = 0;
%     end

%     % NACA 00xx
%     for i = nLE(k):nTE(k)
%         yBT(i,k) = 10*tau*chord*(0.2969*sqrt((x(i)+chord/2)/chord) - ...
%                     0.1260*((x(i)+chord/2)/chord) - ...
%                     0.3537*((x(i)+chord/2)/chord)^2 + ...
%                     0.2843*((x(i)+chord/2)/chord)^3 - ...
%                     0.1015*((x(i)+chord/2)/chord)^4);
%         yBB(i,k) = -yBT(i,k);
%     end

%     for i = nTE(k)+1:nx
%         yBT(i,k) = 0;
%         yBB(i,k) = 0;
%     end
% end

for k = nRGHT+1:nz
    for i = 1:nx
        yBT(i,k) = 0;
        yBB(i,k) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:nz
    for i = 1:nx
        for j = 1:ny
            Beta(i,j,k) = sqrt(1-Ma^2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = zeros(nx,ny,nz);
uold = zeros(nx,ny,nz);

res = 1;
iter = 1;

while (iter<maxiter && res>resmax)

    % Planes before left wingtip
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
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2 - ...
                        (uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
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
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2 - ...
                        (uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
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
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2 - ...
                        (uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
            end

            % Just before body
            for j = ny/2;
                A(j) = 1/dy^2;
                B(j) = -1/dy^2-2/dx^2-2/dz^2;
                C(j) = 0;
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2 - ...
                        (uold(i,j,k-1)+uold(i,j,k+1))/dz^2 - ...
                        (yBB(i+1,k)-2*yBB(i,k)+yBB(i-1,k))/(dx^2*dy);
            end

            % Just after body
            for j = ny/2+1;
                A(j) = 0;
                B(j) = -1/dy^2-2/dx^2-2/dz^2;
                C(j) = 1/dy^2;
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2 - ...
                        (uold(i,j,k-1)+uold(i,j,k+1))/dz^2 + ...
                        (yBT(i+1,k)-2*yBT(i,k)+yBT(i-1,k))/(dx^2*dy);
            end

            % Region above body
            for j = ny/2+2:ny-1
                A(j) = 1/dy^2;
                B(j) = -2/dy^2-2/dx^2-2/dz^2;
                C(j) = 1/dy^2;
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2 - ...
                        (uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
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
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2 - ...
                        (uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
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

    % Planes after right wingtip
    for k = nRGHT + 1 : nz - 1
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
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2 - ...
                        (uold(i,j,k-1)+uold(i,j,k+1))/dz^2;
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
                resid(i,j,k) = abs( ...
                    (u(i+1,j,k)-2*u(i,j,k)+u(i-1,j,k))/dx^2 + ...
                    (u(i,j+1,k)-2*u(i,j,k)+u(i,j-1,k))/dy^2 + ...
                    (u(i,j,k+1)-2*u(i,j,k)+u(i,j,k-1))/dz^2);
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

    res = max(max(max(resid)));
    resplot(iter) = res;
    iter = iter + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for k = nLEFT:nRGHT
%     for i = nLE:nTE
%         xcp(i-nLE+1) = x(i);
%         cpbot(i-nLE+1) = -2*u(i,ny/2)/uinf;
%         cptop(i-nLE+1) = -2*u(i,ny/2+1)/uinf;
%     end
% end

for k = 1:nz
    for i = 1:nx
        cpbot(i,k) = -2*u(i,ny/2,k)/uinf;
        cptop(i,k) = -2*u(i,ny/2+1,k)/uinf;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)

subplot(2, 2, 1)
plot(x,yBB(:,nz/2),'--b')
title('Airfoil Profile')
hold on
plot(x,yBT(:,nz/2),'-b')
daspect([1 1 1])
axis([xmin xmax -4 4])
legend('Bottom of Profile','Top of Profile')
grid on
hold off

subplot(2, 2, 2)
surf(x,z,-cptop')
title('-C_P Across Top of Wing')

subplot(2, 2, 3)
semilogy(resplot)
title('Residual versus Iteration Count')
xlabel('Number of Iterations')
ylabel('Residual')
grid on
