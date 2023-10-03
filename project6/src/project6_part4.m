%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 5 - Part 4
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uinf = 1;
aoa_deg = 5;
alpha = aoa_deg*pi/180;
Ma = 0.0;
% AR = 6;

% ny must be even
nx = 30;
ny = 20;
nz = 100;

% xmin = -nx/2;
% xmax = nx/2;
xmin = -20;
xmax = 10;
ymin = -ny/2;
ymax = ny/2;
zmin = -nz/2;
zmax = nz/2;

dx = (xmax-xmin)/(nx-1);
dy = (ymax-ymin)/(ny-1);
dz = (zmax-zmin)/(nz-1);

x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);
z = linspace(zmin,zmax,nz);

% nLEFT = round(2*nz/5);
% nRGHT = round(3*nz/5)+1;
nLEFT = 22;
nRGHT = 78;

for k = 1:nLEFT-1
    nLE(k) = 1;
    nTE(k) = 1;
end

% for k = nLEFT:nRGHT
%     nLE(k) = round(2*nx/5);
%     nTE(k) = round(3*nx/5)+1;
% end

for k = nLEFT:nRGHT
    nLE(k) = 11;
    nTE(k) = 20;
end

for k = nRGHT+1:nz
    nLE(k) = 1;
    nTE(k) = 1;
end

maxiter = 300;
omega = 1.95;
resmax = 10^-6;

% Define chordlength, span, and aspect ratio
chord = (x(nTE(nz/2))-x(nLE(nz/2)));
span = z(nRGHT)-z(nLEFT);
AR = span/chord;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:nLEFT-1
    for i = 1:nx
        yBT(i,k) = 0;
        yBB(i,k) = 0;
    end
end

%vFlat Wing
for k = nLEFT:nRGHT
    for i = 1:nLE(k)-1
        yBT(i,k) = 0;
        yBB(i,k) = 0;
    end

    for i = nLE(k):nTE(k)
        yBT(i,k) = 0;
        yBB(i,k) = 0;
    end

    for i = nTE(k)+1:nx
        yBT(i,k) = 0;
        yBB(i,k) = 0;
    end
end

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

for i = 1
    u(i,:,:) = uinf*cos(alpha);
    uold(i,:,:) = uinf*cos(alpha);
end

for i = nx
    u(i,:,:) = uinf*cos(alpha);
    uold(i,:,:) = uinf*cos(alpha);
end

for j = 1
    u(:,j,:) = uinf*cos(alpha);
    uold(:,j,:) = uinf*cos(alpha);
end

for j = ny
    u(:,j,:) = uinf*cos(alpha);
    uold(:,j,:) = uinf*cos(alpha);
end

for k = 1
    u(:,:,k) = uinf*cos(alpha);
    uold(:,:,k) = uinf*cos(alpha);
end

for k = nz
    u(:,:,k) = uinf*cos(alpha);
    uold(:,:,k) = uinf*cos(alpha);
end

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
                D(j) = uinf*cos(alpha);
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
                D(j) = uinf*cos(alpha);
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
                D(j) = uinf*cos(alpha);
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
                D(j) = uinf*cos(alpha);
            end

            u(i,:,k) = tridiagscalar(A,B,C,D);
            u(i,:,k) = uold(i,:,k)+omega*(u(i,:,k)-uold(i,:,k));
            uold(i,:,k) = u(i,:,k);
        end

        % Points on vertical LE line
        for i = nLE(k)
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = uinf*cos(alpha);
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
                        (yBB(i+1,k)-2*yBB(i,k)+yBB(i-1,k))/(dx^2*dy)+0.30/(2*dx);
            end

            % Just after body
            for j = ny/2+1;
                A(j) = 0;
                B(j) = -1/dy^2-2/dx^2-2/dz^2;
                C(j) = 1/dy^2;
                D(j) = -(uold(i-1,j,k)+uold(i+1,j,k))/dx^2 - ...
                        (uold(i,j,k-1)+uold(i,j,k+1))/dz^2 + ...
                        (yBT(i+1,k)-2*yBT(i,k)+yBT(i-1,k))/(dx^2*dy)+0.30/(2*dx);
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
                D(j) = uinf*cos(alpha);
            end

            u(i,:,k) = tridiagscalar(A,B,C,D);
            u(i,:,k) = uold(i,:,k)+omega*(u(i,:,k)-uold(i,:,k));
            uold(i,:,k) = u(i,:,k);
        end

        % Points where airfoil is
        for i = nLE(k)+1:nTE(k)-1
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = uinf*cos(alpha);
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
                D(j) = uinf*cos(alpha);
            end

            u(i,:,k) = tridiagscalar(A,B,C,D);
            u(i,:,k) = uold(i,:,k)+omega*(u(i,:,k)-uold(i,:,k));
            uold(i,:,k) = u(i,:,k);
        end

        % Points on vertical TE line
        for i = nTE(k)
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = uinf*cos(alpha);
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
                D(j) = uinf*cos(alpha);
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
                D(j) = uinf*cos(alpha);
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
                D(j) = uinf*cos(alpha);
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
                D(j) = uinf*cos(alpha);
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
                D(j) = uinf*cos(alpha);
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

    res = max(max(max(resid)))
    resplot(iter) = res;
    iter = iter+1
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
subplot(2,2,1)
semilogy(resplot)
title('Residual versus Iteration Count')
xlabel('Number of Iterations')
ylabel('Residual')
grid on

subplot(2,2,2)
surf(x,z,-cptop')
title('-C_P Across Top of Wing')
% daspect([1 1 1])

subplot(2,2,3)
plot(x,yBB(:,nz/2),'--b')
title('Airfoil Profile')
hold on
plot(x,yBT(:,nz/2),'-b')
daspect([1 1 1])
axis([xmin xmax -4 4])
legend('Bottom of Profile','Top of Profile')
grid on
hold off
