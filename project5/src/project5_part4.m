%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 5 - Part 1 - Problem 4: Swept Flat Wing at aoa_deg
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uinf = 1;
aoa_deg = 0;
Ma = 0.0;

% nx, ny, and nz all must be even numbers
nx = 40;
ny = 20;
nz = 40;

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
zmin = -10;
zmax = 10;

dx = (xmax-xmin) / (nx-1);
dy = (ymax-ymin) / (ny-1);
dz = (zmax-zmin) / (nz-1);

x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);
z = linspace(zmin, zmax, nz);

span2 = 10;
chord2 = 5;

for k = 1:nz/2
    [rspan] = find(z(k+1)>-(span2/2) & z(k)<-span2/2);
    rlens(k) = length(rspan);
end

nLEFT = find(rlens);
nRGHT = nz-nLEFT+1;

for i = 1:nx/2
    [rchord] = find(x(i+1)>-(chord2/2) & x(i)<-chord2/2);
    rlenc(i) = length(rchord);
end

nLE(nz/2) = find(rlenc);
nRGHT = nz-nLEFT+1;

% nLEFT = round(2*nz/5);
% nRGHT = round(3*nz/5)+1;
nSWEP = nz/2-nLEFT+1;
nCHRD = round(3*nx/5)+1-round(2*nx/5);

maxiter = 300;
omega = 1.95;
resmax = 10^-6;

% Chordlength of airfoil POINTS ON GRID
chord = (x(round(3*nx/5)+1)-x(round(2*nx/5)));
span = z(nRGHT)-z(nLEFT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:nLEFT-1
    nLE(k) = 0;
    nTE(k) = 0;
end

for k = nz/2:nz/2+1
    nLE(k) = round(2*nx/5);
    nTE(k) = nLE(k)+nCHRD;
end

for k = nz/2+2:nRGHT
    nLE(k) = nLE(nz/2+1)+k-nz/2-1;
    nTE(k) = nLE(k)+nCHRD;
end

for k = nLEFT:nz/2-1
    nLE(k) = nLE(nz-k+1);
    nTE(k) = nLE(k)+nCHRD;
end

for k = nRGHT+1:nz
    nLE(k) = 0;
    nTE(k) = 0;
end

figure()
plot(nLE)
hold on
plot(nTE)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:nz
    for i = 1:nx
        yBT(i,k) = 0;
        yBB(i,k) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 1:nLE-1
%     alpha(i) = 0;
% end
% for i = nLE:nTE
%     alpha(i) = aoa_deg*pi/180;
% end
% for i = nTE+1:nx
%     alpha(i) = 0;
% end
% Gamma = 0;

% for k = 1:nz
%     for i = 1:nx
%         for j = 1:ny
%             Beta(i,j,k) = sqrt(1-Ma);
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = zeros(nx,ny,nz);
uold = zeros(nx,ny,nz);
res = 1;
iter = 1;

while (iter<maxiter && res>resmax)

    % Planes before left tip of wing
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
    %     for j = ny/2; %just before body
    %         uBot = uBot+u(i,j);
    %     end
    %     for j = ny/2+1; %just after body
    %         uTop = uTop+u(i,j);
    %     end
    % end
    Gamma = dx*(uBot-uTop);

    res  =  max(max(max(resid)))
    resplot(iter)  =  res;
    iter  =  iter + 1
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
        cpbot(i, k) = -2 * u(i, ny/2, k) / uinf;
        cptop(i, k) = -2 * u(i, ny/2+1, k) / uinf;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(x,yBB(:,nz/2),'--b')
title('Airfoil Profile')
hold on
plot(x,yBT(:,nz/2),'-b')
daspect([1 1 1])
axis([xmin xmax -4 4])
legend('Bottom of Profile','Top of Profile')
grid on
hold off

figure()
semilogy(resplot)
title('Residual versus Iteration Count')
xlabel('Number of Iterations')
ylabel('Residual')
grid on

figure()
surf(x,z,-cptop')
title('-C_P Across Top of Wing')

% for i = 2:nz-1
%     figure(i+2)
%     contour(x,y,u(:,:,i)',100)
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

% figure()
% plot(xcp,-cptop,'-b')
% title('-C_P Across Airfoil')
% xlabel('Location Along Chord')
% ylabel('-C_P')
% hold on
% plot(xcp,-cpbot,'--b')
% legend('C_P Top','C_P Bottom')
% grid on
% hold off

% figure()
% plot(x,alpha)
% title('Angle of Attack on Airfoil')
% grid on

% figure()
% plot(x,u(:,ny/2),'--b')
% hold on
% plot(x,u(:,ny/2+1),'-b')
% title('Perturbation Velocity Above and Below Airfoil')
% legend('Perturbation Velocity on Bottom','Perturbation Velocity on Top')
% grid on
% hold off

% figure()
% plot(x,d2yBBdx2,'--b')
% hold on
% plot(x,d2yBTdx2,'-b')
% grid on
% title('Second Derivative of Airfoil')
% legend('Second Derivative on Bottom','Second Derivative on Top')
% hold off
