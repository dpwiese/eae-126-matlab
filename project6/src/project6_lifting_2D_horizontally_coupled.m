%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 6 - Lifting Problem: 2-D
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uinf = 1;
Ma = 0;
aoa_deg = 5;
alpha = aoa_deg*pi/180;

% ny must be even
nx = 100;
ny = 100;

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;

dx = (xmax - xmin) / (nx - 1);
dy = (ymax - ymin) / (ny - 1);
dt = dx;
eps = dx;

x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);

nLE = round(2*nx/5);
nTE = round(3*nx/5)+1;

maxiter = 30;
omega = 1.91;
resmax = 10^-6;

% tau - thickness ratio of parabolic profiles and NACA
% sh - vertical shift for parabolic profiles points on grid
% coe - 'A' coefficient for parabolic profiles points on grid
% chord - chordlength of airfoil points on grid
tau = 0.1;
sh = tau*(x(nTE)+dx/2);
coe = sh/x(nLE)^2;
chord = (x(nTE)-x(nLE));

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

% Biconvex
for i = nLE:nTE
    yBT(i) = -coe*x(i)^2+sh;
    yBB(i) = coe*x(i)^2-sh;
end

% % NACA 00xx
% for i = nLE:nTE
%     yBT(i) = 10*tau*chord*(0.2969*sqrt((x(i)+chord/2)/chord)-0.1260*((x(i)+chord/2)/chord)-0.3537*((x(i)+chord/2)/chord)^2+0.2843*((x(i)+chord/2)/chord)^3-0.1015*((x(i)+chord/2)/chord)^4);
%     yBB(i) = -yBT(i);
% end

for i = nTE+1:nx
    yBT(i) = 0;
    yBB(i) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nx
    for j = 1:ny
        Beta(i,j) = sqrt(1-Ma^2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = zeros(nx,ny);
uold = zeros(nx,ny);
v = zeros(nx,ny);
vold = zeros(nx,ny);
res = 1;
iter = 1;

while (iter<maxiter && res>resmax)

    % Region below airfoil
    for j = 2:ny/2-1

        % Solve for u
        for i = 1
            Au(i) = 0;
            Bu(i) = 1;
            Cu(i) = 0;
            Du(i) = 0;
        end

        for i = 2:ny-1
            Au(i) = Beta(i,j)/(2*dx)-(eps*Beta(i,j))/dx^2;
            Bu(i) = 1/dt+(2*eps*Beta(i,j))/dx^2+(2*eps)/dy^2;
            Cu(i) = -Beta(i,j)/(2*dx)-(Beta(i,j)*eps)/dx^2;
            Du(i) = uold(i,j)/dt+(vold(i,j+1)-vold(i,j-1))/(2*dy)+eps*(uold(i,j-1)+uold(i,j+1))/dy^2;
        end

        for i = ny
            Au(i) = 0;
            Bu(i) = 1;
            Cu(i) = 0;
            Du(i) = 0;
        end

        u(j,:) = tridiagscalar(Au,Bu,Cu,Du);
        u(j,:) = uold(j,:)+omega*(u(j,:)-uold(j,:));
        uold(j,:) = u(j,:);

        % Solve for v
        for i = 1
            Av(i) = 0;
            Bv(i) = 1;
            Cv(i) = 0;
            Dv(i) = uinf*sin(alpha);
        end

        for i = 2:ny-1
            Av(i) = -eps/dx^2-1/(2*dx);
            Bv(i) = 1/dt+(2*eps)/dx^2+(2*eps)/dy^2;
            Cv(i) = 1/(2*dx)-eps/dx^2;
            Dv(i) = vold(i,j)/dt+(uold(i,j+1)-uold(i,j-1))/(2*dy)+eps*(vold(i,j+1)+vold(i,j-1))/dy^2;
        end

        for i = ny
            Av(i) = 0;
            Bv(i) = 1;
            Cv(i) = 0;
            Dv(i) = uinf*sin(alpha);
        end

        v(j,:) = tridiagscalar(Av,Bv,Cv,Dv);
        v(j,:) = vold(j,:)+omega*(v(j,:)-vold(j,:));
        vold(j,:) = v(j,:);
    end

    % Line just underneath airfoil
    for j = ny/2

        % Solve for u
        for i = 1
            Au(i) = 0;
            Bu(i) = 1;
            Cu(i) = 0;
            Du(i) = 0;
        end

        % Region in front of airfoil
        for i = 2:nLE-1
            Au(i) = Beta(i,j)/(2*dx)-(eps*Beta(i,j))/dx^2;
            Bu(i) = 1/dt+(2*eps*Beta(i,j))/dx^2+(2*eps)/dy^2;
            Cu(i) = -Beta(i,j)/(2*dx)-(Beta(i,j)*eps)/dx^2;
            Du(i) = uold(i,j)/dt+(vold(i,j+1)-vold(i,j-1))/(2*dy)+eps*(uold(i,j-1)+uold(i,j+1))/dy^2;
        end

        % LE point
        for i = nLE;
            Au(i) = 0;
            Bu(i) = 0;
            Cu(i) = 0;
            Du(i) = 0;
        end

        % Region across airfoil
        for i = nLE+1:nTE-1;
            Au(i) = 0;
            Bu(i) = 0;
            Cu(i) = 0;
            Du(i) = 0;
        end

        % TE point
        for i = nTE;
            Au(i) = 0;
            Bu(i) = 0;
            Cu(i) = 0;
            Du(i) = 0;
        end

        % Just after TE point NO KUTTA CONDITION YET - SAME AS EVERYTHING AFTER THIS POINT
        for i = nTE+1;
            Au(i) = 0;
            Bu(i) = 0;
            Cu(i) = 0;
            Du(i) = 0;
        end

        % Region behind airfoil
        for i = nTE+1:nx-1;
            Au(i) = Beta(i,j)/(2*dx)-(eps*Beta(i,j))/dx^2;
            Bu(i) = 1/dt+(2*eps*Beta(i,j))/dx^2+(2*eps)/dy^2;
            Cu(i) = -Beta(i,j)/(2*dx)-(Beta(i,j)*eps)/dx^2;
            Du(i) = uold(i,j)/dt+(vold(i,j+1)-vold(i,j-1))/(2*dy)+eps*(uold(i,j-1)+uold(i,j+1))/dy^2;
        end

        for i = nx
            Au(i) = 0;
            Bu(i) = 1;
            Cu(i) = 0;
            Du(i) = 0;
        end

        u(j,:) = tridiagscalar(Au,Bu,Cu,Du);
        u(j,:) = uold(j,:)+omega*(u(j,:)-uold(j,:));
        uold(j,:) = u(j,:);

        % Solve for v
        for i = 1
            Av(i) = 0;
            Bv(i) = 1;
            Cv(i) = 0;
            Dv(i) = uinf*sin(alpha);
        end

        % Region in front of airfoil
        for i = 2:nLE-1;
            Av(i) = -eps/dx^2-1/(2*dx);
            Bv(i) = 1/dt+(2*eps)/dx^2+(2*eps)/dy^2;
            Cv(i) = 1/(2*dx)-eps/dx^2;
            Dv(i) = vold(i,j)/dt+(uold(i,j+1)-uold(i,j-1))/(2*dy)+eps*(vold(i,j+1)+vold(i,j-1))/dy^2;
        end

        % LE point
        for i = nLE;
            Av(i) = -3/(8*dx)-eps/dx^2;
            Bv(i) = 1/dt+(2*eps)/dx^2+(3*eps)/dy^2;
            Cv(i) = 1/(2*dx)-eps/dx^2;
            Dv(i) = vold(i,j)/dt+(uold(i,j)-uold(i,j-1))/(2*dy)+eps*vold(i,j-1)/dy^2-vold(i-1,j+1)/(8*dx);
        end

        % Region across airfoil
        for i = nLE+1:nTE-1;
            Av(i) = -eps/dx^2-1/(2*dx);
            Bv(i) = 1/dt+(2*eps)/dx^2+(3*eps)/dy^2;
            Cv(i) = 1/(2*dx)-eps/dx^2;
            Dv(i) = vold(i,j)/dt+(uold(i,j)-uold(i,j-1))/(2*dy)+eps*vold(i,j-1)/dy^2;
        end

        % TE point
        for i = nTE;
            Av(i) = -eps/dx^2-1/(2*dx);
            Bv(i) = 1/dt+(2*eps)/dx^2+(3*eps)/dy^2;
            Cv(i) = 3/(8*dx)-eps/dx^2;
            Dv(i) = vold(i,j)/dt+(uold(i,j)-uold(i,j-1))/(2*dy)+eps*vold(i,j-1)/dy^2-vold(i+1,j+1)/(8*dx);
        end

        % Just after TE point NO KUTTA CONDITION YET - SAME AS EVERYTHING AFTER THIS POINT
        for i = nTE+1;
            Av(i) = -eps/dx^2-1/(2*dx);
            Bv(i) = 1/dt+(2*eps)/dx^2+(2*eps)/dy^2;
            Cv(i) = 1/(2*dx)-eps/dx^2;
            Dv(i) = vold(i,j)/dt+(uold(i,j+1)-uold(i,j-1))/(2*dy)+eps*(vold(i,j+1)+vold(i,j-1))/dy^2;
        end

        % Region behind airfoil
        for i = nTE+2:nx-1;
            Av(i) = -eps/dx^2-1/(2*dx);
            Bv(i) = 1/dt+(2*eps)/dx^2+(2*eps)/dy^2;
            Cv(i) = 1/(2*dx)-eps/dx^2;
            Dv(i) = vold(i,j)/dt+(uold(i,j+1)-uold(i,j-1))/(2*dy)+eps*(vold(i,j+1)+vold(i,j-1))/dy^2;
        end

        for i = nx
            Av(i) = 0;
            Bv(i) = 1;
            Cv(i) = 0;
            Dv(i) = uinf*sin(alpha);
        end

        v(j,:) = tridiagscalar(Av,Bv,Cv,Dv);
        v(j,:) = vold(j,:)+omega*(v(j,:)-vold(j,:));
        vold(j,:) = v(j,:);
    end

    % LINE JUST ON TOP OF AIRFOIL
    for j = ny/2+1

        % Solve for u
        for i = 1
            Au(i) = 0;
            Bu(i) = 1;
            Cu(i) = 0;
            Du(i) = 0;
        end

        % Region in front of airfoil
        for i = 2:nLE-1;
            Au(i) = Beta(i,j)/(2*dx)-(eps*Beta(i,j))/dx^2;
            Bu(i) = 1/dt+(2*eps*Beta(i,j))/dx^2+(2*eps)/dy^2;
            Cu(i) = -Beta(i,j)/(2*dx)-(Beta(i,j)*eps)/dx^2;
            Du(i) = uold(i,j)/dt+(vold(i,j+1)-vold(i,j-1))/(2*dy)+eps*(uold(i,j-1)+uold(i,j+1))/dy^2;
        end

        % LE point
        for i = nLE;
            Au(i) = 0;
            Bu(i) = 0;
            Cu(i) = 0;
            Du(i) = 0;
        end

        % Region across airfoil
        for i = nLE+1:nTE-1;
            Au(i) = 0;
            Bu(i) = 0;
            Cu(i) = 0;
            Du(i) = 0;
        end

        % TE point
        for i = nTE;
            Au(i) = 0;
            Bu(i) = 0;
            Cu(i) = 0;
            Du(i) = 0;
        end

        % Just after TE point NO KUTTA CONDITION YET - SAME AS EVERYTHING AFTER THIS POINT
        for i = nTE+1;
            Au(i) = 0;
            Bu(i) = 0;
            Cu(i) = 0;
            Du(i) = 0;
        end

        % Region behind airfoil
        for i = nTE+2:nx-1;
            Au(i) = Beta(i,j)/(2*dx)-(eps*Beta(i,j))/dx^2;
            Bu(i) = 1/dt+(2*eps*Beta(i,j))/dx^2+(2*eps)/dy^2;
            Cu(i) = -Beta(i,j)/(2*dx)-(Beta(i,j)*eps)/dx^2;
            Du(i) = uold(i,j)/dt+(vold(i,j+1)-vold(i,j-1))/(2*dy)+eps*(uold(i,j-1)+uold(i,j+1))/dy^2;
        end

        for i = nx
            Au(i) = 0;
            Bu(i) = 1;
            Cu(i) = 0;
            Du(i) = 0;
        end

        u(j,:) = tridiagscalar(Au,Bu,Cu,Du);
        u(j,:) = uold(j,:)+omega*(u(j,:)-uold(j,:));
        uold(j,:) = u(j,:);

        % Solve for v
        for i = 1
            Av(i) = 0;
            Bv(i) = 1;
            Cv(i) = 0;
            Dv(i) = uinf*sin(alpha);
        end

        % Region in front of airfoil
        for i = 2:nLE-1;
            Av(i) = -eps/dx^2-1/(2*dx);
            Bv(i) = 1/dt+(2*eps)/dx^2+(2*eps)/dy^2;
            Cv(i) = 1/(2*dx)-eps/dx^2;
            Dv(i) = vold(i,j)/dt+(uold(i,j+1)-uold(i,j-1))/(2*dy)+eps*(vold(i,j+1)+vold(i,j-1))/dy^2;
        end

        % LE Point
        for i = nLE;
            Av(i) = -3/(8*dx)-eps/dx^2;
            Bv(i) = 1/dt+(2*eps)/dx^2+(3*eps)/dy^2;
            Cv(i) = 1/(2*dx)-eps/dx^2;
            Dv(i) = vold(i,j)/dt+(uold(i,j+1)-uold(i,j))/(2*dy)+eps*vold(i,j+1)/dy^2-vold(i-1,j-1)/(8*dx);
        end

        % Region across airfoil
        for i = nLE+1:nTE-1;
            Av(i) = -1/(2*dx)-eps/dx^2;
            Bv(i) = 1/dt+(2*eps)/dx^2+(3*eps)/dy^2;
            Cv(i) = 1/(2*dx)-eps/dx^2;
            Dv(i) = vold(i,j)/dt+(uold(i,j+1)-uold(i,j))/(2*dy)+eps*vold(i,j+1)/dy^2;
        end

        % TE Point
        for i = nTE;
            Av(i) = -1/(2*dx)-eps/dx^2;
            Bv(i) = 1/dt+(2*eps)/dx^2+(3*eps)/dy^2;
            Cv(i) = 3/(8*dx)-eps/dx^2;
            Dv(i) = vold(i,j)/dt+(uold(i,j+1)-uold(i,j))/(2*dy)+eps*vold(i,j+1)/dy^2-vold(i+1,j-1)/(8*dx);
        end

        % TE Point%just after TE point NO KUTTA CONDITION YET - SAME AS EVERYTHING AFTER THIS POINT
        for i = nTE+1;
            Av(i) = -eps/dx^2-1/(2*dx);
            Bv(i) = 1/dt+(2*eps)/dx^2+(2*eps)/dy^2;
            Cv(i) = 1/(2*dx)-eps/dx^2;
            Dv(i) = vold(i,j)/dt+(uold(i,j+1)-uold(i,j-1))/(2*dy)+eps*(vold(i,j+1)+vold(i,j-1))/dy^2;
        end

        % Region behind airfoil
        for i = nTE+2:nx-1;
            Av(i) = -eps/dx^2-1/(2*dx);
            Bv(i) = 1/dt+(2*eps)/dx^2+(2*eps)/dy^2;
            Cv(i) = 1/(2*dx)-eps/dx^2;
            Dv(i) = vold(i,j)/dt+(uold(i,j+1)-uold(i,j-1))/(2*dy)+eps*(vold(i,j+1)+vold(i,j-1))/dy^2;
        end

        for i = nx
            Av(i) = 0;
            Bv(i) = 1;
            Cv(i) = 0;
            Dv(i) = uinf*sin(alpha);
        end

        v(j,:) = tridiagscalar(Av,Bv,Cv,Dv);
        v(j,:) = vold(j,:)+omega*(v(j,:)-vold(j,:));
        vold(j,:) = v(j,:);
    end

    % Region above airfoil
    for j = ny/2+2:ny-1

        % Solve for u
        for i = 1
            Au(i) = 0;
            Bu(i) = 1;
            Cu(i) = 0;
            Du(i) = 0;
        end

        for i = 2:ny-1
            Au(i) = Beta(i,j)/(2*dx)-(eps*Beta(i,j))/dx^2;
            Bu(i) = 1/dt+(2*eps*Beta(i,j))/dx^2+(2*eps)/dy^2;
            Cu(i) = -Beta(i,j)/(2*dx)-(Beta(i,j)*eps)/dx^2;
            Du(i) = uold(i,j)/dt+(vold(i,j+1)-vold(i,j-1))/(2*dy)+eps*(uold(i,j-1)+uold(i,j+1))/dy^2;
        end

        for i = ny
            Au(i) = 0;
            Bu(i) = 1;
            Cu(i) = 0;
            Du(i) = 0;
        end

        u(j,:) = tridiagscalar(Au,Bu,Cu,Du);
        u(j,:) = uold(j,:)+omega*(u(j,:)-uold(j,:));
        uold(j,:) = u(j,:);

        % Solve for v
        for i = 1
            Av(i) = 0;
            Bv(i) = 1;
            Cv(i) = 0;
            Dv(i) = uinf*sin(alpha);
        end

        for i = 2:ny-1
            Av(i) = -eps/dx^2-1/(2*dx);
            Bv(i) = 1/dt+(2*eps)/dx^2+(2*eps)/dy^2;
            Cv(i) = 1/(2*dx)-eps/dx^2;
            Dv(i) = vold(i,j)/dt+(uold(i,j+1)-uold(i,j-1))/(2*dy)+eps*(vold(i,j+1)+vold(i,j-1))/dy^2;
        end

        for i = ny
            Av(i) = 0;
            Bv(i) = 1;
            Cv(i) = 0;
            Dv(i) = uinf*sin(alpha);
        end

        v(j,:) = tridiagscalar(Av,Bv,Cv,Dv);
        v(j,:) = vold(j,:)+omega*(v(j,:)-vold(j,:));
        vold(j,:) = v(j,:);
    end

    iter = iter+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(1)
% subplot(2,2,1)
% plot(x,yBB,'--b')
% title('Airfoil Profile')
% hold on
% plot(x,yBT,'-b')
% daspect([1 1 1])
% axis([xmin xmax -4 4])
% legend('Bottom of Profile','Top of Profile')
% grid on

% subplot(2,2,2)
% contour(x,y,u',100)
% title('Pressure Contours Around Airfoil')
% xlabel('u-direction')
% ylabel('v-direction')
% colorbar

% subplot(2,2,3)
% plot(xcp,-cptop,'-b')
% title('-C_P Across Airfoil')
% xlabel('Location Along Chord')
% ylabel('-C_P')
% hold on
% plot(xcp,-cpbot,'--b')
% legend('C_P Top','C_P Bottom')
% grid on

% subplot(2,2,4)
% plot(x,u(:,ny/2),'--b')
% hold on
% plot(x,u(:,ny/2+1),'-b')
% title('Perturbation Velocity Above and Below Airfoil')
% legend('Perturbation Velocity on Bottom','Perturbation Velocity on Top')
% grid on

% figure(2)
% semilogy(resplot)
% hold on
% title('Residual versus Iteration Count')
% xlabel('Number of Iterations')
% ylabel('Residual')
% grid on
