%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 7 - Part 1 - 2-D Supersonic Flow NONLINEAR MACH NUMBER
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aoa_deg = 5;
aoa_rad = aoa_deg*pi/180;
uinf = 1;
rhoinf = 1;

% Set gamma to 1.4 for nonlinear, -1 for linear and cap itermax to a single iteration
gamma = 1.4;
itermax = 100;

% ny must be even
nx = 200;
ny = 100;

xmin = -10;
xmax = 10;
ymin = -20;
ymax = 20;

dx = (xmax-xmin) / (nx-1);
dy = (ymax-ymin) / (ny-1);

x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);

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

sh2 = tau*(x(nTE)+dx/2);
coe2 = sh2/(x(nTE)+dx/2)^2;
chord2 = (x(nTE)-x(nLE))+dx;

% For diamond airfoil
diamond_slope = tau;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nLE-1
    yBT(i) = 0;
    yBB(i) = 0;
end

% Flat Plate
% for i = nLE:nTE
%     yBT(i) = 0;
%     yBB(i) = 0;
% end

% Biconvex
% for i = nLE:nTE
%     yBT(i) = -coe*x(i)^2+sh;
%     yBB(i) = coe*x(i)^2-sh;
% end

% Diamond Airfoil (first half)
for i = nLE:nx/2
    yBT(i) = diamond_slope*x(i)+tau*chord/2;
    yBB(i) = -diamond_slope*x(i)-tau*chord/2;
end

% Diamond Airfoil (second half)
for i = nx/2+1:nTE
    yBT(i) = -diamond_slope*x(i)+tau*chord/2;
    yBB(i) = diamond_slope*x(i)-tau*chord/2;
end

for i = nTE+1:nx
    yBT(i) = 0;
    yBB(i) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nLE-1
    d2yBBdx2(i) = 0;
end

for i = nLE
    d2yBBdx2(i) = (yBB(i+2)-2*yBB(i+1)+yBB(i))/dx^2;
end

for i = nLE+1:nTE-1
    d2yBBdx2(i) = (yBB(i+1)-2*yBB(i)+yBB(i-1))/dx^2;
end

for i = nTE
    d2yBBdx2(i) = (yBB(i)-2*yBB(i-1)+yBB(i-2))/dx^2;
end

for i = nTE+1:nx
    d2yBBdx2(i) = 0;
end

for i = 1:nLE-1
    d2yBTdx2(i) = 0;
end

for i = nLE
    d2yBTdx2(i) = (yBT(i+2)-2*yBT(i+1)+yBT(i))/dx^2;
end

for i = nLE+1:nTE-1
    d2yBTdx2(i) = (yBT(i+1)-2*yBT(i)+yBT(i-1))/dx^2;
end

for i = nTE
    d2yBTdx2(i) = (yBT(i)-2*yBT(i-1)+yBT(i-2))/dx^2;
end

for i = nTE+1:nx
    d2yBTdx2(i) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for Maloop = 1:1
    Minf = 1.6+(Maloop-1)*0.1;

    for i = 1:nx
        for j = 1:ny
            BS(i,j) = (1-Minf^2);
        end
    end

    u = zeros(nx,ny);

    for i = 1:2
        for j = 1:ny
            u(i,j) = uinf*cos(aoa_rad);
        end
    end
    for i = 1:nx
        for j = 1
            u(i,j) = uinf*cos(aoa_rad);
        end
    end
    for i = 1:nx
        for j = ny
            u(i,j) = uinf*cos(aoa_rad);
        end
    end

    % Points upstream of airfoil
    for i = 3:nLE-1
        iter = 1;
        while iter<itermax
            for j = 1
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(aoa_rad);
            end
            for j = 2:ny-1
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2+BS(i,j)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -BS(i,j)*(-2*u(i-1,j)+u(i-2,j))/dx^2;
            end
            for j = ny
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(aoa_rad);
            end
            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);
            for j = 1:ny
                Ma(i,j) = sqrt(Minf^2+(gamma+1)*Minf^2*(u(i,j)-uinf));
                BS(i,j) = (1-Ma(i,j)^2);
            end
            iter = iter+1;
        end
    end

    % Points on vertical LE line
    for i = nLE
        iter = 1;
        while iter<itermax
            for j = 1
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(aoa_rad);
            end

            % Region below body
            for j = 2:ny/2-1;
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2+BS(i,j)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -BS(i,j)*(-2*u(i-1,j)+u(i-2,j))/dx^2;
            end

            % Just before body Bottom
            for j = ny/2;
                Au(j) = 1/dy^2;
                Bu(j) = -1/dy^2+BS(i,j)/dx^2;
                Cu(j) = 0;
                Du(j) = -BS(i,j)*(-2*u(i-1,j)+u(i-2,j))/dx^2 + ...
                            uinf*sin(aoa_rad)/(dx*dy) - ...
                            (yBB(nLE+1)-yBB(nLE))/(dy*dx^2);
            end

            % Just after body Top
            for j = ny/2+1;
                Au(j) = 0;
                Bu(j) = -1/dy^2+BS(i,j)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -BS(i,j)*(-2*u(i-1,j)+u(i-2,j))/dx^2 - ...
                            uinf*sin(aoa_rad)/(dx*dy) + ...
                            (yBT(nLE+1)-yBT(nLE))/(dy*dx^2);
            end

            % Region above body
            for j = ny/2+2:ny-1
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2+BS(i,j)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -BS(i,j)*(-2*u(i-1,j)+u(i-2,j))/dx^2;
            end

            for j = ny
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(aoa_rad);
            end

            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);

            for j = 1:ny
                Ma(i,j) = sqrt(Minf^2+(gamma+1)*Minf^2*(u(i,j)-uinf));
                BS(i,j) = (1-Ma(i,j)^2);
            end

            iter = iter + 1;
        end
    end

    % Points where airfoil is
    for i = nLE+1:nTE-1
        iter = 1;
        while iter<itermax
            for j = 1
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(aoa_rad);
            end

            % Region below body
            for j = 2:ny/2-1;
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2+BS(i,j)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -BS(i,j)*(-2*u(i-1,j)+u(i-2,j))/dx^2;
            end

            % Just before body Bottom
            for j = ny/2;
                Au(j) = 1/dy^2;
                Bu(j) = -1/dy^2+BS(i,j)/dx^2;
                Cu(j) = 0;
                Du(j) = -BS(i,j)*(-2*u(i-1,j)+u(i-2,j))/dx^2-d2yBBdx2(i)/dy;
            end

            % Just after body Top
            for j = ny/2+1;
                Au(j) = 0;
                Bu(j) = -1/dy^2+BS(i,j)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -BS(i,j)*(-2*u(i-1,j)+u(i-2,j))/dx^2+d2yBTdx2(i)/dy;
            end

            % Region above body
            for j = ny/2+2:ny-1
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2+BS(i,j)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -BS(i,j)*(-2*u(i-1,j)+u(i-2,j))/dx^2;
            end

            for j = ny
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(aoa_rad);
            end

            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);

            for j = 1:ny
                Ma(i,j) = sqrt(Minf^2+(gamma+1)*Minf^2*(u(i,j)-uinf));
                BS(i,j) = (1-Ma(i,j)^2);
            end
            iter = iter+1;
        end
    end

    % Points on vertical TE line
    for i = nTE
        iter = 1;
        while iter<itermax
            for j = 1
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(aoa_rad);
            end

            % Region below body
            for j = 2:ny/2-1;
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2+BS(i,j)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -BS(i,j)*(-2*u(i-1,j)+u(i-2,j))/dx^2;
            end

            % Just before body
            for j = ny/2;
                Au(j) = 1/dy^2;
                Bu(j) = -1/dy^2+BS(i,j)/dx^2;
                Cu(j) = 0;
                Du(j) = -BS(i,j)*(-2*u(i-1,j)+u(i-2,j))/dx^2 - ...
                            uinf*sin(aoa_rad)/(dx*dy) + ...
                            (yBB(nTE)-yBB(nTE-1))/(dy*dx^2);
            end

            % Just after body
            for j = ny/2+1;
                Au(j) = 0;
                Bu(j) = -1/dy^2+BS(i,j)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -BS(i,j)*(-2*u(i-1,j)+u(i-2,j))/dx^2 + ...
                            uinf*sin(aoa_rad)/(dx*dy) - ...
                            (yBT(nTE)-yBT(nTE-1))/(dy*dx^2);
            end

            % Region above body
            for j = ny/2+2:ny-1
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2+BS(i,j)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -BS(i,j)*(-2*u(i-1,j)+u(i-2,j))/dx^2;
            end

            for j = ny
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(aoa_rad);
            end

            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);

            for j = 1:ny
                Ma(i,j) = sqrt(Minf^2+(gamma+1)*Minf^2*(u(i,j)-uinf));
                BS(i,j) = (1-Ma(i,j)^2);
            end
            iter = iter+1;
        end
    end

    % Points downstream of airfoil
    for i = nTE+1:nx
        iter = 1;
        while iter<itermax
            for j = 1
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(aoa_rad);
            end

            for j = 2:ny-1
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2+BS(i,j)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -BS(i,j)*(-2*u(i-1,j)+u(i-2,j))/dx^2;
            end

            for j = ny
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(aoa_rad);
            end

            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);

            for j = 1:ny
                Ma(i,j) = sqrt(Minf^2+(gamma+1)*Minf^2*(u(i,j)-uinf));
                BS(i,j) = (1-Ma(i,j)^2);
            end
            iter = iter+1;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = nLE:nTE
        xcp(i-nLE+1) = x(i);
        cpbot(i-nLE+1) = -2*(u(i,ny/2)-uinf)/uinf;
        cptop(i-nLE+1) = -2*(u(i,ny/2+1)-uinf)/uinf;
    end

    uBot = 0;
    uTop = 0;

    for i = nLE:nTE
        % Just before body
        for j = ny/2;
            uBot = uBot+u(i,j);
        end

        % Just after body
        for j = ny/2+1;
            uTop = uTop+u(i,j);
        end
    end

    CM = 0;

    for i = 1:length(xcp)
        cm(i) = (xcp(i)-x(nLE))*(cptop(i)-cpbot(i));
        CM = CM+cm(i);
    end

    Gamma = dx*(uBot-uTop);
    LIFT = rhoinf*uinf*Gamma;
    CL = LIFT/(0.5*rhoinf*uinf^2);
    XCP = CM/CL;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure()
    subplot(2,2,1)
    contour(x,y,u',100)
    title('u Countours Around Airfoil')
    xlabel('u-direction')
    ylabel('v-direction')
    colorbar
    grid on

    subplot(2,2,2)
    plot(x,u(:,ny/2),'--b')
    hold on
    plot(x,u(:,ny/2+1),'-b')
    title('u Velocity Above and Below Airfoil')
    legend('u Velocity on Bottom','u Velocity on Top')
    grid on

    subplot(2,2,3)
    plot(xcp,-cptop,'-b')
    title('-C_P Across Airfoil')
    xlabel('Location Along Chord')
    ylabel('-C_P')
    hold on
    plot(xcp,-cpbot,'--b')
    legend('C_P Top','C_P Bottom')
    grid on

    subplot(2,2,4)
    plot(x,yBB,'--b')
    title('Airfoil Profile')
    hold on
    plot(x,yBT,'-b')
    daspect([1 1 1])
    axis([x(nLE)-1 x(nTE)+1 -sh-1 sh+1])
    legend('Bottom of Profile','Top of Profile')
    grid on
    hold off

end
