%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 5 - Part 1 - Flat Plate at Angle of Attack
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aoa_deg = 5;
alpha = aoa_deg * pi/180;
uinf = 1;

% ny must be even
nx = 100;
ny = 100;

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;

dx = (xmax-xmin)/(nx-1);
dy = (ymax-ymin)/(ny-1);

x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);

nLE = round(2*nx/5);
nTE = round(3*nx/5)+1;

maxiter = 1000;
omega = 1.91;
resmax = 10^-6;

% Chordlength of airfoil
chord = (x(nTE)-x(nLE));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nLE-1
    yBT(i) = 0;
    yBB(i) = 0;
end

% Flat Plate
for i = nLE:nTE
    yBT(i) = 0;
    yBB(i) = 0;
end

for i = nTE+1:nx
    yBT(i) = 0;
    yBB(i) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for loop = 1:7
    mach_num = (loop-1)*0.1;

    for i = 1:nx
        for j = 1:ny
            Beta(i,j) = sqrt(1-mach_num^2);
        end
    end

    u = zeros(nx,ny);
    uold = zeros(nx,ny);

    for i = 1
        for j = 1:ny
            u(i,j) = uinf*cos(alpha);
            uold(i,j) = uinf*cos(alpha);
        end
    end

    for i = nx
        for j = 1:ny
            u(i,j) = uinf*cos(alpha);
            uold(i,j) = uinf*cos(alpha);
        end
    end

    for i = 1:nx
        for j = 1
            u(i,j) = uinf*cos(alpha);
            uold(i,j) = uinf*cos(alpha);
        end
    end

    for i = 1:nx
        for j = ny
            u(i,j) = uinf*cos(alpha);
            uold(i,j) = uinf*cos(alpha);
        end
    end

    resu = 1;
    iteru = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while (iteru<maxiter && resu>resmax)
        % Points upstream of airfoil
        for i = 2:nLE-1
            for j = 1
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(alpha);
            end

            for j = 2:ny-1
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2-2*(Beta(i,j)^2)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2;
            end

            for j = ny
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(alpha);
            end

            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);
            u(i,:) = uold(i,:)+omega*(u(i,:)-uold(i,:));
            uold(i,:) = u(i,:);
        end

        % Points on vertical LE line
        for i = nLE
            for j = 1
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(alpha);
            end

            % Region below body
            for j = 2:ny/2-1;
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2-2*(Beta(i,j)^2)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2;
            end

            % Just before body
            for j = ny/2;
                Au(j) = 1/dy^2;
                Bu(j) = -1/dy^2-2*(Beta(i,j)^2)/dx^2;
                Cu(j) = 0;
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2-(yBB(i+1)-2*yBB(i)+yBB(i-1))/(dx^2*dy)+0.30/(2*dx);
            end

            % Just after body
            for j = ny/2+1;
                Au(j) = 0;
                Bu(j) = -1/dy^2-2*(Beta(i,j)^2)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2+(yBT(i+1)-2*yBT(i)+yBT(i-1))/(dx^2*dy)-0.30/(2*dx);
            end

            % Region above body
            for j = ny/2+2:ny-1
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2-2*(Beta(i,j)^2)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2;
            end
            for j = ny
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(alpha);
            end
            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);
            u(i,:) = uold(i,:)+omega*(u(i,:)-uold(i,:));
            uold(i,:) = u(i,:);
        end

        % Points where airfoil is
        for i = nLE+1:nTE
            for j = 1
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(alpha);
            end

            % Region below body
            for j = 2:ny/2-1;
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2-2*(Beta(i,j)^2)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2;
            end

            % Just before body
            for j = ny/2;
                Au(j) = 1/dy^2;
                Bu(j) = -1/dy^2-2*(Beta(i,j)^2)/dx^2;
                Cu(j) = 0;
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2-(yBB(i+1)-2*yBB(i)+yBB(i-1))/(dx^2*dy);
            end

            % Just after body
            for j = ny/2+1;
                Au(j) = 0;
                Bu(j) = -1/dy^2-2*(Beta(i,j)^2)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2+(yBT(i+1)-2*yBT(i)+yBT(i-1))/(dx^2*dy);
            end

            % Region above body
            for j = ny/2+2:ny-1
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2-2*(Beta(i,j)^2)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2;
            end

            for j = ny
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(alpha);
            end

            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);
            u(i,:) = uold(i,:)+omega*(u(i,:)-uold(i,:));
            uold(i,:) = u(i,:);
        end

        % Points downstream of airfoil
        for i = nTE+1:nx-1
            for j = 1
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(alpha);
            end

            for j = 2:ny-1
                Au(j) = 1/dy^2;
                Bu(j) = -2/dy^2-2*(Beta(i,j)^2)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2;
            end

            for j = ny
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = uinf*cos(alpha);
            end
            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);
            u(i,:) = uold(i,:)+omega*(u(i,:)-uold(i,:));
            uold(i,:) = u(i,:);
        end

        % Residual
        for i = 2:nx-1
            for j = 2:ny-1
                residu(i,j) = abs((Beta(i,j)^2)*(u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2 + (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2);
            end
        end
        for i = nLE-1:nTE+1
            for j = ny/2-1:ny/2+2
                residu(i,j) = 0;
            end
        end

        resu = max(max(residu));
        resplotu(iteru) = resu;
        iteru = iteru+1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = nLE:nTE
        xcp(i-nLE+1) = x(i);
        cpbot(i-nLE+1) = -2*u(i,ny/2)/uinf;
        cptop(i-nLE+1) = -2*u(i,ny/2+1)/uinf;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(1)
    subplot(2,2,1)
    contour(x,y,u',100)
    title('u Countours Around Airfoil')
    xlabel('u-direction')
    ylabel('v-direction')
    colorbar

    subplot(2,2,2)
    plot(x,u(:,ny/2),'--b')
    hold on
    plot(x,u(:,ny/2+1),'-b')
    title('u Perturbation Velocity Above and Below Airfoil for mach_num = 0:0.6')
    legend('u Perturbation Velocity on Bottom','u Perturbation Velocity on Top')
    grid on

    subplot(2,2,3)
    plot(xcp,-cptop,'-b')
    title('-C_P Across Airfoil for mach_num = 0:0.6')
    xlabel('Location Along Chord')
    ylabel('-C_P')
    hold on
    plot(xcp,-cpbot,'--b')
    legend('C_P Top','C_P Bottom')
    grid on

    subplot(2,2,4)
    semilogy(resplotu,'-b')
    title('Residual versus Iteration Count')
    xlabel('Number of Iterations')
    ylabel('Residual')
    grid on
    hold off

end
