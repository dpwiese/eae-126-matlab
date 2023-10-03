%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 5 - Part 3 - Cross-Flow Problem
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aoa_deg = 5;
alpha = aoa_deg*pi/180;
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

for loop = 1:1

    mach_num = (loop-1)*0.1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:nx
        for j = 1:ny
            Beta(i,j) = sqrt(1-mach_num^2);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u = zeros(nx,ny);
    uold = zeros(nx,ny);
    v = zeros(nx,ny);
    vold = zeros(nx,ny);

    for i = 1
        for j = 1:ny
            v(i,j) = uinf*sin(alpha);
            vold(i,j) = uinf*sin(alpha);
        end
    end

    for i = nx
        for j = 1:ny
            v(i,j) = uinf*sin(alpha);
            vold(i,j) = uinf*sin(alpha);
        end
    end

    for i = 1:nx
        for j = 1
            v(i,j) = uinf*sin(alpha);
            vold(i,j) = uinf*sin(alpha);
        end
    end

    for i = 1:nx
        for j = ny
            v(i,j) = uinf*sin(alpha);
            vold(i,j) = uinf*sin(alpha);
        end
    end

    resu = 1;
    iteru = 1;
    resv = 1;
    iterv = 1;

    %%%%%%%%%%%%%%%%%%%%%%% SOLVE FOR U %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while (iteru<maxiter && resu>resmax)
        % Points upstream of airfoil
        for i = 2:nLE-1
            for j = 1
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = 0;
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
                Du(j) = 0;
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
                Du(j) = 0;
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
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2 - ...
                        (yBB(i+1)-2*yBB(i)+yBB(i-1))/(dx^2*dy) + 0.15/(2*dx);
            end

            % Just after body
            for j = ny/2+1;
                Au(j) = 0;
                Bu(j) = -1/dy^2-2*(Beta(i,j)^2)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2 + ...
                        (yBT(i+1)-2*yBT(i)+yBT(i-1))/(dx^2*dy) - 0.15/(2*dx);
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
                Du(j) = 0;
            end

            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);
            u(i,:) = uold(i,:)+omega*(u(i,:)-uold(i,:));
            uold(i,:) = u(i,:);
        end

        % Points where airfoil is
        for i = nLE+1:nTE-1
            for j = 1
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = 0;
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
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2 - ...
                        (yBB(i+1)-2*yBB(i)+yBB(i-1))/(dx^2*dy);
            end

            % Just after body
            for j = ny/2+1;
                Au(j) = 0;
                Bu(j) = -1/dy^2-2*(Beta(i,j)^2)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2 + ...
                        (yBT(i+1)-2*yBT(i)+yBT(i-1))/(dx^2*dy);
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
                Du(j) = 0;
            end

            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);
            u(i,:) = uold(i,:)+omega*(u(i,:)-uold(i,:));
            uold(i,:) = u(i,:);
        end

        % Points on TE vertical line
        for i = nTE
            for j = 1
                Au(j) = 0;
                Bu(j) = 1;
                Cu(j) = 0;
                Du(j) = 0;
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
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2 - ...
                        (yBB(i+1)-2*yBB(i)+yBB(i-1))/(dx^2*dy)-0.15/(2*dx);
            end

            % Just after body
            for j = ny/2+1;
                Au(j) = 0;
                Bu(j) = -1/dy^2-2*(Beta(i,j)^2)/dx^2;
                Cu(j) = 1/dy^2;
                Du(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2 + ...
                        (yBT(i+1)-2*yBT(i)+yBT(i-1))/(dx^2*dy)+0.15/(2*dx);
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
                Du(j) = 0;
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
                Du(j) = 0;
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
                Du(j) = 0;
            end

            u(i,:) = tridiagscalar(Au,Bu,Cu,Du);
            u(i,:) = uold(i,:)+omega*(u(i,:)-uold(i,:));
            uold(i,:) = u(i,:);
        end

        % Residual
        for i = 2:nx-1
            for j = 2:ny-1
                residu(i,j) = abs( ...
                    (Beta(i,j)^2)*(u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2 + ...
                    (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2);
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

    %%%%%%%%%%%%%%%%%%%%%%% SOLVE FOR V %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while (iterv<maxiter && resv>resmax)
        % Points upstream of airfoil
        for i = 2:nLE-2
            for j = 1
                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = uinf*sin(alpha);
            end

            for j = 2:ny-1
                Av(j) = 1/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 1/dy^2;
                Dv(j) = -(vold(i-1,j)+vold(i+1,j))/dx^2;
            end

            for j = ny
                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = uinf*sin(alpha);
            end

            v(i,:) = tridiagscalar(Av,Bv,Cv,Dv);
            v(i,:) = vold(i,:)+omega*(v(i,:)-vold(i,:));
            vold(i,:) = v(i,:);
        end

        % Line just before LE
        for i = nLE-1
            for j = 1
                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = uinf*sin(alpha);
            end

            for j = 2:ny/2-1
                Av(j) = 1/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 1/dy^2;
                Dv(j) = -(vold(i-1,j)+vold(i+1,j))/dx^2;
            end

            for j = ny/2
                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = 0.15;
            end

            for j = ny/2+1:ny-1
                Av(j) = 1/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 1/dy^2;
                Dv(j) = -(vold(i-1,j)+vold(i+1,j))/dx^2;
            end

            for j = ny
                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = uinf*sin(alpha);
            end

            v(i,:) = tridiagscalar(Av,Bv,Cv,Dv);
            v(i,:) = vold(i,:)+omega*(v(i,:)-vold(i,:));
            vold(i,:) = v(i,:);
        end

        % Points Across Airfoil
        for i = nLE:nTE
            for j = 1
                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = uinf*sin(alpha);
            end

            % Region below body
            for j = 2:ny/2-1;
                Av(j) = 1/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 1/dy^2;
                Dv(j) = -(vold(i-1,j)+vold(i+1,j))/dx^2;
            end

            % The body
            for j = ny/2;
                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = 0;
            end

            % Region above the body
            for j = ny/2+1:ny-1;
                Av(j) = 1/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 1/dy^2;
                Dv(j) = -(vold(i-1,j)+vold(i+1,j))/dx^2;
            end

            for j = ny
                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = uinf*sin(alpha);
            end

            v(i,:) = tridiagscalar(Av,Bv,Cv,Dv);
            v(i,:) = vold(i,:)+omega*(v(i,:)-vold(i,:));
            vold(i,:) = v(i,:);
        end

        % Line just after TE
        for i = nTE+1
            for j = 1
                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = uinf*sin(alpha);
            end

            for j = 2:ny/2-1
                Av(j) = 1/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 1/dy^2;
                Dv(j) = -(vold(i-1,j)+vold(i+1,j))/dx^2;
            end

            for j = ny/2
                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = 0.15;
            end

            for j = ny/2+1:ny-1
                Av(j) = 1/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 1/dy^2;
                Dv(j) = -(vold(i-1,j)+vold(i+1,j))/dx^2;
            end

            for j = ny
                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = uinf*sin(alpha);
            end

            v(i,:) = tridiagscalar(Av,Bv,Cv,Dv);
            v(i,:) = vold(i,:)+omega*(v(i,:)-vold(i,:));
            vold(i,:) = v(i,:);
        end

        % Points downstream of airfoil
        for i = nTE+2:nx-1
            for j = 1
                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = uinf*sin(alpha);
            end

            for j = 2:ny-1
                Av(j) = 1/dy^2;
                Bv(j) = -2/dy^2-2/dx^2;
                Cv(j) = 1/dy^2;
                Dv(j) = -(vold(i-1,j)+vold(i+1,j))/dx^2;
            end

            for j = ny
                Av(j) = 0;
                Bv(j) = 1;
                Cv(j) = 0;
                Dv(j) = uinf*sin(alpha);
            end

            v(i,:) = tridiagscalar(Av,Bv,Cv,Dv);
            v(i,:) = vold(i,:)+omega*(v(i,:)-vold(i,:));
            vold(i,:) = v(i,:);
        end

        % Residual
        for i = 2:nx-1
            for j = 2:ny-1
                residv(i,j) = abs( ...
                    (v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2 + (v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2);
            end
        end
        for i = nLE-1:nTE+1
            for j = ny/2-1:ny/2+2
                residv(i,j) = 0;
            end
        end

        resv = max(max(residv));
        resplotv(iterv) = resv;
        iterv = iterv+1;
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
    title('u Perturbation Velocity Above and Below Airfoil')
    legend('u Perturbation Velocity on Bottom','u Perturbation Velocity on Top')
    grid on

    subplot(2,2,3)
    contour(x,y,v',100)
    title('v Contours Around Airfoil')
    xlabel('u-direction')
    ylabel('v-direction')
    colorbar

    subplot(2,2,4)
    plot(x,v(:,ny/2),'--b')
    hold on
    plot(x,v(:,ny/2+1),'-b')
    title('v Perturbation Velocity Above and Below Airfoil')
    legend('v Perturbation Velocity on Bottom','v Perturbation Velocity on Top')
    grid on

    figure(2)
    plot(xcp,-cptop,'-b')
    title('-C_P Across Airfoil')
    xlabel('Location Along Chord')
    ylabel('-C_P')
    hold on
    plot(xcp,-cpbot,'--b')
    legend('C_P Top','C_P Bottom')
    grid on

    figure(3)
    semilogy(resplotv,'--b')
    hold on
    semilogy(resplotu,'-b')
    title('Residual versus Iteration Count')
    xlabel('Number of Iterations')
    ylabel('Residual')
    legend('v-Residual','u-Residual')
    grid on
    hold off

end
