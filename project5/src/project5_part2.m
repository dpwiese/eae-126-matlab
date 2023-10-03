%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 5 - THICKNESS - Problem 2: Body of Revolution
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for loop = 1

    Ma = (loop-1)*0.1;

    Ma = 0;
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
    omega = 1.8;
    resmax = 10^-6;

    % tau - thickness ratio of parabolic profiles and NACA
    % sh - vertical shift for parabolic profiles
    % coe - 'A' coefficient for parabolic profiles
    % chord - chordlength of airfoil points on grid
    tau = 0.10;
    sh = tau*(x(nTE)+dx/2);
    coe = sh/x(nLE)^2;
    chord = (x(nTE)-x(nLE));
    aell = chord/2;
    bell = 0.40;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:nLE-1
        yBT(i) = 0;
        yBB(i) = 0;
    end

    % % Biconvex
    % for i = nLE:nTE
    %     yBT(i) = -coe*x(i)^2+sh;
    %     yBB(i) = -yBT(i)
    % end

    % Ellipse
    for i = nLE:nTE
        yBT(i) = bell*sqrt(1-x(i)^2/aell^2);
        yBB(i) = -yBT(i);
    end

    for i = nTE+1:nx
        yBT(i) = 0;
        yBB(i) = 0;
    end

    for i = 1:nLE-1
        S(i) = 0;
     end

    for i = nLE:nTE
        S(i) = pi*yBT(i)^2;
    end

    for i = nTE+1:nx
        S(i) = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:nx
        for j = 1:ny
            Beta(i,j) = sqrt(1-Ma^2);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u = zeros(nx,ny);
    uold = zeros(nx,ny);
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
                A(j) = (y(j)-dy/2)/(y(j)*dy^2);
                B(j) = -(y(j)-dy/2)/(y(j)*dy^2)-(y(j)+dy/2)/(y(j)*dy^2)-2*(Beta(i,j)^2)/dx^2;
                C(j) = (y(j)+dy/2)/(y(j)*dy^2);
                D(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2;
            end
            for j = ny
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end
            u(i,:) = tridiagscalar(A,B,C,D);
            u(i,:) = uold(i,:)+omega*(u(i,:)-uold(i,:));
            uold(i,:) = u(i,:);
        end

        % Points where airfoil is
        for i = nLE:nTE
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end

            % Region below body
            for j = 2:ny/2-1;
                A(j) = (y(j)-dy/2)/(y(j)*dy^2);
                B(j) = -(y(j)-dy/2)/(y(j)*dy^2)-(y(j)+dy/2)/(y(j)*dy^2)-2*(Beta(i,j)^2)/dx^2;
                C(j) = (y(j)+dy/2)/(y(j)*dy^2);
                D(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2;
            end

            % Just before body (bottom)
            for j = ny/2;
                A(j) = (y(j)-dy/2)/(y(j)*dy^2);
                B(j) = -(y(j)-dy/2)/(y(j)*dy^2)-2*(Beta(i,j)^2)/dx^2;
                C(j) = 0;
                D(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2 - ...
                            (S(i+1)-2*S(i)+S(i-1))/(2*pi*dx^2*dy*y(j));
            end

            % Just after body (top)
            for j = ny/2+1;
                A(j) = 0;
                B(j) = -(y(j)+dy/2)/(y(j)*dy^2)-2*(Beta(i,j)^2)/dx^2;
                C(j) = (y(j)+dy/2)/(y(j)*dy^2);
                D(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2 + ...
                            (S(i+1)-2*S(i)+S(i-1))/(2*pi*dx^2*dy*y(j));
            end

            % Region above body
            for j = ny/2+2:ny-1
                A(j) = (y(j)-dy/2)/(y(j)*dy^2);
                B(j) = -(y(j)-dy/2)/(y(j)*dy^2)-(y(j)+dy/2)/(y(j)*dy^2)-2*(Beta(i,j)^2)/dx^2;
                C(j) = (y(j)+dy/2)/(y(j)*dy^2);
                D(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2;
            end

            for j = ny
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end

            u(i,:) = tridiagscalar(A,B,C,D);
            u(i,:) = uold(i,:)+omega*(u(i,:)-uold(i,:));
            uold(i,:) = u(i,:);
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
                A(j) = (y(j)-dy/2)/(y(j)*dy^2);
                B(j) = -(y(j)-dy/2)/(y(j)*dy^2)-(y(j)+dy/2)/(y(j)*dy^2)-2*(Beta(i,j)^2)/dx^2;
                C(j) = (y(j)+dy/2)/(y(j)*dy^2);
                D(j) = -(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2;
            end
            for j = ny
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end
            u(i,:) = tridiagscalar(A,B,C,D);
            u(i,:) = uold(i,:)+omega*(u(i,:)-uold(i,:));
            uold(i,:) = u(i,:);
        end

        % Residual
        for i = 2:nx-1
            for j = 2:ny-1
                resid(i,j) = abs(((y(j)-dy/2)/(y(j)*dy^2))*u(i,j-1)+(-(y(j)-dy/2)/(y(j)*dy^2)-(y(j)+dy/2)/(y(j)*dy^2)-2*(Beta(i,j)^2)/dx^2)*u(i,j)+((y(j)+dy/2)/(y(j)*dy^2))*u(i,j+1)-(-(Beta(i,j)^2)*(uold(i-1,j)+uold(i+1,j))/dx^2));
            end
        end
        for i = nLE-1:nTE+1
            for j = ny/2-1:ny/2+2
                resid(i,j) = 0;
            end
        end

        res = max(max(resid));
        resplot(iter) = res;
        iter = iter+1;
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
    plot(x,yBB,'--b')
    title('Airfoil Profile')
    hold on
    plot(x,yBT,'-b')
    daspect([1 1 1])
    axis([xmin xmax -4 4])
    legend('Bottom of Profile','Top of Profile')
    grid on

    subplot(2,2,2)
    contour(x,y,u',100)
    title('Pressure Contours Around Airfoil')
    xlabel('u-direction')
    ylabel('v-direction')
    colorbar

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
    plot(x,u(:,ny/2),'--b')
    hold on
    plot(x,u(:,ny/2+1),'-b')
    title('Perturbation Velocity Above and Below Airfoil')
    legend('Perturbation Velocity on Bottom','Perturbation Velocity on Top')
    grid on

    figure(2)
    semilogy(resplot)
    hold on
    title('Residual versus Iteration Count')
    xlabel('Number of Iterations')
    ylabel('Residual')
    grid on

end
