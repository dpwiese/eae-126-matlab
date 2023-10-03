%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 1 - Problem 5 - Flow Over Rotating Cylinder
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section sets up geometric parameters of the cylinder, flow parameters and grid spacing

% Omega is a relaxation parameter
omega = 1.87;
maxiter = 100;
vinf = 1;
rho = 1;

% Specify min and max radius, and number of points in the radial direction
rmin = 1;
rmax = 3;
nr = 12;
dr = (rmax-rmin)/(nr-1);

% Set number of points and grid spacing in the angular direction
thetamin = 0;
ntheta = 50;
dtheta = 2*pi/ntheta;
thetamax = thetamin + ntheta * dtheta;

% Change Gamma
Gamma = -12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r(1) = rmin;
for j=1:nr-1;
    r(j+1) = r(j) + dr;
end

theta(1) = thetamin;

for j=1:ntheta-1;
    theta(j+1) = theta(j) + dtheta;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rubar = ones(ntheta,nr);
rubarold = ones(ntheta,nr);
rvbar = ones(ntheta,nr);
rvbarold = ones(ntheta,nr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U BAR

iter = 1;
while (iter < maxiter)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1;
        for j = 1
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = 0;
        end

        for j = 2:nr-1;
            A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
            B(j) = -(((r(j+1)+2*r(j)+r(j-1))/(2*dr^2))+(2/(r(j)*dtheta^2)));
            C(j) = ((r(j)+r(j+1))/(2*dr^2));
            D(j) = -(rubarold(i+1,j)+rubarold(ntheta,j))/(r(j)*dtheta^2);
        end

        for j = nr;
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = r(j)*vinf*cos(theta(i));
        end

        rubar(i,:) = tridiagscalar(A,B,C,D);
        rubar(i,:) = rubarold(i,:)+omega*(rubar(i,:)-rubarold(i,:));
        rubarold(i,:) = rubar(i,:);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 2:ntheta-1;
        for j = 1
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = 0;
        end
        for j = 2:nr-1;
            A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
            B(j) = -(((r(j+1)+2*r(j)+r(j-1))/(2*dr^2))+(2/(r(j)*dtheta^2)));
            C(j) = ((r(j)+r(j+1))/(2*dr^2));
            D(j) = -(rubarold(i+1,j)+rubarold(i-1,j))/(r(j)*dtheta^2);
        end
        for j = nr;
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = r(j)*vinf*cos(theta(i));
        end
        rubar(i,:) = tridiagscalar(A,B,C,D);
        rubar(i,:) = rubarold(i,:)+omega*(rubar(i,:)-rubarold(i,:));
        rubarold(i,:) = rubar(i,:);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = ntheta;
        for j = 1
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = 0;
        end
        for j = 2:nr-1;
            A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
            B(j) = -(((r(j+1)+2*r(j)+r(j-1))/(2*dr^2))+(2/(r(j)*dtheta^2)));
            C(j) = ((r(j)+r(j+1))/(2*dr^2));
            D(j) = -(rubarold(1,j)+rubarold(i-1,j))/(r(j)*dtheta^2);
        end
        for j = nr;
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = r(j)*vinf*cos(theta(i));
        end
        rubar(i,:) = tridiagscalar(A,B,C,D);
        rubar(i,:) = rubarold(i,:)+omega*(rubar(i,:)-rubarold(i,:));
        rubarold(i,:) = rubar(i,:);
    end

    % for i = 2:ntheta-1
    %     for j = 2:nr-1
    %         residu(i,j) = ((1/(2*dr^2))*(r(j-1)+r(j)))*rubar(i,j-1)-((r(j+1)+2*r(j)+r(j-1))/(2*dr^2)+2/(r(j)*dtheta^2))*rubar(i,j)+((r(j)+r(j+1))/(2*dr^2))*rubar(i,j+1)+(rubar(i+1,j)+rubar(i-1,j))/(r(j)*dtheta^2);
    %     end
    % end
    % erru(iteru) = max(max(abs(residu)));
    % Ru = erru(iteru);

    iter = iter+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V BAR

iter = 1;
while (iter < maxiter)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1;
        for j = 1
            A(j) = 0;
            B(j) = -(((r(j+1)+2*r(j)+(r(j)-dr))/(2*dr^2))+(2/(r(j)*dtheta^2)));
            C(j) = ((r(j)+2*r(j+1)+(r(j)-dr))/(2*dr^2));
            D(j) = -(rvbarold(i+1,j)+rvbarold(ntheta,j))/(r(j)*dtheta^2);
        end

        for j = 2:nr-1;
            A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
            B(j) = -(((r(j+1)+2*r(j)+r(j-1))/(2*dr^2))+(2/(r(j)*dtheta^2)));
            C(j) = ((r(j)+r(j+1))/(2*dr^2));
            D(j) = -(rvbarold(i+1,j)+rvbarold(ntheta,j))/(r(j)*dtheta^2);
        end

        for j = nr;
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = -r(j)*vinf*sin(theta(i))+(Gamma/(2*pi))*(1/1);
        end

        rvbar(i,:) = tridiagscalar(A,B,C,D);
        rvbar(i,:) = rvbarold(i,:)+omega*(rvbar(i,:)-rvbarold(i,:));
        rvbarold(i,:) = rvbar(i,:);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 2:ntheta-1;
        for j = 1
            A(j) = 0;
            B(j) = -(((r(j+1)+2*r(j)+(r(j)-dr))/(2*dr^2))+(2/(r(j)*dtheta^2)));
            C(j) = ((r(j)+2*r(j+1)+(r(j)-dr))/(2*dr^2));
            D(j) = -(rvbarold(i+1,j)+rvbarold(i-1,j))/(r(j)*dtheta^2);
        end

        for j = 2:nr-1;
            A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
            B(j) = -(((r(j+1)+2*r(j)+r(j-1))/(2*dr^2))+(2/(r(j)*dtheta^2)));
            C(j) = ((r(j)+r(j+1))/(2*dr^2));
            D(j) = -(rvbarold(i+1,j)+rvbarold(i-1,j))/(r(j)*dtheta^2);
        end

        for j = nr;
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = -r(j)*vinf*sin(theta(i))+(Gamma/(2*pi))*(1/1);
        end

        rvbar(i,:) = tridiagscalar(A,B,C,D);
        rvbar(i,:) = rvbarold(i,:)+omega*(rvbar(i,:)-rvbarold(i,:));
        rvbarold(i,:) = rvbar(i,:);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = ntheta;
        for j = 1
            A(j) = 0;
            B(j) = -(((r(j+1)+2*r(j)+(r(j)-dr))/(2*dr^2))+(2/(r(j)*dtheta^2)));
            C(j) = ((r(j)+2*r(j+1)+(r(j)-dr))/(2*dr^2));
            D(j) = -(rvbarold(1,j)+rvbarold(ntheta-1,j))/(r(j)*dtheta^2);
        end

        for j = 2:nr-1;
            A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
            B(j) = -(((r(j+1)+2*r(j)+r(j-1))/(2*dr^2))+(2/(r(j)*dtheta^2)));
            C(j) = ((r(j)+r(j+1))/(2*dr^2));
            D(j) = -(rvbarold(1,j)+rvbarold(ntheta-1,j))/(r(j)*dtheta^2);
        end

        for j = nr;
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = -r(j)*vinf*sin(theta(i))+(Gamma/(2*pi))*(1/1);
        end
        
        rvbar(i,:) = tridiagscalar(A,B,C,D);
        rvbar(i,:) = rvbarold(i,:)+omega*(rvbar(i,:)-rvbarold(i,:));
        rvbarold(i,:) = rvbar(i,:);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for i = 2:ntheta-1
    %     for j = 2:nr-1
    %         residv(i,j) = ((1/dr^2)*(r(j-1)+r(j))*0.5)*rvbar(i,j-1)-((r(j+1)+2*r(j)+r(j-1))/(2*dr^2)+2/(r(j)*dtheta^2))*rvbar(i,j)+((r(j)+r(j+1))/(2*dr^2))*rvbar(i,j+1)+(rvbar(i+1,j)+rvbar(i-1,j))/(r(j)*dtheta^2);
    %     end
    % end
    % errv(iterv) = max(max(abs(residv)));
    % Rv = errv(iterv);

    iter = iter+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through angular values
for i = 1:ntheta;
    % Loop through radial values
    for j = 1:nr;
        x(i,j) = r(j)*cos(theta(i));
        y(i,j) = r(j)*sin(theta(i));
    end
end

for i = 1:ntheta
    for j = 1:nr
        ubar(i,j) = rubar(i,j)/r(j);
        vbar(i,j) = rvbar(i,j)/r(j);
    end
end

for i = 1:ntheta
    for j = 1:nr
        ux(i,j) = ubar(i,j)*cos(theta(i));
        uy(i,j) = ubar(i,j)*sin(theta(i));
        vx(i,j) = -vbar(i,j)*sin(theta(i));
        vy(i,j) = vbar(i,j)*cos(theta(i));
    end
end

for i = 1:ntheta
    for j = 1:nr
        xvel(i,j) = ux(i,j)+vx(i,j);
        yvel(i,j) = uy(i,j)+vy(i,j);
        vtot(i,j) = sqrt(xvel(i,j)^2+yvel(i,j)^2);
        ptot(i,j) = 1-((vtot(i,j)^2)/(vinf^2));
    end
end

figure(1)
subplot(2,2,1)
plot(x(:,1),y(:,1),'-k','linewidth',2)
hold on
fill(x(:,1),y(:,1),'r')
quiver(x,y,u,v,1)
title('Velocity Vectors Around Cylinder')
axis([-3,3,-3,3])
xlabel('x-axis: u flow direction')
ylabel('y-axis: v flow direction')
daspect([1 1 1])
hold off

subplot(2,2,2)
contour(x,y,vtot,200)
hold on
plot(x(:,1),y(:,1),'-k','linewidth',2)
fill(x(:,1),y(:,1),'r')
title('Velocity Contours Around Cylinder')
axis([-3,3,-3,3])
xlabel('x-axis: u flow direction')
ylabel('y-axis: v flow direction')
daspect([1 1 1])
hold off

subplot(2,2,3)
contour(x,y,cp,200)
hold on
plot(x(:,1),y(:,1),'-k','linewidth',2)
fill(x(:,1),y(:,1),'r')
title('C_p Surface Contours Around Cylinder')
axis([-3,3,-3,3])
xlabel('x-axis: u flow direction')
ylabel('y-axis: v flow direction')
daspect([1 1 1])
hold off

subplot(2,2,4)
semilogy(erru)
hold on
semilogy(errv,'--k')
grid on
legend('u residual','v residual')
