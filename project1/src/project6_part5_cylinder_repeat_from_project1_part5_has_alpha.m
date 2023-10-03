%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 6 - Problem 5 - Flow Around Cylinder
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section sets up geometric parameters of the cylinder, flow parameters and grid spacing

nr = 30;
rmin = 1;
rmax = 10;
r = linspace(rmin,rmax,nr);
dr = (rmax-rmin)/(nr-1);

ntheta = 40;
thetamin = 0;
thetamax = 2*pi;
theta = linspace(thetamin,thetamax,ntheta);
dtheta = (thetamax-thetamin)/(ntheta-1);

vinf = 1;
rho = 1;
alphadeg = 0;
alpha = alphadeg*pi/180;
a = rmin;

Gamma = -5*pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega = 1.8;
maxiter = 300;

maxerr = 10^-6;
Ru = 1;
Rv = 1;
residu = zeros(ntheta,nr);
residv = zeros(ntheta,nr);

rubar = ones(ntheta,nr);
rubarold = ones(ntheta,nr);
rvbar = ones(ntheta,nr);
rvbarold = ones(ntheta,nr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          U BAR          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iteru = 1;
while (iteru<maxiter && Ru>maxerr)

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
            B(j) = -((r(j+1)+2*r(j)+r(j-1))/(2*dr^2)+2/(r(j)*dtheta^2));
            C(j) = (r(j)+r(j+1))/(2*dr^2);
            D(j) = -(rubarold(i+1,j)+rubarold(ntheta,j))/(r(j)*dtheta^2);
        end

        for j = nr;
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = r(j)*vinf*cos(theta(i)-alpha);
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
            D(j) = r(j)*vinf*cos(theta(i)-alpha);
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
            D(j) = r(j)*vinf*cos(theta(i)-alpha);
        end

        rubar(i,:) = tridiagscalar(A,B,C,D);
        rubar(i,:) = rubarold(i,:)+omega*(rubar(i,:)-rubarold(i,:));
        rubarold(i,:) = rubar(i,:);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 2:ntheta-1
        for j = 2:nr-1
            residu(i,j) = ((1/(2*dr^2))*(r(j-1)+r(j)))*rubar(i,j-1)-((r(j+1)+2*r(j)+r(j-1))/(2*dr^2)+2/(r(j)*dtheta^2))*rubar(i,j)+((r(j)+r(j+1))/(2*dr^2))*rubar(i,j+1)+(rubar(i+1,j)+rubar(i-1,j))/(r(j)*dtheta^2);
        end
    end

    erru(iteru) = max(max(abs(residu)));
    Ru = erru(iteru);
iteru = iteru+1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          V BAR          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iterv = 1;
while (iterv<maxiter && Rv>maxerr)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1;
        for j = 1
            A(j) = 0;
            B(j) = -(((r(j+1)+2*r(j)+(r(j)-dr))/(2*dr^2))+(2/(r(j)*dtheta^2)));
            C(j) = ((r(j+1)+2*r(j)+(r(j)-dr))/(2*dr^2));
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
            D(j) = -r(j)*vinf*sin(theta(i)-alpha)+(Gamma/(2*pi))*(1/1);
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
            C(j) = ((r(j+1)+2*r(j)+(r(j)-dr))/(2*dr^2));
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
            D(j) = -r(j)*vinf*sin(theta(i)-alpha)+(Gamma/(2*pi))*(1/1);
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
            C(j) = ((r(j+1)+2*r(j)+(r(j)-dr))/(2*dr^2));
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
            D(j) = -r(j)*vinf*sin(theta(i)-alpha)+(Gamma/(2*pi))*(1/1);
        end

        rvbar(i,:) = tridiagscalar(A,B,C,D);
        rvbar(i,:) = rvbarold(i,:)+omega*(rvbar(i,:)-rvbarold(i,:));
        rvbarold(i,:) = rvbar(i,:);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 2:ntheta-1
        for j = 2:nr-1
            residv(i,j) = ((1/dr^2)*(r(j-1)+r(j))*0.5)*rvbar(i,j-1)-((r(j+1)+2*r(j)+r(j-1))/(2*dr^2)+2/(r(j)*dtheta^2))*rvbar(i,j)+((r(j)+r(j+1))/(2*dr^2))*rvbar(i,j+1)+(rvbar(i+1,j)+rvbar(i-1,j))/(r(j)*dtheta^2);
        end
    end

    errv(iterv) = max(max(abs(residv)));
    Rv = errv(iterv);
    iterv = iterv+1;
end

for i = 1:ntheta
    for j = 1:nr
        ubar(i,j) = rubar(i,j)/r(j);
        vbar(i,j) = rvbar(i,j)/r(j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through angular values
for i = 1:ntheta ;
    % Loop through radial values
    for j = 1:nr;
        x(i,j) = r(j)*cos(theta(i));
        y(i,j) = r(j)*sin(theta(i));
    end
end

for i = 1:ntheta
    for j = 1:nr
        ubarx(i,j) = ubar(i,j)*cos(theta(i));
        ubary(i,j) = ubar(i,j)*sin(theta(i));
        vbarx(i,j) = -vbar(i,j)*sin(theta(i));
        vbary(i,j) = vbar(i,j)*cos(theta(i));
    end
end

for i = 1:ntheta
    for j = 1:nr
        u(i,j) = ubarx(i,j)+vbarx(i,j);
        v(i,j) = ubary(i,j)+vbary(i,j);
        vtot(i,j) = sqrt(u(i,j)^2+v(i,j)^2);
        cp(i,j) = 1-((vtot(i,j)^2)/(vinf^2));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
