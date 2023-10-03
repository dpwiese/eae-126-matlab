%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 1 - Problem 6 - Flow Over Cylinder with Fins
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section sets up geometric parameters of the cylinder, flow parameters and grid spacing

omega = 1.00;
maxiter = 500;
Gamma = 0;
uinf = 1;

rmin = 1;
rmax = 4*rmin;
nr = 22;
nf = round(2*nr/3);
dr = (rmax-rmin)/(nr-1);
thetamin = 0;
ntheta = 60;
dtheta = 2*pi/ntheta;
thetamax = thetamin+ntheta*dtheta;
thetaftop = thetamax/4;
thetafbot = 3*thetamax/4;
lf = 2*rmin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r(1) = rmin;

for j = 1:nr-1;
    r(j+1) = r(j)+dr;
end

theta(1) = thetamin;

for j = 1:ntheta-1;
    theta(j+1) = theta(j)+dtheta;
end

rubar = ones(ntheta,nr);
rubarold = ones(ntheta,nr);
rvbar = ones(ntheta,nr);
rvbarold = ones(ntheta,nr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          U BAR          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iter = 1;
while (iter<maxiter)
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
            D(j) = r(j)*uinf*cos(theta(i));
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
            D(j) = r(j)*uinf*cos(theta(i));
        end

        % Top fin
        if i == round(ntheta/4)
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end

            for j = 2:nf;
                A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
                B(j) = -(((r(j+1)+2*r(j)+r(j-1))/(2*dr^2))+(2/(r(j)*dtheta^2)));
                C(j) = ((r(j)+r(j+1))/(2*dr^2));
                D(j) = -(rubarold(i-1,j)+rubarold(i+1,j))/(r(j)*dtheta^2);
            end

            for j = nf+1:nr;
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = r(j)*uinf*cos(theta(i));
            end
        end

        if i == round(ntheta/4)+1
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end

            for j = 2:nf;
                A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
                B(j) = -(((r(j+1)+2*r(j)+r(j-1))/(2*dr^2))+(2/(r(j)*dtheta^2)));
                C(j) = ((r(j)+r(j+1))/(2*dr^2));
                D(j) = -(2*rubarold(i-1,j))/(r(j)*dtheta^2); % -(rubarold(i+1,j)-rubarold(i-1,j))/(r(j)*dtheta^2);%
            end

            for j = nf+1:nr;
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = r(j)*uinf*cos(theta(i));
            end
        end

        if i == round(ntheta/4)+2
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end

            for j = 2:nf;
                A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
                B(j) = -(((r(j+1)+2*r(j)+r(j-1))/(2*dr^2))+(2/(r(j)*dtheta^2)));
                C(j) = ((r(j)+r(j+1))/(2*dr^2));
                D(j) = -(rubarold(i+1,j)+rubarold(i-1,j))/(r(j)*dtheta^2);
            end

            for j = nf+1:nr;
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = r(j)*uinf*cos(theta(i));
            end
        end

        % Bottom fin
        if i == round(3*ntheta/4)
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end

            for j = 2:nf;
                A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
                B(j) = -(((r(j+1)+2*r(j)+r(j-1))/(2*dr^2))+(2/(r(j)*dtheta^2)));
                C(j) = ((r(j)+r(j+1))/(2*dr^2));
                D(j) = -(rubarold(i-1,j)+rubarold(i+1,j))/(r(j)*dtheta^2);
            end

            for j = nf+1:nr;
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = r(j)*uinf*cos(theta(i));
            end
        end

        if i == round(3*ntheta/4)+1
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end

            for j = 2:nf;
                A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
                B(j) = -(((r(j+1)+2*r(j)+r(j-1))/(2*dr^2))+(2/(r(j)*dtheta^2)));
                C(j) = ((r(j)+r(j+1))/(2*dr^2));
                D(j) = -(2*rubarold(i-1,j))/(r(j)*dtheta^2); % -(rubarold(i+1,j)-rubarold(i-1,j))/(r(j)*dtheta^2);%
            end

            for j = nf+1:nr;
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = r(j)*uinf*cos(theta(i));
            end
        end

        if i == round(3*ntheta/4)+2
            for j = 1
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end

            for j = 2:nf;
                A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
                B(j) = -(((r(j+1)+2*r(j)+r(j-1))/(2*dr^2))+(2/(r(j)*dtheta^2)));
                C(j) = ((r(j)+r(j+1))/(2*dr^2));
                D(j) = -(rubarold(i+1,j)+rubarold(i-1,j))/(r(j)*dtheta^2);
            end

            for j = nf+1:nr;
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = r(j)*uinf*cos(theta(i));
            end
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
            D(j) = r(j)*uinf*cos(theta(i));
        end

        rubar(i,:) = tridiagscalar(A,B,C,D);
        rubar(i,:) = rubarold(i,:)+omega*(rubar(i,:)-rubarold(i,:));
        rubarold(i,:) = rubar(i,:);
    end

    iter = iter+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          V BAR          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iter = 1;
while (iter<maxiter)
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
            D(j) = -r(j)*uinf*sin(theta(i))+(Gamma/(2*pi))*(1/1);
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
            D(j) = -r(j)*uinf*sin(theta(i))+(Gamma/(2*pi))*(1/1);
        end

        % Top fin
        if i == round(ntheta/4)+1
            for j = 1:nf;
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end

            for j = nf+1:nr;
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = -r(j)*uinf*sin(theta(i))+(Gamma/(2*pi))*(1/1);
            end
        end

        % Bottom fin
        if i == round(3*ntheta/4)+1
            for j = 1:nf;
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = 0;
            end

            for j = nf+1:nr;
                A(j) = 0;
                B(j) = 1;
                C(j) = 0;
                D(j) = -r(j)*uinf*sin(theta(i))+(Gamma/(2*pi))*(1/1);
            end
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
            D(j) = -r(j)*uinf*sin(theta(i))+(Gamma/(2*pi))*(1/1);
        end
        
        rvbar(i,:) = tridiagscalar(A,B,C,D);
        rvbar(i,:) = rvbarold(i,:)+omega*(rvbar(i,:)-rvbarold(i,:));
        rvbarold(i,:) = rvbar(i,:);
    end

    iter = iter+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:ntheta ;
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
        ptot(i,j) = 1-((vtot(i,j)^2)/(uinf^2));
    end
end

% quiver(x,y,xvel,yvel,1)
% axis([-rmax,rmax,-rmax,rmax])
% axis square

figure(1)
subplot(2,2,1)
quiver(x,y,xvel,yvel,1)
title('Velocity Vectors')
colorbar
axis square
axis([-4,4,-4,4])

subplot(2,2,2)
contour(x,y,ptot,200)
title('Pressure Contour')
colorbar
axis square
axis([-4,4,-4,4])

subplot(2,2,3)
contour(x,y,vtot,200)
title('Velocity Contour')
colorbar
axis square
axis([-4,4,-4,4])

subplot(2,2,4)
quiver(x,y,xvel,yvel,1)
hold on
contour(x,y,ptot,200)
title('Velocity Vectors and Pressure Contours')
colorbar
axis square
axis([-4,4,-4,4])
hold off
