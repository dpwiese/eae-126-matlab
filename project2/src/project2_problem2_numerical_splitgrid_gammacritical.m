%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 2 - Problem 2 - Cp, Velocity Vectors, Etc.
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section sets up geometric parameters of the cylinder, flow parameters and grid spacing

nr = 20;
rmin = 1;
rmax = 4;
r = linspace(rmin,rmax,nr);
dr = (rmax-rmin)/(nr-1);

ntheta = 4*35;
dtheta = 2*pi/ntheta;
thetamin = 0.5*dtheta;
thetamax = 2*pi-dtheta;
theta = linspace(thetamin,thetamax,ntheta);

vinf = 1;
rho = 1;
alphadeg = 0;
alpha = alphadeg*pi/180;
a = rmin;
pinf = 0;

epsilon = -0.1;
mu = 0.1;
b = sqrt(a^2-mu^2)+epsilon;

% epsilon = 0;
% mu = 0;
% b = 0.9; %use this b with mu and epsilon zero to make ellipse with b<a

D = vinf*rmin^2;
Gamma = -4*pi;
thetasep = -alpha-asin(mu/a);
% Gamma = 4*pi*vinf*rmin*sin(thetasep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega = 1.87;
maxiter = 500;
rubar = ones(ntheta,nr);
rubarold = ones(ntheta,nr);
rvbar = ones(ntheta,nr);
rvbarold = ones(ntheta,nr);
maxerr = 10^-6;
Ru = 1;
Rv = 1;
% rubartemp = ones(ntheta,nr);
% rvbartemp = ones(ntheta,nr);
residu = zeros(ntheta,nr);
residv = zeros(ntheta,nr);

for i = 1:ntheta
    rubar(i,nr) = r(nr)*vinf*cos(theta(i)-alpha);
    rvbar(i,nr) = -r(nr)*vinf*sin(theta(i)-alpha)+Gamma/(2*pi);
    rubarold(i,nr) = r(nr)*vinf*cos(theta(i)-alpha);
    rvbarold(i,nr) = -r(nr)*vinf*sin(theta(i)-alpha)+Gamma/(2*pi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          U BAR          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iteru = 1;
while (iteru<maxiter && Ru>maxerr);

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

        dd(1) = D(1);
        bb(1) = B(1);

        for j = 2:nr
            bb(j) = B(j)-A(j)*C(j-1)/bb(j-1);
            dd(j) = D(j)-A(j)*dd(j-1)/bb(j-1);
        end

        gs(nr) = dd(nr)/bb(nr);

        for j = nr-1:-1:1
            gs(j) = (dd(j)-rubar(j+1)*C(j))/bb(j);
        end

        for j = 1:nr-1
            rubar(i,j) = gs(j);
        end

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
            B(j) = -((r(j+1)+2*r(j)+r(j-1))/(2*dr^2)+2/(r(j)*dtheta^2));
            C(j) = (r(j)+r(j+1))/(2*dr^2);
            D(j) = -(rubarold(i+1,j)+rubarold(i-1,j))/(r(j)*dtheta^2);
        end

        for j = nr;
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = r(j)*vinf*cos(theta(i)-alpha);
        end

        dd(1) = D(1);
        bb(1) = B(1);

        for j = 2:nr
            bb(j) = B(j)-A(j)*C(j-1)/bb(j-1);
            dd(j) = D(j)-A(j)*dd(j-1)/bb(j-1);
        end

        gs(nr) = dd(nr)/bb(nr);

        for j = nr-1:-1:1
            gs(j) = (dd(j)-rubar(j+1)*C(j))/bb(j);
        end

        for j = 1:nr-1
            rubar(i,j) = gs(j);
        end

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
            B(j) = -((r(j+1)+2*r(j)+r(j-1))/(2*dr^2)+2/(r(j)*dtheta^2));
            C(j) = (r(j)+r(j+1))/(2*dr^2);
            D(j) = -(rubarold(1,j)+rubarold(i-1,j))/(r(j)*dtheta^2);
        end

        for j = nr;
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = r(j)*vinf*cos(theta(i)-alpha);
        end

        dd(1) = D(1);
        bb(1) = B(1);

        for j = 2:nr
            bb(j) = B(j)-A(j)*C(j-1)/bb(j-1);
            dd(j) = D(j)-A(j)*dd(j-1)/bb(j-1);
        end

        gs(nr) = dd(nr)/bb(nr);

        for j = nr-1:-1:1
            gs(j) = (dd(j)-rubar(j+1)*C(j))/bb(j);
        end

        for j = 1:nr-1
            rubar(i,j) = gs(j);
        end

        rubarold(i,:) = rubar(i,:);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 2:ntheta-1
        for j = 2:nr-1
            residu(i,j) = [((1/dr^2)*(r(j-1)+r(j))*0.5)]*rubar(i,j+1)-[(r(j+1)+2*r(j)+r(j-1))/(2*dr^2)+2/(r(j)*dtheta^2)]*rubar(i,j)+((r(j)+r(j+1))/(2*dr^2))*rubar(i,j-1)+(rubar(i+1,j)+rubar(i-1,j))/(r(j)*dtheta^2);
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
            B(j) = -((r(j+1)+2*r(j)+(r(j)-dr))/(2*dr^2)+2/(r(j)*dtheta^2));
            C(j) = (r(j)+2*r(j+1)+(r(j)-dr))/(2*dr^2);
            D(j) = -(rvbarold(i+1,j)+rvbarold(ntheta,j))/(r(j)*dtheta^2);
        end

        for j = 2:nr-1;
            A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
            B(j) = -((r(j+1)+2*r(j)+r(j-1))/(2*dr^2)+2/(r(j)*dtheta^2));
            C(j) = (r(j)+r(j+1))/(2*dr^2);
            D(j) = -(rvbarold(i+1,j)+rvbarold(ntheta,j))/(r(j)*dtheta^2);
        end

        for j = nr;
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = -r(j)*vinf*sin(theta(i)-alpha)+Gamma/(2*pi);
        end

        dd(1) = D(1);
        bb(1) = B(1);

        for j = 2:nr
            bb(j) = B(j)-A(j)*C(j-1)/bb(j-1);
            dd(j) = D(j)-A(j)*dd(j-1)/bb(j-1);
        end

        gs(nr) = dd(nr)/bb(nr);

        for j = nr-1:-1:1
            gs(j) = (dd(j)-rvbar(j+1)*C(j))/bb(j);
        end

        for j = 1:nr-1
            rvbar(i,j) = gs(j);
        end

        rvbarold(i,:) = rvbar(i,:);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 2:ntheta-1;
        for j = 1
            A(j) = 0;
            B(j) = -((r(j+1)+2*r(j)+(r(j)-dr))/(2*dr^2)+2/(r(j)*dtheta^2));
            C(j) = (r(j)+2*r(j+1)+(r(j)-dr))/(2*dr^2);
            D(j) = -(rvbarold(i+1,j)+rvbarold(i-1,j))/(r(j)*dtheta^2);
        end

        for j = 2:nr-1;
            A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
            B(j) = -((r(j+1)+2*r(j)+r(j-1))/(2*dr^2)+2/(r(j)*dtheta^2));
            C(j) = (r(j)+r(j+1))/(2*dr^2);
            D(j) = -(rvbarold(i+1,j)+rvbarold(i-1,j))/(r(j)*dtheta^2);
        end

        for j = nr;
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = -r(j)*vinf*sin(theta(i)-alpha)+Gamma/(2*pi);
        end

        dd(1) = D(1);
        bb(1) = B(1);

        for j = 2:nr
            bb(j) = B(j)-A(j)*C(j-1)/bb(j-1);
            dd(j) = D(j)-A(j)*dd(j-1)/bb(j-1);
        end

        gs(nr) = dd(nr)/bb(nr);

        for j = nr-1:-1:1
            gs(j) = (dd(j)-rvbar(j+1)*C(j))/bb(j);
        end

        for j = 1:nr-1
            rvbar(i,j) = gs(j);
        end

        rvbarold(i,:) = rvbar(i,:);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = ntheta;
        for j = 1
            A(j) = 0;
            B(j) = -((r(j+1)+2*r(j)+(r(j)-dr))/(2*dr^2)+2/(r(j)*dtheta^2));
            C(j) = (r(j)+2*r(j+1)+(r(j)-dr))/(2*dr^2);
            D(j) = -(rvbarold(1,j)+rvbarold(i-1,j))/(r(j)*dtheta^2);
        end

        for j = 2:nr-1;
            A(j) = (1/dr^2)*(r(j-1)+r(j))*0.5;
            B(j) = -((r(j+1)+2*r(j)+r(j-1))/(2*dr^2)+2/(r(j)*dtheta^2));
            C(j) = (r(j)+r(j+1))/(2*dr^2);
            D(j) = -(rvbarold(1,j)+rvbarold(i-1,j))/(r(j)*dtheta^2);
        end

        for j = nr;
            A(j) = 0;
            B(j) = 1;
            C(j) = 0;
            D(j) = -r(j)*vinf*sin(theta(i)-alpha)+Gamma/(2*pi);
        end

        dd(1) = D(1);
        bb(1) = B(1);

        for j = 2:nr
            bb(j) = B(j)-A(j)*C(j-1)/bb(j-1);
            dd(j) = D(j)-A(j)*dd(j-1)/bb(j-1);
        end

        gs(nr) = dd(nr)/bb(nr);

        for j = nr-1:-1:1
            gs(j) = (dd(j)-rvbar(j+1)*C(j))/bb(j);
        end

        for j = 1:nr-1
            rvbar(i,j) = gs(j);
        end
        
        rvbarold(i,:) = rvbar(i,:);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 2:ntheta-1
        for j = 2:nr-1
            residv(i,j) = ((1/dr^2)*(r(j-1)+r(j))*0.5)*rvbar(i,j+1)-((r(j+1)+2*r(j)+r(j-1))/(2*dr^2)+2/(r(j)*dtheta^2))*rvbar(i,j)+((r(j)+r(j+1))/(2*dr^2))*rvbar(i,j-1)+(rvbar(i+1,j)+rvbar(i-1,j))/(r(j)*dtheta^2);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        p(i,j) = 0.5*cp(i,j)*rho*vinf^2+pinf;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This sets up transformed coordinate system
for i = 1:ntheta
    for j = 1:nr
        xair(i,j) = (x(i,j)+epsilon);
        yair(i,j) = (y(i,j)+mu);
    end
end

% This sets up transformed coordinate system
for i = 1:ntheta
    for j = 1:nr
        X(i,j) = (xair(i,j))*(1+(b^2)/((xair(i,j))^2+(yair(i,j))^2));
        Y(i,j) = (yair(i,j))*(1-(b^2)/((xair(i,j))^2+(yair(i,j))^2));
        dXdx(i,j) = 1-2*(xair(i,j)^2)*(b^2)/(xair(i,j)^2+yair(i,j)^2)^2+(b^2)/(xair(i,j)^2+yair(i,j)^2);
        dYdy(i,j) = 1+2*(yair(i,j)^2)*(b^2)/(xair(i,j)^2+yair(i,j)^2)^2-(b^2)/(xair(i,j)^2+yair(i,j)^2);
        dXdy(i,j) = -2*(b^2)*xair(i,j)*yair(i,j)/(xair(i,j)^2+yair(i,j)^2)^2;
        dYdx(i,j) = 2*(b^2)*xair(i,j)*yair(i,j)/(xair(i,j)^2+yair(i,j)^2)^2;
    end
end

% Solve matrix using cramers rule
for i = 1:ntheta
    for j = 1:nr
        U(i,j) = (u(i,j)*dYdy(i,j)-dYdx(i,j)*v(i,j))/(dXdx(i,j)*dYdy(i,j)-dYdx(i,j)*dXdy(i,j));
        V(i,j) = (dXdx(i,j)*v(i,j)-u(i,j)*dXdy(i,j))/(dXdx(i,j)*dYdy(i,j)-dYdx(i,j)*dXdy(i,j));
    end
end

for i = 1:ntheta
    for j = 1:nr
        VTOT(i,j) = sqrt(U(i,j)^2+V(i,j)^2);
        CP(i,j) = 1-VTOT(i,j)^2/vinf^2;
        P(i,j) = 0.5*CP(i,j)*rho*vinf^2+pinf;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1
     for i = 1:ntheta-1
        DX(i,j) = X(i+1,j)-X(i,j);
        DY(i,j) = Y(i+1,j)-Y(i,j);
        DS(i,j) = sqrt(DX(i,j)^2+DY(i,j)^2);
        dx(i,j) = x(i+1,j)-x(i,j);
        dy(i,j) = y(i+1,j)-y(i,j);
        ds(i,j) = sqrt(dx(i,j)^2+dy(i,j)^2);
    end
    for i = ntheta
        DX(i,j) = X(1,j)-X(i,j);
        DY(i,j) = Y(1,j)-Y(i,j);
        DS(i,j) = sqrt(DX(i,j)^2+DY(i,j)^2);
        dx(i,j) = x(1,j)-x(i,j);
        dy(i,j) = y(1,j)-y(i,j);
        ds(i,j) = sqrt(dx(i,j)^2+dy(i,j)^2);
    end
end

for j = 1
    for i = 1:ntheta
        F(i,j) = -P(i,j)*DS(i,j);
        f(i,j) = -p(i,j)*ds(i,j);
    end
end

for j = 1
    for i = 1:ntheta
        FX(i,j) = F(i,j)*(DY(i,j)/DS(i,j));
        FY(i,j) = -F(i,j)*(DX(i,j)/DS(i,j));
        fx(i,j) = f(i,j)*(dy(i,j)/ds(i,j));
        fy(i,j) = -f(i,j)*(dx(i,j)/ds(i,j));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:ntheta
    CLSUM(i) = -rmin*CP(i,1)*sin(theta(i))*dtheta;
end

CL = 0;

for i = 1:ntheta
    CL = CL+CLSUM(i);
end

LIFT = 0.5*CL*rho*vinf^2;

for i = 1:ntheta
    M0MSUM(i,1) = -FX(i,1)*(Y(i,1)-Y(round(ntheta/2),1))+FY(i,1)*(X(i,1)-X(round(ntheta/2),1));
end

MOM = 0;

for i = 1:ntheta
    MOM = MOM+M0MSUM(i,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(x(:,1),y(:,1),'-k','linewidth',2)
hold on
fill(x(:,1),y(:,1),'r')
quiver(x,y,u,v,1)
title('Velocity Vectors Around Cylinder')
axis([-3,3,-3,3])
xlabel('x-axis: u flow direction')
ylabel('y-axis: v flow direction')
daspect([1 1 1])
text(0,2.7,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
hold off

figure(2)
contour(x,y,vtot,200)
hold on
plot(x(:,1),y(:,1),'-k','linewidth',2)
fill(x(:,1),y(:,1),'r')
title('Velocity Contours Around Cylinder')
axis([-3,3,-3,3])
xlabel('x-axis: u flow direction')
ylabel('y-axis: v flow direction')
daspect([1 1 1])
text(0,2.7,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
hold off

figure(3)
contour(x,y,cp,200)
hold on
plot(x(:,1),y(:,1),'-k','linewidth',2)
fill(x(:,1),y(:,1),'r')
title('C_p Surface Contours Around Cylinder')
axis([-3,3,-3,3])
xlabel('x-axis: U flow direction')
ylabel('y-axis: V flow direction')
daspect([1 1 1])
text(0,2.7,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
hold off

% figure(4)
% plot(x(:,1),-cp(:,1),'linewidth',2)
% title('C_p Distribution Versus X For Cylinder')
% xlabel('x-axis: X')
% ylabel('y-axis: CP')
% axis([-1,1,-1,6])
% text(0,5.5,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
% grid on
%
% figure(5)
% quiver(x(:,1),y(:,1),fx(:,1),fy(:,1))
% hold on
% plot(x(:,1),y(:,1),'-k','linewidth',2)
% title('Forces on Cylinder Surface')
% axis([-1.5,1.5,-1.5,1.5])
% xlabel('x-axis: u flow direction')
% ylabel('y-axis: v flow direction')
% text(0,1.35,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
% daspect([1 1 1])
% grid on
% hold off
%
% figure(6)
% plot(X(:,1),Y(:,1),'-k','linewidth',2)
% hold on
% fill(X(:,1),Y(:,1),'-r')
% quiver(X,Y,U,V,1)
% title('Velovity Vectors Around Airfoil')
% axis([-3,3,-2,2])
% xlabel('x-axis: U flow direction')
% ylabel('y-axis: V flow direction')
% text(0,1.7,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
% daspect([1 1 1])
% hold off
%
% figure(7)
% contour(X,Y,VTOT,200)
% hold on
% plot(X(:,1),Y(:,1),'-k','linewidth',2)
% fill(X(:,1),Y(:,1),'-r')
% title('Velocity Contours Around Airfoil')
% axis([-3,3,-2,2])
% xlabel('x-axis: u flow direction')
% ylabel('y-axis: v flow direction')
% text(0,1.7,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
% daspect([1 1 1])
% hold off
%
% figure(8)
% contour(X,Y,CP,200)
% hold on
% plot(X(:,1),Y(:,1),'-k','linewidth',2)
% fill(X(:,1),Y(:,1),'-r')
% title('C_p Surface Contours Around Airfoil')
% axis([-3,3,-2,2])
% xlabel('x-axis: U flow direction')
% ylabel('y-axis: V flow direction')
% text(0,1.7,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
% daspect([1 1 1])
% hold off
%
% figure(9)
% plot(X(:,1),-CP(:,1),'linewidth',2)
% title('C_p Distribution Versus X For Airfoil')
% xlabel('x-axis: X')
% ylabel('y-axis: CP')
% axis([-2,2,-1,3])
% text(0,2.75,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
% grid on
%
% figure(10)
% quiver(X(:,1),Y(:,1),FX(:,1),FY(:,1))
% hold on
% plot(X(:,1),Y(:,1),'-k','linewidth',2)
% title('Forces on Airfoil Surface')
% axis([-3,3,-1,1])
% xlabel('x-axis: U flow direction')
% ylabel('y-axis: V flow direction')
% text(0,-0.7,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
% grid on
% daspect([1 1 1])
% hold off

figure(11)
plot(erru)

figure(12)
plot(errv,'--k')
