%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 2 - Problem 2 - Cp, Velocity Vectors, Etc.
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section sets up geometric parameters of the cylinder, flow parameters and grid spacing

nr = 50;
rmin = 1;
rmax = 7;
r = linspace(rmin,rmax,nr);
dr = (rmax-rmin)/(nr-1);

ntheta = 4*20;
dtheta = 2*pi/ntheta;
thetamin = 0.5*dtheta;
thetamax = 2*pi-dtheta;
theta = linspace(thetamin,thetamax,ntheta);

vinf = 1;
pinf = 0;
rho = 1;
alphadeg = 5;
alpha = alphadeg*pi/180;
a = rmin;

epsilon = -0.1;
mu = 0.1;
b = sqrt(a^2-mu^2)+epsilon;

% epsilon = 0;
% mu = 0;
% b = 0.9; %use this b with mu and epsilon zero to make ellipse with b<a

D = vinf*rmin^2;
thetasep = -alpha-asin(mu/a);
% Gamma = 4*pi*vinf*rmin*sin(thetasep);
Gamma = -4*pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:ntheta
    for j = 1:nr
        ubar(i,j) = vinf*cos(theta(i)-alpha)-D*cos(theta(i)-alpha)/(r(j)^2);
        vbar(i,j) = -vinf*sin(theta(i)-alpha)-D*sin(theta(i)-alpha)/(r(j)^2)+Gamma/(2*pi*r(j));
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
        CF(i,j) = -CP(i,j)*DS(i,j);
        cf(i,j) = -cp(i,j)*ds(i,j);
    end
end

for j = 1
    for i = 1:ntheta
        FX(i,j) = F(i,j)*(DY(i,j)/DS(i,j));
        FY(i,j) = -F(i,j)*(DX(i,j)/DS(i,j));
        fx(i,j) = f(i,j)*(dy(i,j)/ds(i,j));
        fy(i,j) = -f(i,j)*(dx(i,j)/ds(i,j));
        CFX(i,j) = CF(i,j)*(DY(i,j)/DS(i,j));
        CFY(i,j) = -CF(i,j)*(DX(i,j)/DS(i,j));
        cfx(i,j) = cf(i,j)*(dy(i,j)/ds(i,j));
        cfy(i,j) = -cf(i,j)*(dx(i,j)/ds(i,j));
    end
end

% FY(42) = FY(43);
% FX(42) = FX(43);
% FX(80) = FX(79);
% FY(80) = FY(79);
% U(1,1) = U(2,1);
% V(1,1) = V(2,1);
% U(80,1) = U(79,1);
% V(80,1) = V(79,1);
% CP(42,1) = CP(43,1);
% CP(78,1) = CP(77,1);
% CP(79,1) = CP(78,1);
% CP(80,1) = CP(79,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:ntheta
    CLSUM(i) = -rmin*CP(i,1)*sin(theta(i))*dtheta;
    % CDSUM(i) = -0.5*rmin*CP(i,1)*cos(theta(i))*dtheta;
    CDSUM(i) = -0.5*rmin*CP(i,1)*(X(i,1)/(sqrt(X(i,1)^2+Y(i,1)^2)))*dtheta;
end

CL = 0;
CD = 0;

for i = 1:ntheta
    CL = CL+CLSUM(i);
    CD = CD+CDSUM(i);
end

LIFT = 0.5*CL*rho*vinf^2;

for i = 1:ntheta
    M0MSUM(i,1) = -FX(i,1)*(Y(i,1)-Y(round(ntheta/2),1))+FY(i,1)*(X(i,1)-X(round(ntheta/2),1));
    CMSUM(i,1) = -CFX(i,1)*(Y(i,1)-Y(round(ntheta/2),1))+CFY(i,1)*(X(i,1)-X(round(ntheta/2),1));
end

MOM = 0;
CM = 0;

for i = 1:ntheta
    MOM = MOM+M0MSUM(i,1);
    CM = CM+CMSUM(i,1);
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
colorbar
text(0,2.7,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
hold off

figure(3)
contour(x,y,cp,200)
hold on
plot(x(:,1),y(:,1),'-k','linewidth',2)
fill(x(:,1),y(:,1),'r')
title('C_p Surface Contours On Cylinder')
axis([-3,3,-3,3])
xlabel('x-axis: U flow direction')
ylabel('y-axis: V flow direction')
daspect([1 1 1])
colorbar
text(0,2.7,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
hold off

figure(4)
plot(x(:,1),-cp(:,1),'linewidth',2)
title('C_p Distribution Versus X For Cylinder')
xlabel('x-axis: X')
ylabel('y-axis: CP')
axis([-1,1,-1,6])
text(0,5.5,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
grid on

figure(5)
quiver(x(:,1),y(:,1),fx(:,1),fy(:,1))
hold on
plot(x(:,1),y(:,1),'-k','linewidth',2)
title('Forces on Cylinder Surface')
axis([-1.5,1.5,-1.5,1.5])
xlabel('x-axis: u flow direction')
ylabel('y-axis: v flow direction')
text(0,1.35,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
daspect([1 1 1])
grid on
hold off

figure(6)
plot(X(:,1),Y(:,1),'-k','linewidth',2)
hold on
fill(X(:,1),Y(:,1),'-r')
quiver(X,Y,U,V,1)
title('Velocity Vectors Around Airfoil')
axis([-3,3,-2,2])
xlabel('x-axis: U flow direction')
ylabel('y-axis: V flow direction')
text(0,1.7,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
daspect([1 1 1])
hold off

figure(7)
contour(X,Y,VTOT,2000)
hold on
plot(X(:,1),Y(:,1),'-k','linewidth',2)
fill(X(:,1),Y(:,1),'-r')
title('Velocity Contours Around Airfoil')
axis([-3,3,-2,2])
xlabel('x-axis: u flow direction')
ylabel('y-axis: v flow direction')
text(0,1.7,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
daspect([1 1 1])
colorbar
hold off

figure(8)
contour(X,Y,CP,2000)
hold on
plot(X(:,1),Y(:,1),'-k','linewidth',2)
fill(X(:,1),Y(:,1),'-r')
title('C_p Surface Contours On Airfoil')
axis([-3,3,-2,2])
xlabel('x-axis: U flow direction')
ylabel('y-axis: V flow direction')
text(0,1.7,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
daspect([1 1 1])
colorbar
hold off

figure(9)
plot(X(:,1),-CP(:,1),'linewidth',2)
title('C_p Distribution Versus X For Airfoil')
xlabel('x-axis: X')
ylabel('y-axis: CP')
axis([-2,2,-2,3])
text(0,2.7,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
grid on

figure(10)
quiver(X(:,1),Y(:,1),FX(:,1),FY(:,1))
hold on
plot(X(:,1),Y(:,1),'-k','linewidth',2)
title('Forces on Airfoil Surface')
axis([-3,3,-1,1])
xlabel('x-axis: U flow direction')
ylabel('y-axis: V flow direction')
text(0,-0.7,['\alpha  =  ',sprintf('%0.2f',alphadeg),'^\circ, v_\infty =  ',sprintf('%0.2f',vinf),', \rho_\infty =  ',sprintf('%0.2f',rho),', a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
grid on
daspect([1 1 1])
hold off
