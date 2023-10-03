%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 4 - Problem 1 - Rectangular Wing
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% b is span
b = 12;
AR = 12;
cbar = b/AR;
tau = 0.10;
uinf = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ymin = -b/2;
ymax = b/2;
ny = 200;
dy = (ymax-ymin)/(ny-1);
y = linspace(ymin,ymax,ny);

xmin = -cbar/2;
xmax = cbar/2;
nx = 20;
dx = (xmax-xmin)/(nx-1);
x = linspace(xmin,xmax,nx);

ntx = nx-1;
tx = linspace(xmin+dx/2,xmax-dx/2,ntx);
dtx = dx;

nty = ny-1;
ty = linspace(ymin+dy/2,ymax-dx/2,nty);
dty = dy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:ny
    chord(i) = cbar;
end

for i = 1:ny
    for j = 1:nx
        z(i,j) = -(tau/(cbar/2))*x(j)^2+tau*(cbar/2);
    end
end

for i = 1:ntx
    dzdx(i) = -4*tau*tx(i)/cbar;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:ny
    for j = 1:nx
        u(i,j) = 0;
        for k = 1:ntx
            for m = 1:nty
                u(i,j) = u(i,j)+((uinf)/(2*pi))*(((z(m,k+1)-z(m,k))+(z(m+1,k+1)-z(m+1,k)))/(2*dx))*((x(j)-tx(k))*dx*dy)/((x(j)-tx(k))^2+(y(i)-ty(m))^2)^(3/2); %numerical derivative
                % u(i,j) = u(i,j)+((uinf)/(2*pi))*dzdx(k)*((x(j)-tx(k))*dx*dy)/((x(j)-tx(k))^2+(y(i)-ty(m))^2)^(3/2); %analytical derivative
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(2,1,1)
plot(y,chord/2,'linewidth',2)
hold on
grid on
plot(y,-chord/2,'linewidth',2)
axis([-b/2 b/2 -cbar/2-1 cbar/2+1])
axis square
xlabel('y-axis: Spanwise Direction')
ylabel('x-axis: Chordwise Direction')
title('Wing Planform')
daspect([1 1 1])
hold off

subplot(2,1,2)
plot(y,u(:,nx/2),'linewidth',2)
title('u_T_,_m_a_x On Upper Surface of Wing')
xlabel('y-axis: Spanwise Direction')
ylabel('u_T')
axis([-b/2 b/2 0 1.1*max(u(:,nx/2))])
grid on

figure(2)
subplot(2,1,1)
plot(x,z(nx/2,:),'linewidth',2)
hold on
grid on
plot(x,-z(nx/2,:),'linewidth',2)
axis([-1.0*cbar/2 1.0*cbar/2 -1.1*tau*cbar/2 1.1*tau*cbar/2])
xlabel('x-axis: Chordwise Direction')
ylabel('z-axis: Vertical Direction')
title('Wing Profile Cross Section')
daspect([1 1 1])
hold off

subplot(2,1,2)
hold on
for i = 1:nx
    plot(x,u(i,:))
end
axis([-1.0*cbar/2 1.0*cbar/2 -1.5*max(u(:,nx/2)) 1.5*max(u(:,nx/2))])
title('u_T On Upper Surface of Wing')
xlabel('x-axis: Chordwise Direction')
ylabel('u_T')
grid on
hold off
