%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 4 - Problem 1 - Elliptic Wing
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 6;
% AR = 12;          % For elliptic only
% b = 4*a/(pi*AR);  % For elliptic only
b = 0.5;
tau = 0.10;
uinf = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ymin = -a;
ymax = a;
ny = 200;
dy = (ymax-ymin)/(ny-1);
y = linspace(ymin,ymax,ny);

xmin = -b;
xmax = b;
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
    % chord(i) = 2*(b*sqrt(1-y(i)^2/a^2));  % Ellipse
    chord(i) = 2*b;                         % Rectangle
    aeff(i) = chord(i)/2;
end

for i = 1:ny
    for j = 1:nx
        z(i,j) = -(tau/aeff(i))*x(j)^2+tau*aeff(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:ny
    for j = 1:nx
        u(i,j) = 0;
        for k = 1:ntx
            for m = 1:nty
                u(i,j) = u(i,j)+((tau*uinf)/(2*pi))*(((z(m,k+1)-z(m,k))+(z(m+1,k+1)-z(m+1,k)))/(2*dx))*((x(j)-tx(k))*dx*dy)/((x(j)-tx(k))^2+(y(i)-ty(m))^2)^(3/2);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(1,2,1)
plot(y,chord/2,'linewidth',2)
hold on
grid on
plot(y,-chord/2,'linewidth',2)
axis([-a-1 a+1 -b-1 b+1])
axis square
xlabel('y-axis')
ylabel('x-axis')
title('Wing Planform')
hold off

subplot(1,2,2)
plot(x,z(10,:))
hold on
grid on
plot(x,-z(10,:))
axis([-b-1 b+1 -1 1])
axis square
xlabel('y-axis')
ylabel('x-axis')
title('Wing Cross Section')
hold off

figure(2)
subplot(2,1,1)
plot(x,u(1,:))
hold on
plot(x,u(3,:))
plot(x,u(5,:))
plot(x,u(7,:))
plot(x,u(9,:))
plot(x,u(11,:))
plot(x,u(13,:))
plot(x,u(15,:))
plot(x,u(17,:))
plot(x,u(19,:))
title('U_T On Upper Surface of Wing')
xlabel('x-axis: Chordwise Direction')
ylabel('U_T')
grid on
hold off

subplot(2,1,2)
plot(y,u(:,nx/2))
title('U_T_,_m_a_x On Upper Surface of Wing')
xlabel('y-axis: Spanwise Direction')
ylabel('U_T')
grid on
