%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 3 - Problem 2A_3 - 3D ELLIPSOID/BICONVEX
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uinf = 1;
a = 1;
tau = 0.1;
b = tau * a;

xmin = -a;
xmax = a;
nx = 300;
dx = (xmax-xmin)/(nx-1);
x = linspace(xmin,xmax,nx);

tmin = xmin+dx/2;
tmax = xmax-dx/2;
nt = nx-1;
t = linspace(tmin, tmax, nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nx
    y(i) = sqrt(b^2*(1-x(i)^2/a^2));    % Ellipsoid
    % y(i) = -(0.1/a)*x(i)^2+0.1*a;     % Biconvex body of revolution
end

for j = 1:nt
    dSdt(j) = ((pi*y(j+1)^2)-(pi*y(j)^2))/dx;
end

for i = 1:nx
    for j = 1:nt
        ut(j) = (a/(2*pi)) * uinf * dSdt(j) * (x(i)-t(j)) / ((x(i)-t(j))^2+y(i)^2)^(3/2);
    end

    utsum = 0;

    for k = 1:nt-1
        temp(k) = 0.5*(ut(k)+ut(k+1))*dx;
        utsum = utsum+temp(k);
    end

    ut_num(i) = utsum;
end

for i = 1:nx
    cpx(i) = 2*ut_num(i)/uinf;
end

cpx(1) = cpx(2);
cpx(nx) = cpx(nx-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(1,2,1)
plot(x,y,'linewidth',2)
hold on
plot(x,-y,'linewidth',2)
daspect([1 1 1])
axis([-a, a, -a, a])
title('Profile Shape')
grid on
xlabel('x-axis')
ylabel('y-axis')
hold off

subplot(1,2,2)
plot(x,cpx)
title('C_p Distribution Over Profile Surface')
xlabel('Chord Position')
ylabel('-C_p')
% legend('Numerical Solution','Analytical Solution')
grid on
