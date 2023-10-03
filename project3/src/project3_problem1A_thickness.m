%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 3 - Problem 1A - Thickness Problem
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uinf = 1;
a = 1;
tau = 0.1;
b = tau*a;
alpha = 0;

xmin = -a;
xmax = a;
nx = 100;
dx = (xmax-xmin)/(nx-1);
x = linspace(xmin,xmax,nx);

tmin = xmin+dx/2;
tmax = xmax-dx/2;
nt = nx-1;
t = linspace(tmin,tmax,nt);

for i = 1:nx
    % y(i) = sqrt(b^2*(1-x(i)^2/a^2));  % Ellipse
    y(i) = -(0.1/a)*x(i)^2+0.1*a;       % Biconvex
end

for j = 1:nt
    dydt(j) = (y(j+1)-y(j))/dx;
end

for i = 1:nx
    for j = 1:nt
        ut(j) = 2*a*(uinf/pi)*dydt(j)/(x(i)-t(j));
    end
    utsum = 0;
    for k = 1:nt-1
        temp(k) = 0.5*(ut(k)+ut(k+1))*dx;
        utsum = utsum+temp(k);
    end
    ut_num(i) = utsum;
end

for i = 1:nx
    ut_ana(i) = tau*uinf;
end

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
plot(x,ut_num,'--')
hold on
plot(x,ut_ana)
title('U_T Across Profile Chord')
grid on
xlabel('x-axis')
ylabel('y-axis')
legend('Numerical Solution','Analytical Solution')
hold off
