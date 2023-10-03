%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 4 - Problem 1 - Elliptic Wing - STRIPS
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 6;
AR = 12;
b = 4*a/(pi*AR);
tau = 0.10;
uinf = 1;

ymin = -a;
ymax = a;
ny = 100;
dy = (ymax-ymin)/ny;
y = linspace(ymin,ymax,ny);
nx = 20;
nt = nx-1;

for i = 1:ny
    x(i) = (b*sqrt(1-y(i)^2/a^2));
    chord(i) = 2*x(i);
end

for i = 1:ny
    X(i,:) = linspace(-x(i),x(i),nx);
    dx(i) = 2*x(i)/(nx-1);
    T(i,:) = linspace(-x(i)+dx(i),x(i)-dx(i),nt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:ny
    for j = 1:nx
        z(i,j) = -(tau/x(i))*X(i,j)^2+tau*x(i);
    end
end

for i = 1:ny
    for k = 1:nt
        dzdt(i,k) = (z(2,k+1)-z(2,k))/dx(2); %changed i to 2 since dzdt same across anywhere
    end
end

for i = 1:ny
    for j = 1:nx
        for k = 1:nt
            ut(i,k) = 2*abs(x(i))*(uinf/pi)*dzdt(i,k)/(X(i,j)-T(i,k));
        end
        utsum = 0;
        for m = 1:nt-1
            temp(m) = 0.5*(ut(i,m)+ut(i,m+1))*dx(i);
            utsum = utsum+temp(m);
        end
        ut_num(i,j) = utsum;
    end
end

for i = 1:ny
    for j = 1:nx
        ut_ana(i,j) = tau*uinf;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(y,chord/2,'linewidth',2)
hold on
grid on
plot(y,-chord/2,'linewidth',2)
axis([-8 8 -8 8])
axis square
xlabel('y-axis')
ylabel('x-axis')
title('Elliptic Wing')
hold off

figure(2)
plot(X(50,:),z(50,:))
hold on
grid on
plot(X(50,:),-z(50,:))
axis([-1 1 -1 1])
axis square
xlabel('y-axis')
ylabel('x-axis')
title('Cross Section')
hold off

figure(3)
plot(X(50,:),ut_num(50,:),'--')
hold on
plot(X(50,:),ut_ana(50,:))
title('U_T Across Profile Chord')
grid on
xlabel('x-axis')
ylabel('y-axis')
legend('Numerical Solution','Analytical Solution')
hold off
