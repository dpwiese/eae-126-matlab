%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 3 - Problem 1B - LIFTING PROBLEM
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uinf = 1;
rhoinf = 1;
a = 1;
tau = 0.1;
b = tau*a;
alphadeg = 0;
alpha = alphadeg*pi/180;

xmin = -a;
xmax = a;
nx = 200;
dx = (xmax-xmin)/(nx-1);
x = linspace(xmin,xmax,nx);

tmin = xmin+dx/2;
tmax = xmax-dx/2;
nt = nx-1;
t = linspace(tmin,tmax,nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nx
    y(i) = -(0.1/a)*x(i)^2+0.1*a;   % Circular arc
    % y(i) = 0;                     % Flat plate
end

for i = 1
    dydx(i) = (y(i+1)-y(i))/dx;
end

for i = 2:nx-1
    dydx(i) = (y(i+1)-y(i-1))/(2*dx);
end

for i = nx
    dydx(i) = (y(i)-y(i-1))/dx;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = nx-1;

for i = 1:nx-2
    for j = 1:nt
        A(i,j) = dx/(x(i+1)-t(j));
    end
    B(i) = (alpha-dydx(i))*2*pi*uinf;
end

for i = nx-1
    for j = 1:nt-2
        A(i,j) = 0;
    end
        A(i,nt-1) = -1/2;
        A(i,nt) = 2/2;
        B(i) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Gaussian Elimination

% Column
for j = 1:n-1
    % Rows
    for i = j:n-1
        bar  =  A(i+1,j)/A(j,j);
        A(i+1,:)  =  A(i+1,:) - A(j,:)*bar;
        B(i+1)  =  B(i+1) - bar*B(j);
    end
end

% Perform back substitution
gamma  =  zeros(n,1);
gamma(n)  =  B(n)/A(n,n);

for j = n-1:-1:1,
    gamma(j)  =  (B(j)-A(j,j+1:n)*gamma(j+1:n))/A(j,j);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gammax(1) = gamma(1)-0.5*(gamma(2)-gamma(1));
for i = 2:nx-1
    gammax(i) = (gamma(i)+gamma(i-1))/2;
end
gammax(nx) = 0;

gammasum = 0;
for i = 1:nx-1
    gammasum = gammasum+(gammax(i)+gammax(i+1))*0.5*dx;
end
GAMMA = gammasum;

numer = 0;
denom = 0;
for i = 1:nx-1
    numer = numer+(x(i)-xmin)*(gammax(i)+gammax(i+1))*0.5*dx;
    denom = denom+(gammax(i)+gammax(i+1))*0.5*dx;
end
xcp = 0.5*numer/denom;

LIFT = rhoinf*uinf*GAMMA;
CL = LIFT/(0.5*rhoinf*uinf^2);

for i = 1:nx
    uc_x(i) = 0.5*gammax(i);
end

for i = 1:nx
    cp(i) = 2*uc_x(i)/uinf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(1,2,1)
plot(x,y,'linewidth',2)
hold on
daspect([1 1 1])
axis([-a, a, -a, a])
title('Profile Shape')
grid on
xlabel('x-axis')
ylabel('y-axis')
hold off

subplot(1,2,2)
plot(x,cp)
title('C_p Distribution Over Profile Surface')
xlabel('Chord Position')
ylabel('-C_p')
grid on
