%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 4 - Problem 5 - Flow Over Slender Delta Wing
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uinf = 1;
alpha_deg = 5;
alpha = alpha_deg*pi/180;

b = 5;
chord(1) = 8;
ny = 21;
nx = 21;
nt = nx-1;
y = linspace(0,b/2,ny);
m = b/chord(1);

for i = 1:ny
    xTE(i) = b/(2*m)-y(i)/m;
    xLE(i) = -xTE(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:ny
    x = linspace(xLE(i),xTE(i),ny);
    dx(i) = 2*xTE(i)/(nt);
    t = linspace(xLE(i)+dx(i)/2,xTE(i)-dx(i)/2,nt);

    A = zeros(nx,nx);

    A(1,1) = 1;
    A(nx,nx) = 1;

        A(i,j) = 
    for j = 2:nt-1
        A(i,j) = -1/(y(i-1)-t(j-1))-1/(y(i-1)-t(j));
    end

    % TODO@dpwiese - this is wrong
    A(i,i) = 0

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(xTE,y,'--','linewidth',2)
hold on
plot(xLE,y,'-','linewidth',2)
daspect([1 1 1])
legend('Trailing Edge','Leading Edge')
grid on
hold off
