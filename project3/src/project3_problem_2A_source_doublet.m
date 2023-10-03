%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 3 - Problem 2A_1 - 3D SOURCE/DOUBLET
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q = 1;
D = 1;
rho = 1;
min = -4;
max = 4;
n = 10;

x = linspace(min, max, n);
y = linspace(min, max, n);
z = linspace(min, max, n);

[X, Y, Z] = meshgrid(x, y, z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOURCE

for i = 1:n
    for j = 1:n
        for k = 1:n
            vx(i, j, k) = ((Q/(4*pi*rho))*X(i,j,k)/(X(i,j,k)^2+Y(i,j,k)^2+Z(i,j,k)^2)^(3/2));
            vy(i, j, k) = ((Q/(4*pi*rho))*Y(i,j,k)/(X(i,j,k)^2+Y(i,j,k)^2+Z(i,j,k)^2)^(3/2));
            vz(i, j, k) = ((Q/(4*pi*rho))*Z(i,j,k)/(X(i,j,k)^2+Y(i,j,k)^2+Z(i,j,k)^2)^(3/2));
        end
    end
end

figure(1)
quiver3(X, Y, Z, vx, vy, vz, 3)
title('Source')
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOUBLET

for i = 1:n
    for j = 1:n
        for k = 1:n
            ux(i, j, k) = -2*D*X(i, j, k)*Z(i,j,k)/(X(i,j,k)^2+Y(i,j,k)^2+Z(i,j,k)^2)^(5/2);
            vr(i, j, k) = -3*D*Z(i, j, k)*sqrt(Z(i,j,k)^2+Y(i,j,k)^2)/(X(i,j,k)^2+Y(i,j,k)^2+Z(i,j,k)^2)^(5/2)+D*(Z(i,j,k)/sqrt(Y(i,j,k)^2+Z(i,j,k)^2))+1/(X(i,j,k)^2+Y(i,j,k)^2+Z(i,j,k)^2)^(3/2);
            wT(i, j, k) = D*(Y(i, j, k)/sqrt(Z(i,j,k)^2+Y(i,j,k)^2))/(X(i,j,k)^2+Y(i,j,k)^2+Z(i,j,k)^2)^(3/2);
        end
    end
end

for i = 1:n
    for j = 1:n
        for k = 1:n
            vx(i, j, k) = ux(i, j, k);
            vz(i, j, k) = vr(i, j, k)*Z(i,j,k)/sqrt(X(i,j,k)^2+Y(i,j,k)^2)-wT(i,j,k)*Y(i,j,k)/sqrt(X(i,j,k)^2+Y(i,j,k)^2)+1;
            vy(i, j, k) = vr(i, j, k)*Y(i,j,k)/sqrt(X(i,j,k)^2+Y(i,j,k)^2)+wT(i,j,k)*Z(i,j,k)/sqrt(X(i,j,k)^2+Y(i,j,k)^2);
        end
    end
end

figure(2)
quiver3(X, Y, Z, vx, vy, vz, 3)
title('Sphere')
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
