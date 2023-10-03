%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 1 - Problem 4 - Theta Separation versus Gamma
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section sets up coordinates

vinf = 1;
R = 1;
gamma = linspace(0,4*pi,20);

for i = 1:20
    theta1(i) = asin(-gamma(i)/(4*pi*vinf*R));
    theta2(i) = pi+theta1(i);
end

figure(1)
plot(gamma,theta1,'--')
hold on
plot(gamma,theta2)
title('Theta Separation versus Gamma')
xlabel('Gamma')
ylabel('Theta Separation')
hold off
