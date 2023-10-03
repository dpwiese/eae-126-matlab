%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 2 - Problem 1 - Joukowski Transformations
% Daniel Wiese
% Page 66 Moran

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 1;
ntheta = 50;
theta = linspace(0,2*pi,ntheta);

epsilon = 0;
mu = 0;
b = sqrt(a^2-mu^2)+epsilon;

% Use this b with mu and epsilon zero to make ellipse with b < a
b = 0.9;

for i = 1:ntheta
    x(i) = epsilon+a*cos(theta(i));
    y(i) = mu+a*sin(theta(i));
    X(i) = x(i)*(1+(b^2)/(x(i)^2+y(i)^2));
    Y(i) = y(i)*(1-(b^2)/(x(i)^2+y(i)^2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(1,2,1)
plot(x,y,'-k','linewidth',2)
fill(x,y,'r')
title('Circle Before Transformation')
text(0,1.7,['a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
axis([-2,2,-2,2])
grid on
daspect([1 1 1])

subplot(1,2,2)
plot(X,Y,'-k','linewidth',2)
fill(X,Y,'r')
title('Circle After Transformation: Airfoil')
text(0,2.7,['a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
axis([-3,3,-3,3])
grid on
daspect([1 1 1])

figure(2)
plot(X,Y,'-k','linewidth',2)
%fill(X,Y,'r')
title('Joukowski Airfoil','fontsize',20)
text(0,-0.5,['a =  ',sprintf('%0.2f',a),', \epsilon =  ',sprintf('%0.2f',epsilon),', \mu =  ',sprintf('%0.2f',mu),', b =  ',sprintf('%0.2f',b)],'fontsize',14,'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
axis([-2,2,-0.6,0.6])
grid on
daspect([1 1 1])
