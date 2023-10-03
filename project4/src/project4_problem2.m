%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 4 - Problem 2 - Flow Over an Elliptic Wing
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 18;
AR = 12;
b = (4*a)/(pi*AR);

uinf = 1;
alphadeg = 5;
alpha = alphadeg*pi/180;

ymin = -a;
ymax = a;
n = 100;
dy = (ymax-ymin)/(n-1);
tmin = ymin+dy/2;
tmax = ymax-dy/2;
nt = n-1;
y = linspace(ymin,ymax,n);
t = linspace(tmin,tmax,nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:n
    xplan(i) = (b*sqrt(1-y(i)^2/a^2));
    chord(i) = 2*xplan(i);
end

A = zeros(n,n);

for i = 1:n-1
    for j = 1
        A(i,j) = 1;
    end
end

for i = 2:n
    for j = n
        A(i,j) = 1;
    end
end

for i = 2:n-1
    for j = 2:n-1
            A(i,j) = (1/(4*uinf*pi))*((1/(y(i)-t(j-1)))-(1/(y(i)-t(j))));
    end
end

for i = 2:n-1
    A(i,i) = 1/(pi*uinf*chord(i))+(1/(4*uinf*pi))*(1/(y(i)-t(i-1))-1/(y(i)-t(i)));
end

B(1) = 0;
for i = 2:n-1
    B(i) = alpha;
end
B(n) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Gaussian Elimination

% Columns
for j = 1:n-1
    % Rows
    for i = j:n-1
        bar = A(i+1,j)/A(j,j);
        A(i+1,:) = A(i+1,:)- A(j,:)*bar;
        B(i+1)  =  B(i+1) - bar*B(j);
    end
end

% Perform back substitution
xsol = zeros(n,1);
xsol(n) = B(n)/A(n,n);

for j = n-1:-1:1,
    xsol(j) = (B(j)-A(j,j+1:n)*xsol(j+1:n))/A(j,j);
end

Gamma = xsol;
GammaS = (8*b*pi*uinf*a*alpha)/(2*b*pi+4*a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(2,1,1)
plot(y,xplan,'linewidth',2)
hold on
plot(y,-xplan,'linewidth',2)
grid on
axis([ymin ymax -2*b 2*b])
title('Wing Planform')
xlabel('y-axis: Spanwise Direction')
ylabel('x-axis: Chordwise Direction')
daspect([1 1 1])
hold off

subplot(2,1,2)
plot(y,Gamma,'linewidth',2)
axis([ymin ymax 0 1.1*GammaS])
title('Spanwise Circulation Distribution')
xlabel('y-axis: Spanwise Direction')
ylabel('Circulation \Gamma(y)')
text(0,0.2,['\Gamma_a_n_a_l_y_t_i_c_a_l =  ',sprintf('%0.2f',GammaS)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
grid on
