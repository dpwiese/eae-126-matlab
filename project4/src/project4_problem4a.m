%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 4 - Problem 4a - Elliptic Crescent Wing
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for alpha_deg = 1:11

    % b is span
    b = 5;
    ae1 = b/2;
    be1 = 0.2*b;
    ae2 = 0.95*ae1;
    be2 = 0.9*be1;
    uinf = 1;
    rhoinf = 1;
    alpha_rad = alpha_deg*pi/180;

    ny = 401;
    num_p = ny-1;
    dy = (b/2)/(ny-1);
    y = linspace(0,b/2,ny);
    sh = 0.5;
    p = linspace(dy/2,(b-dy)/2,num_p);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:ny
        xnLE(i) = -be1*sqrt(1-(y(i)/ae1)^2);
        xnLE2(i) = sh-be2*sqrt(1-(y(i)/ae2)^2);
        xnTE2(i) = sh+be2*sqrt(1-(y(i)/ae2)^2);
    end

    for i = 1:ny
        xnTE(i) = be1*sqrt(1-(y(i)/ae1)^2);
        if xnTE(i)>xnLE2(i)
            xnTE(i) = sh-be2*sqrt(1-(y(i)/ae2)^2);
        end
    end

    for i = 1:ny
        chord(i) = xnTE(i)-xnLE(i);
    end

    for i = 1:ny
        xn(i) = xnLE(i)+0.25*chord(i);
    end

    for i = 1:num_p
        xm(i) = (xnTE(i)+xnTE(i+1))/2;
    end

    for i = 1:num_p
        ym(i) = y(i)+dy/2;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:num_p
        B(i) = -uinf*alpha_rad;
    end

    % Right Side aka Starboard
    for i = 1:num_p
        for j = 1:num_p
            temp1 = 1/((xm(i)-xn(j))*(ym(i)-y(j+1))-(xm(i)-xn(j+1))*(ym(i)-y(j)));
            temp2 = ((xn(j+1)-xn(j))*(xm(i)-xn(j))+(y(j+1)-y(j))*(ym(i)-y(j)))/sqrt((xm(i)-xn(j))^2+(ym(i)-y(j))^2);
            temp3 = -((xn(j+1)-xn(j))*(xm(i)-xn(j+1))+(y(j+1)-y(j))*(ym(i)-y(j+1)))/sqrt((xm(i)-xn(j+1))^2+(ym(i)-y(j+1))^2);
            temp4 = (1/(y(j)-ym(i)))*(1+(xm(i)-xn(j))/sqrt((xm(i)-xn(j))^2+(ym(i)-y(j))^2));
            temp5 = -(1/(y(j+1)-ym(i)))*(1+(xm(i)-xn(j+1))/sqrt((xm(i)-xn(j+1))^2+(ym(i)-y(j+1))^2));
            wmnR(i,j) = (1/(4*pi))*(temp1*(temp2+temp3)+temp4+temp5);
        end
    end

    % Left Side aka Port
    for i = 1:num_p
        for j = 1:num_p
            temp1 = 1/((xm(i)-xn(j))*(ym(i)+y(j+1))-(xm(i)-xn(j+1))*(ym(i)+y(j)));
            temp2 = ((xn(j+1)-xn(j))*(xm(i)-xn(j))+(-y(j+1)+y(j))*(ym(i)+y(j)))/sqrt((xm(i)-xn(j))^2+(ym(i)+y(j))^2);
            temp3 = -((xn(j+1)-xn(j))*(xm(i)-xn(j+1))+(-y(j+1)+y(j))*(ym(i)+y(j+1)))/sqrt((xm(i)-xn(j+1))^2+(ym(i)+y(j+1))^2);
            temp4 = (1/(-y(j)-ym(i)))*(1+(xm(i)-xn(j))/sqrt((xm(i)-xn(j))^2+(ym(i)+y(j))^2));
            temp5 = -(1/(-y(j+1)-ym(i)))*(1+(xm(i)-xn(j+1))/sqrt((xm(i)-xn(j+1))^2+(ym(i)+y(j+1))^2));
            wmnL(i,j) = -(1/(4*pi))*(temp1*(temp2+temp3)+temp4+temp5);
        end
    end

    for i = 1:num_p
        for j = 1:num_p
            A(i,j) = wmnR(i,j)+wmnL(i,j);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform Gaussian Elimination

    n = num_p;

    % Columns
    for j = 1:n-1,

        % Rows
        for i = j:n-1,
            bar = A(i+1,j)/A(j,j);
            A(i+1,:) = A(i+1,:)-A(j,:)*bar;
            B(i+1) = B(i+1)-bar*B(j);
        end
    end

    % Perform back substitution
    xsol = zeros(n,1);
    xsol(n) = B(n)/A(n,n);

    for j = n-1:-1:1,
        xsol(j) = (B(j)-A(j,j+1:n)*xsol(j+1:n))/A(j,j);
    end

    Gamma = xsol;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    L = 0;
    S = 0;

    for i = 1:num_p
        L = L+2*rhoinf*uinf*Gamma(i)*dy;
        S = S+((chord(i)+chord(i+1))/2)*dy;
    end

    S = 2*S;
    AR = b^2/S;
    CL(il) = L/(S*0.5*rhoinf*uinf^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(1)
    hold on
    plot(p,Gamma,'-','linewidth',2)
    title('\Gamma(y) Across Half Wingspan')
    xlabel('y-axis: Span Location')
    ylabel('\Gamma(y)')
    grid on
    hold off

end

figure(2)
plot(y,xnLE,'-','linewidth',2)
hold on
axis ij
plot(y,xnTE,'--','linewidth',2)
plot(-y,xnLE,'-','linewidth',2)
plot(-y,xnTE,'--','linewidth',2)
daspect([1 1 1])
axis([-1.1*ae1 1.1*ae1 -2.0*be1 2.0*be1])
title('Wing Planform')
xlabel('y-axis: Spanwise Direction')
ylabel('x-axis: Chordwise Direction')
legend('Wing Leading Edge','Wing Trailing Edge')
grid on
text(0,0.5*be1,[sprintf('%0.2f',num_p),' Panels, S =  ',sprintf('%0.3f',S),', AR =  ',sprintf('%0.1f',AR)],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
hold off

figure(3)
plot(alpha_deg,CL,'-','linewidth',2)
title('C_L versus \alpha')
grid on
xlabel('\alpha (\circ)')
ylabel('C_L')
