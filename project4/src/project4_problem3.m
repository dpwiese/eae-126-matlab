%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAE 126 Computational Aerodynamics (Spring 2011)
% Project 4 - Problem 3 - Swept Wings
% Daniel Wiese

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for alpha_deg = 1:2

    b = 5;
    chord(1) = .2*b;
    LambdaLEdeg = -45;
    LambdaLE = LambdaLEdeg*pi/180;
    LambdaTEdeg = -45;
    LambdaTE = LambdaTEdeg*pi/180;
    uinf = 1;
    rhoinf = 1;
    alphadeg(alpha_deg) = alpha_deg-1;
    alpha = alphadeg(alpha_deg)*pi/180;

    ny = 401;
    np = ny-1;
    dy = (b/2)/(ny-1);
    y = linspace(0,b/2,ny);
    p = linspace(dy/2,(b-dy)/2,np);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:ny
        xnLE(i) = y(i)*tan(LambdaLE);
    end

    for i = 1:ny
        xnTE(i) = chord(1)+y(i)*tan(LambdaTE);
    end

    for i = 1:ny
        chord(i) = xnTE(i)-xnLE(i);
    end

    for i = 1:ny
        xn(i) = xnLE(i)+0.25*chord(i);
    end

    for i = 1:np
        xm(i) = ((xnLE(i)+xnLE(i+1))/2)+0.75*((chord(i)+chord(i+1))/2);
    end

    for i = 1:np
        ym(i) = y(i)+dy/2;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:np
        B(i) = -uinf*alpha;
    end

    % Right Side aka Starboard
    for i = 1:np
        for j = 1:np
            temp1 = 1/((xm(i)-xn(j))*(ym(i)-y(j+1))-(xm(i)-xn(j+1))*(ym(i)-y(j)));
            temp2 = ((xn(j+1)-xn(j))*(xm(i)-xn(j))+(y(j+1)-y(j))*(ym(i)-y(j)))/sqrt((xm(i)-xn(j))^2+(ym(i)-y(j))^2);
            temp3 = -((xn(j+1)-xn(j))*(xm(i)-xn(j+1))+(y(j+1)-y(j))*(ym(i)-y(j+1)))/sqrt((xm(i)-xn(j+1))^2+(ym(i)-y(j+1))^2);
            temp4 = (1/(y(j)-ym(i)))*(1+(xm(i)-xn(j))/sqrt((xm(i)-xn(j))^2+(ym(i)-y(j))^2));
            temp5 = -(1/(y(j+1)-ym(i)))*(1+(xm(i)-xn(j+1))/sqrt((xm(i)-xn(j+1))^2+(ym(i)-y(j+1))^2));
            wmnR(i,j) = (1/(4*pi))*(temp1*(temp2+temp3)+temp4+temp5);
        end
    end

    % Left Side aka Port
    for i = 1:np
        for j = 1:np
            temp1 = 1/((xm(i)-xn(j))*(ym(i)+y(j+1))-(xm(i)-xn(j+1))*(ym(i)+y(j)));
            temp2 = ((xn(j+1)-xn(j))*(xm(i)-xn(j))+(-y(j+1)+y(j))*(ym(i)+y(j)))/sqrt((xm(i)-xn(j))^2+(ym(i)+y(j))^2);
            temp3 = -((xn(j+1)-xn(j))*(xm(i)-xn(j+1))+(-y(j+1)+y(j))*(ym(i)+y(j+1)))/sqrt((xm(i)-xn(j+1))^2+(ym(i)+y(j+1))^2);
            temp4 = (1/(-y(j)-ym(i)))*(1+(xm(i)-xn(j))/sqrt((xm(i)-xn(j))^2+(ym(i)+y(j))^2));
            temp5 = -(1/(-y(j+1)-ym(i)))*(1+(xm(i)-xn(j+1))/sqrt((xm(i)-xn(j+1))^2+(ym(i)+y(j+1))^2));
            wmnL(i,j) = -(1/(4*pi))*(temp1*(temp2+temp3)+temp4+temp5);
        end
    end

    for i = 1:np
        for j = 1:np
            A(i,j) = wmnR(i,j)+wmnL(i,j);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform Gaussian Elimination

    n = np;

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

    for j = 1:np
        L = L+2*rhoinf*uinf*Gamma(j)*dy;
    end

    S = (chord(1)+chord(ny))*y(ny);
    CL(alpha_deg) = L/(S*0.5*rhoinf*uinf^2);

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
plot(y,xnTE,'--','linewidth',2)
plot(-y,xnLE,'-','linewidth',2)
plot(-y,xnTE,'--','linewidth',2)
daspect([1 1 1])
title('Wing Planform')
xlabel('y-axis: Spanwise Direction')
ylabel('x-axis: Chordwise Direction')
legend('Wing Leading Edge','Wing Trailing Edge')
grid on
% axis([-y(ny) y(ny) -min(xnLE) max(abs(xnTE))])
text(0,0.75*max(abs(xnTE)),[sprintf('%0.2f',np),' Panels'],'HorizontalAlignment','center','BackgroundColor','w','Edgecolor','k');
axis ij
hold off

figure(3)
plot(alphadeg,CL,'-','linewidth',2)
title('C_L versus \alpha')
grid on
xlabel('\alpha (\circ)')
ylabel('C_L')
