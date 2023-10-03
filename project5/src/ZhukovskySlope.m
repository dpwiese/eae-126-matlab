function [ dydxtop dydxbot ytop ybot ] = ZhukovskySlope( xvect, iLe, iTe, eps, mu )
%ZhukovskySlope
%   This function computes the slope of a Zhukovsky airfoil.
%
%   Give it your entire x-vector, the values iLe & iTe, and
%     the Zhukovsky Transformation values eps & mu.  It will kindly
%     give you the slope in the vectors dydxtop & dydxbot in the
%     range iLe:iTe

%PARAMETERS
r   = 1.;                                   %radius of cylinder
b   = sqrt(r^2-mu^2)+eps;                   %b

%GRID
theta=linspace(0,2*pi,500);
for i=1:length(theta)
    xC(i)=r.*cos(theta(i))+eps;             %Shifting the cylinder
    yC(i)=r.*sin(theta(i))+mu;
    X(i)=xC(i)*(1+b^2/(xC(i)^2+yC(i)^2));   %Zhukovsky Transformation
    Y(i)=yC(i)*(1-b^2/(xC(i)^2+yC(i)^2));
end

%SHIFT & TRANSFORM
Xmin=min(X);
Xmax=max(X);
scale=1/(Xmax-Xmin);
for i=1:length(theta)
    X(i)=(X(i)-Xmin)*scale;                 %to get airfoil in (0 <= x <= 1)
    Y(i)=Y(i)*scale;
end

%FIND SLOPE
for i=1:length(theta)-1
    dYdX(i)=(Y(i+1)-Y(i))/(X(i+1)-X(i));    %Forwards difference derivative of Y w.r.t X
    XPhalf(i)=(X(i+1)+X(i))/2.;             %X value at the midpoint
end

%INTERPOLATE
x=xvect(iLe:iTe);
for i=1:length(x)
    for k=1:length(XPhalf)-1
        if x(i)<XPhalf(k) && x(i)>XPhalf(k+1)
            break
        end
    end
    dydxtop(i)=(dYdX(k)+dYdX(k+1))/2;       %Slope on top
    for k=length(XPhalf):-1:2
        if x(i)<XPhalf(k) && x(i)>XPhalf(k-1)
            break
        end
    end
    dydxbot(i)=(dYdX(k)+dYdX(k+1))/2;       %Slope on bottom

     for k=1:length(X)-1
        if x(i)<X(k) && x(i)>X(k+1)
            break
        end
    end
    ytop(i)=(Y(k)+Y(k+1))/2;       %Coordinates on top
    for k=length(X):-1:2
        if x(i)<X(k) && x(i)>X(k-1)
            break
        end
    end
    ybot(i)=(Y(k)+Y(k+1))/2;       %Coordinates on bottom
end

% PLOT FOR DEBUGGING
% figure(99),close
% figure(99),hold on,axis([0 1 -1 1])
% plot(X,Y,'k-','linewidth',2)                %Airfoil shape
% plot(XPhalf,dYdX,'r-')                      %Slope from FD approx
% plot(x,dydxtop,'bo',x,dydxbot,'bo')         %Slope from Interpolation
% legend('Airfoil Shape','Finite Diff. Slope','Interpolation Slope')

dydxtop=[zeros(1,iLe-1) dydxtop];
dydxbot=[zeros(1,iLe-1) dydxbot];
ytop=[zeros(1,iLe-1) ytop];
ybot=[zeros(1,iLe-1) ybot];
