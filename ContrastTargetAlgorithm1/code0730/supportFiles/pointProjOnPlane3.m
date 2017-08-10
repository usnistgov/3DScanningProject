function [projPoint,normDist] = pointProjOnPlane(planeEq,point)
% close all; clear all;
N = planeEq; %Equation of the plane. Ax+By+Cz+D = 0; [A B C D]
P = point; %Point that needs to be projected on to the plane



norm1 = rssq2(N(1:3),2);
n = N(1:3)/norm1; %Normal vector
d = N(4)/norm1; %Constant

v = P; %Vector
normDist = dot(n,v) + d; %scalar distance of point to plane along the normal OR the distance of a point closest to the plane
%normDist can be positive or negative
projPoint = P - normDist*n; %Q is the projected point.

% 
% 
% t = (dot(P,n)+d)/dot(n,n)
% Q = P - n*t

FLAG = 0;
if (FLAG == 1)
    RNG = 2;
    X = [P(1)-RNG;P(1)+RNG];
    Y = [P(2)-RNG;P(2)+RNG];
    [X,Y] = meshgrid(X,Y);
    Z = (-N(1)*X -N(2)*Y-N(4))/N(3);
    reOrder = [1 2  4 3];
    figure(1);patch(X(reOrder),Y(reOrder),Z(reOrder),'b');
    grid on;
    alpha(0.3);
    
    
    figure(1); hold on;
    % plot3(X,Y,Z)
    % hold on;
    plot3(P(1),P(2),P(3),'rx')
    plot3(projPoint(1),projPoint(2),projPoint(3),'kx')
    text(projPoint(1),projPoint(2),projPoint(3),'Q')
    
    title(sprintf('%2.4f,',projPoint))
    
    LL = [P;projPoint];
    line(LL(:,1), LL(:,2), LL(:,3));
end

%% METHOD 1
% The closest point is along the normal to the plane. So define a point Q that is offset from P along that normal.
% 
% Q = P - n*t
% Then solve for t that puts Q in the plane:
% 
% dot(Q,n) + d = 0
% dot(P-n*t,n) + d = 0
% dot(P,n) - t*dot(n,n) = -d
% t = (dot(P,n)+d)/dot(n,n)

%% METHOD 2
% http://stackoverflow.com/questions/9605556/how-to-project-a-3d-point-to-a-3d-plane