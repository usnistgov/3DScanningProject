function center1 = findApproxIntersectionPt(data1)
% This function was derived from the following 
% https://blogs.mathworks.com/loren/2011/08/29/intersecting-lines/
%
%data1 is Nx3 data that has two intersecting planes in 3D
FLAG = 0;
data2 = data1(:,1:2);

%Perform a convex hull operation and get the outer edge points of the X
k1 = convhull(data2);
data2H  = data2(k1,:);
         
%Since these 4 sets of points should be uniquely at the 4 corners, they can
%be separated easily using kmeans.
idxH = kmeans(data2H,4);
if FLAG == 1
    color1 = 'rgbk';
    plot(data2(:,1),data2(:,2),'.'); hold on;
    plot(data2H(:,1), data2H(:,2),'ro')
    
    
    for jj = 1:4
        plot(data2H(idxH==jj,1), data2H(idxH==jj,2),'o','Color',color1(jj))
    end
end


%If there are clusters of data for each corner, then get the mean of those,
%so that only one point represents the outer edge of the line. 
for jj = 1:4
    point1(jj,:) = mean(data2H(idxH==jj,:),1);
end

% Perform intersection (instead of taking a mean)
center1 = findIntPoint3(point1);

center1 = [center1 mean(data1(:,3))];



function [finalIntPoint, order1] = findIntPoint3(pt)
finalIntPoint = [NaN NaN];
line1 = pt([1:2],:);
line2 = pt([3:4],:);
idx(1,:) = [1 2 3 4];
[inside(1),intPoint(1,:)] = intersectLines1(line1,line2);

line1 = pt([1,3],:);
line2 = pt([2,4],:);
idx(2,:) = [1 3 2 4];
[inside(2),intPoint(2,:)] = intersectLines1(line1,line2);

line1 = pt([1,4],:);
line2 = pt([2,3],:);
idx(3,:) = [1 4 2 3];
[inside(3),intPoint(3,:)] = intersectLines1(line1,line2);

for jj = 1:3
    if inside(jj)==1
        finalIntPoint = intPoint(jj,:);
        order1 = idx(jj,:);
        break;
    end
end


function [inside,intPoint] = intersectLines1(line1,line2)
FLAG = 0;
%% Example
% Let's start with with some simple line segments.  The first column
% represents to x-values of the line segment, the second column represents
% the y-values.
% line1 = [line1P1;line1P2]
% line2 = [line2P1;line2P2]
% 
% line1 = [.5 1; 1 0];
% line2 = [.5 0; 1 1];


%%
% Plot the data.  
if (FLAG)
h = plot(line1(:,1),line1(:,2));
hold all
h(2) = plot(line2(:,1),line2(:,2));
set(h,'linewidth',2)
axis([0 1 0 1])
end
 
%% First Algorithm
% Using equation y = mx+b, solve for x assuming 2 lines intersect.  Then
% see if that x value is in the necessary range.  Special cases: vertical
% lines (m==inf) and parallel lines (m1 == m2)
%%
% First find slopes and intercepts for both line segments.
% Here are slopes.
slope = @(line) (line(2,2) - line(1,2))/(line(2,1) - line(1,1));
m1 = slope(line1);
m2 = slope(line2);
%%
% Here are the intercepts.
intercept = @(line,m) line(1,2) - m*line(1,1);
b1 = intercept(line1,m1);
b2 = intercept(line2,m2);
xintersect = (b2-b1)/(m1-m2);
yintersect = m1*xintersect + b1;
%%
% Plot results.
if FLAG
plot(xintersect,yintersect,'m*','markersize',8)
legend({'line1','line2','intersection'},'Location','West')
hold off
end
%%  
% Now check that the intersection point is on the line segment and not past
% the end.  
isPointInside = @(xint,myline) ...
    (xint >= myline(1,1) && xint <= myline(2,1)) || ...
    (xint >= myline(2,1) && xint <= myline(1,1));
inside = isPointInside(xintersect,line1) && ...
         isPointInside(xintersect,line2);
intPoint = [xintersect,yintersect];     
   