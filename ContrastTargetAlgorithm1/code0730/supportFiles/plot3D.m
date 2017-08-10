function [hh] = plot3D(data1,color,mk1);

if nargin==1
    color = [0 1 1];
    mk1 = '.';
end
if nargin==2
    mk1 = '.';
end

hh = plot3(data1(:,1),data1(:,2),data1(:,3),'Marker',mk1,'Color',color,'LineStyle','none');
xlabel('X');ylabel('Y'); zlabel('Z');

