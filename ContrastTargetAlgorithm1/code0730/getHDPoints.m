function [data2] = getHDPoints(data1,scaleFactor)
%This function scales the points based on cubic polynomial interpolation

%Extract the 3D data into 3 variables
xx1 = data1(:,1);  yy1 = data1(:,2); zz1 = data1(:,3);

%Get the number of points based on the scale factor
numPoints = sqrt(length(xx1)*scaleFactor);

%Get the evenly spaced grid in X and Y
xx2 = min(xx1):range(xx1)/numPoints:max(xx1);
yy2 = min(yy1):range(yy1)/numPoints:max(yy1);
[xx2,yy2] = meshgrid(xx2,yy2);

%Use griddata to obtain interpolated points
zz2 = griddata(xx1,yy1,zz1,xx2,yy2,'cubic');

%Combine all the 3 axes into a single one and return.
data2 = double([xx2(:) yy2(:) zz2(:)]);

idxN = any(isnan(data2),2);

data2 = data2(~idxN,:);



% figure(1); hold on;
% plot3D([xx1,yy1,zz1])
% plot3D([xx2,yy2,zz2],'r')