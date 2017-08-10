function [params1,stdResids,jj] = lineFit(xx,yy,REPEATS)

if nargin<3
    REPEATS = 100;
end
x1 = xx;
y1 = yy;

for jj = 1:REPEATS
    
    %First perform a fit
    params1 = polyfit(x1,y1,1);
    numPoints = length(xx);
    
    %Calculate the new coordinates based on the polynomial fit
    x2 = x1;
    y2 = polyval(params1,x2);

    %Calculate the residuals
    resids1 = y1-y2;
    
    %Discard points based on 3*sigma
    idx1 = abs(resids1)<3*std(resids1);
    x1 = x1(idx1);
    y1 = y1(idx1);

    
    if (length(x1) == length(x2))
        break;
    end
end

stdResids = std(resids1);


