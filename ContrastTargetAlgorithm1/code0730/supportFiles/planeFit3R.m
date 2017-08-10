function [result1,data2] = planeFit3R(data1,REPEATS)
if nargin<2
    REPEATS = 20;
end
for jj = 1:REPEATS
    result1 = planeFit3(data1);
    resids1 = result1.Residuals;
    idx1 = abs(resids1)<mean(resids1)+3*std(resids1);
    if(sum(idx1) == length(data1)) break; end;
    %figure(1); clf; plot3D(data1); pause(1)
    data1 = data1(idx1,:);
    %length(resids1)
end
data2 = data1;

