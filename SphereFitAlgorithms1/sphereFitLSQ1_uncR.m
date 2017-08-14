function result = sphereFitLSQ1(data,dummyVal1)
% This function takes 3 dimensional data of a sphere in a Cartesian 
% coordinate system and returns the paramters of a sphere (Center and
% radius)
%
%The input data is in the form of Nx3 matrix
%The output is a structure with the sphere parameters.
%The function also takes in an optional parameter "dummyVal1" which does
%not do anything, and is retained for consistency purposes with other
%similar algorithms

%If no "dummyVal" is given, just assign some random number. It is not used anyway. 
if nargin<2
    dummyVal = 1;
end

%Check the dimensions of the data
[m,n] = size(data);
if(n<3 | m<4)
    data = data';
end
[m,n] = size(data);

if(n<3 | m<4)
    size(data)
    error('Not enough parameters');
end

%First get an initial estimate based on a linearized form of the sphere equation
[ro,ao,bo,co,residualo] = sphereFitInitGuess(data);
initGuess = [ao bo co ro];

%Now formulate the objective function and the options for lsqnonlin() 
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxFunEvals',1500,'TolFun',1E-9,'display','off');
objectiveFunction = @(initGuess)calcResiduals(initGuess,data);
[result3] = lsqnonlin(objectiveFunction,initGuess,[],[],options);

%The above part works fine. An additional step can be performed by
%performing multiple initial estimates to (needs an additional toolbox).
%This is to ensure that the solution is not stuck at a local minimum.


% The "residuals" and "Residuals" are the same parameters, but retained for
% legacy purposes
result.Parameters = result3;
result.Residuals = calcResiduals(result3,data);
result.Center = result3(1:3)';
result.Radius = result3(4);



function residuals = calcResiduals(initGuess,data)
%Function to calculate the residuals based on the sphere parameters and the
%data.

ao = initGuess(1);
bo = initGuess(2);
co = initGuess(3);
ro = initGuess(4);

xx = data(:,1);
yy = data(:,2);
zz = data(:,3);
residuals = sqrt( (xx-ao).^2 + (yy-bo).^2 + (zz-co).^2 ) - ro;



function [r,a,b,c,residual] = sphereFitInitGuess(data)
%This is a function to calculate the parameters of a sphere as an initial
%estimate. This calculates both the center and the radius.
xx = data(:,1);
yy = data(:,2);
zz = data(:,3);

AA = [-2*xx, -2*yy , -2*zz , ones(size(xx))];
BB = [ -(xx.^2+yy.^2+zz.^2)];

YY = mldivide(AA,BB);

a = YY(1);
b = YY(2);
c = YY(3);
D = YY(4); % a^2 + b^2 + c^2 -r^2(where a,b,c are centers)
r = sqrt((a^2+b^2+c^2)-D);
residual = AA*YY-BB;


