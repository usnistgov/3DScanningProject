function [imageData,xyzData,intData,ww,hh] = ptx2img(ptx_file_name)
% Input is the .ptx file name
% Output is the imageData in intensity values from 0.0 to 1.0
data1 = dlmread(ptx_file_name);
ww = data1(1,1); %Take the width of the scan from 1st row
hh = data1(2,1); %Take the height of the scan from the 2nd row
intData = data1(11:end,4); %Read the rest of the data from the 11th row
xyzData = data1(11:end,1:3);

maxInt = max(intData(:));
minInt = min(intData(:));

dataImg2 = uint8(255*intData/maxInt);
imageData = (reshape(intData,hh,ww));
% 
% fname2 = strcat(ptx_file_name,'.mat');
% save(fname2,'imageData','xyzData','intData')
