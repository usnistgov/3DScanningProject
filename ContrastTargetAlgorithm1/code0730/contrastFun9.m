function results = contrastFun9(I1,xyzData,intData,FLAG,maskParams)
%maskParams are the center and the radius of a mask that is applied to a
%contrast target so as to exclude the data outside this circular region.
%This helps the code to find the 3D center more reliably.
%
%I1,xyzData,intData are the output of ptximg3. I1 is the reshaped image
%date in mxn form. xyzData and intData are of the form Nx3.

% If the 3D results are not found, iterate by reducing the size of the
% mask. Sometimes other lines creep into the mask and finding the center
% becomes difficult
radius11 = maskParams(3);
for mm  = 1:30
    results = contrast8b_fun(I1,xyzData,intData,FLAG,maskParams);
    pause(0.1);
    if results.FOUND == 0 & maskParams(3)>= floor(0.3*radius11)
        maskParams(3) = maskParams(3)*.75;
        disp(sprintf('Decreasing the radius from %d to %2.1f',radius11,maskParams(3)));
    elseif results.FOUND == 1 
        return;
    else
        break;
    end
end

maskParams(3) = radius11;
frac = 1;
for mm  = 1:30
    results = contrast8b_fun(I1,xyzData,intData,FLAG,maskParams);
    pause(0.1);
    if results.FOUND == 0 & maskParams(3)<=2*radius11
        frac = 1.05^mm;
        disp(sprintf('Increasing the radius from %d to %2.1f',radius11,maskParams(3)));
        maskParams(3) = maskParams(3)*frac;
    else
        return;
    end
    if frac>1.5
        results.FOUND = 0;
        return;
    end
end




function results = contrast8b_fun(I1,xyzData,intData,FLAG,maskParams)
% close all; clear all;
% addpath('C:\Users\prem\Documents\MATLAB\160325-DataProcessingProject\processingCode\processSData2\utilities');
% ptxFile = 'Scan_2016_8_5_156-tt.ptx';


%% Constants
%This is the UL and LL fraction of STD of the intensities
%ulFrac = 1; llFrac = 1;
ulFrac = 1; llFrac = 1;

%inFrac = 0.1; outFrac = 0.1; %Fraction of the data to truncate to remove the inner and outer data of the circular region (inFrac was 0.2 when it worked)
%DATA_DENSITY = 500; %Increase in the data density based on griddata
inFrac = 0.3; outFrac = 0.1; %Fraction of the data to truncate to remove the inner and outer data of the circular region (inFrac was 0.2 when it worked)
DATA_DENSITY = 300; %Increase in the data density based on griddata

%% First rotate the data. The results can be rotated back later
% rotMat = pca(xyzData);
% xyzData = xyzData*rotMat; 

%% First find the approximate center in 2D and 3D from PTX file
% This uses hough transform and hough lines. This function finds the
% intersection of the 2 lines of interest and also the corresponding
% closest 3D point. No 3D center refinement is done here.



resultsT = findImgCenter9(I1,xyzData,intData,0,maskParams);
if resultsT.FOUND <= 0
    results.FOUND = 0;
else
    imgRows = size(I1,1);
    centX = maskParams(1);
    centY = maskParams(2);
    resultsT.maskParams = maskParams;
    iRowIdx = floor(centY)+1;
    iColIdx = round(centX);
    iCenterIdx = (iColIdx-1)*imgRows + iRowIdx;
    iCent3D = xyzData(iCenterIdx,:);
    resultsT.cent3D = iCent3D;
    
end


maskParams = resultsT.maskParams;

cent2D  = resultsT.maskParams(1:2);
cent3DApprox = resultsT.cent3D;
centerR = rssq2(cent3DApprox,2);
radius1 = resultsT.newPixelRadius;


%% Obtain the 3D center by considering the data in AzElIn domain
% Mask the region that is not needed
x1 = cent2D(1);
y1 = cent2D(2);
radius1 = resultsT.newPixelRadius; %Based on the distance, the pixel value will change. But, if scan PPD changes, this may not work
[rNum,cNum,~] = size(I1);

% Generate grid with binary mask representing the circle. This will select
% the region of interest
[xx,yy] = ndgrid((1:rNum)-y1,(1:cNum)-x1);
mask1 = (xx.^2 + yy.^2)>(0.1*radius1)^2 & (xx.^2 + yy.^2)<radius1^2;
%I2 = I1; I2(~mask1) = 0;
idx1 = mask1(:);
xyzData2 = xyzData(idx1,:); intData2 = intData(idx1,:); %Data without points at the origin

% Delete data that is at the center and also points at the origin.
[~,~,tr1] = cart2sph(xyzData(:,1), xyzData(:,2), xyzData(:,3));
idx2 = tr1>0 & idx1;
xyzData3= xyzData(idx2,:); intData3 = intData(idx2,:); %Data not at origin and not at center of the plate.
intData3O = intData3; %Just saving the actual intensity data;

%Scaling the intensity data to go from 0 to 1; 
aa1 = intData3-min(intData3); %Makes the min = 0; 
intData3 = aa1/max(aa1); %Scaling


if isempty(xyzData2) | isempty(xyzData3)
    warning('The masked 3D coordinates are empty. Possible reasons:\n 1. The manual locaiton selection may not be right\n2.The contrast target is curved or has artificial lines\n3.The radius of the mask paramters is too low');
    results.FOUND = 0
    return;
end

figure(106);clf; plot4D(xyzData3,intData3)
set(figure(106),'Position',[2245 15 630 350]);




% Extract data in the HVI domain.
[az3,el3,rr3] = cart2sph(xyzData3(:,1), xyzData3(:,2), xyzData3(:,3));

hviData3 = [az3 el3 intData3];



% Increase the density of the points
oldHVIData3 = hviData3; oldIntData3 = intData3; oldrr3 = rr3;
%keyboard;
[hviData3] = getHDPoints(hviData3,DATA_DENSITY);
intData3 = hviData3(:,3); %Since we are using some cubic extrapolation, intData could -ve values


hviData31 = hviData3;
intData31 = intData3;
% Get the 4 sets of data that belong to the region of transition, in the
% AzElIn domain. Based on a certain threshold, consider only the points
% between the topwhite part and the bottom black part.
tRR1 = rssq2(hviData31(:,1:2),2); tRR1 = tRR1-mean(tRR1);
meanIntVal = mean(oldIntData3);
idxW = oldIntData3>meanIntVal;
idxB = ~idxW;
intUL = meanIntVal + ulFrac*std(oldIntData3(idxW));
intLL = meanIntVal - llFrac*std(oldIntData3(idxB));


idxBW1 = intData31>intLL & intData31<intUL;
figure(104); clf; 
plot3D(oldHVIData3,[0.5 0.5 0.5]); hold on;
plot3D(hviData31(idxBW1,:),[.9 .9 .9],'.');
set(figure(104),'Position',[1700 550 525 425])


%Get the corresponding range values for the transition regions
idxBW2 = oldIntData3>intLL & oldIntData3<intUL;
rr4 = rr3(idxBW2);

hviData32 = hviData31(idxBW1,:);
intData32 = intData31(idxBW1,:);

meanVal = findApproxIntersectionPt(hviData32);
if any(isnan(meanVal)) %Could not find the intersection point of the 3D data in HVI/AzElIntensity domain
    results.FOUND = 0;
    return;
end

% Here the data is split into 4 parts by removing the cylindrical region in 
% the HVI domain. This is done, first by translating the origin to
% the center of the X region then using kmeans() to separate the data.
tData1 = bsxfun(@minus, hviData32(:,1:2), meanVal(1:2));
tRR1 = rssq2(tData1,2);

% A "cylindrical hole" is created by excluding certain data. This way the 
% X shaped data % is now visibly separated. Another step will separate data 
% into different variables. 
tIdx1 = tRR1 > min(tRR1) + inFrac*range(tRR1) & tRR1 < median(tRR1)+1*std(tRR1) - outFrac*range(tRR1);
hviData4 = hviData32(tIdx1,:);
intData4 = intData32(tIdx1,:);




%Use kmeans to separate the data
stdResidThreshold = 3E-4;
TRIES = 20;
for mj = 1:TRIES
    idxP = kmeans(hviData4(:,1:2),4); % Split based on only 2 dimensions
    for kk = 1:4
        lineData1 = hviData4(idxP==kk,1:2);
        [~,stdResids1(kk),repeats1(kk)] = lineFit(lineData1(:,1), lineData1(:,2));
        %         planeData1 = hviData4(idxP==kk,:);
        %         [~,stdResids1(kk),repeats1(kk)] = planeLineFit(planeData1,1)
    end
    %if (mod(mj+1,3)==1) disp(sprintf('%d: kmeansMetric = %2.3E (Should be < %2.3E)',mj, max(stdResids1), stdResidThreshold)); end
    
    %Go through line combinations of the lines to see if any lines intersect
    combo1 = [1 1 1 2 2 3];
    combo2 = [2 3 4 3 4 4];
    for mc = length(combo1)
        lineData1 = hviData4(idxP==combo1(mc),1:2);
        lineData2 = hviData4(idxP==combo2(mc),1:2);
        intersect1(mc) = sum(inhull(lineData1,lineData2));       
    end
    
    if sum(intersect1)>0
        disp('*********************')
        disp(sprintf('Intersections = %d',sum(intersect1)));
        disp('*********************')
    end
        %keyboard;
    if max(stdResids1) < stdResidThreshold & sum(intersect1) == 0 %There is no intersection
        break;
    end
end
if mj == TRIES
    disp('No center found');
    results.FOUND = 0;
    return;
end

disp('KMeans');
for mj = 1:4
    disp(sprintf('Line %d: STD = %2.6f, Repeats = %2.6f',mj,stdResids1(mj), repeats1(mj)));
end



color1 = 'rgbkm';
marker1 = 'sdov';
figure(104);
for kk = 1:4
    segment1 = hviData4(idxP==kk,:);
    hold on; plot3D(segment1,color1(kk),marker1(kk));
end
view(0,90)


%%Below are multiple methods to calculate the centers. Considering two
%%lines and intersecting them works better. Others are not necessary
% 1 = 5 planes. 5th Plane = plane in ranging direction
% 2 = 4 planes
% 3 = 4 lines
% 4 = 2 lines


%% Now fit planes to these 4 data sets and intersect them
pts1 = hviData4(idxP==1,:);
pts2 = hviData4(idxP==2,:);
pts3 = hviData4(idxP==3,:);
pts4 = hviData4(idxP==4,:);
pts5 = hviData3(intData3>meanIntVal-0.01&intData3<meanIntVal+0.01,:);



%Fit a plane
res1(1) = planeFit3R(pts1);
res1(2) = planeFit3R(pts2);
res1(3) = planeFit3R(pts3);
res1(4) = planeFit3R(pts4);
res1(5) = planeFit3R(pts5);


% Intersect the 4 planes and get the center in the HVI domain
for jj = 1:length(res1)
    pParams1(jj,:) = res1(jj).Parameters;
end
AP = pParams1(1:4,1:3);
BP = -pParams1(1:4,4);
XP = AP\BP; % Center in HVI domain
%centerR = mean(rr4);

% From the center in HVI, get centerH, centerV and get the 3D center
[pCent1(1),pCent1(2),pCent1(3)] = sph2cart(XP(1),XP(2),centerR);


AP = pParams1(:,1:3);
BP = -pParams1(:,4);
XP = AP\BP; % Center in HVI domain
[pCent2(1),pCent2(2),pCent2(3)] = sph2cart(XP(1),XP(2),centerR);


%% Use the 4 sets of data and intersect them as lines
pts1L = pts1; pts1L(:,3) = [];
pts2L = pts2; pts2L(:,3) = [];
pts3L = pts3; pts3L(:,3) = [];
pts4L = pts4; pts4L(:,3) = [];



%Fit a Line
tp1 = lineFit(pts1L(:,1),pts1L(:,2),1); lParams1(1,:) = [-tp1(1) 1 -tp1(2)]; clear tp1;
tp2 = lineFit(pts2L(:,1),pts2L(:,2),1); lParams1(2,:) = [-tp2(1) 1 -tp2(2)]; clear tp2;
tp3 = lineFit(pts3L(:,1),pts3L(:,2),1); lParams1(3,:) = [-tp3(1) 1 -tp3(2)]; clear tp3;
tp4 = lineFit(pts4L(:,1),pts4L(:,2),1); lParams1(4,:) = [-tp4(1) 1 -tp4(2)]; clear tp4;

AL = lParams1(:,1:2);
BL = -lParams1(:,3);
XL = AL\BL; % Center in HVI domain
%centerR = mean(rr4);


[pCent3(1),pCent3(2),pCent3(3)] = sph2cart(XL(1),XL(2),centerR);


%% Now, try to combine the 4 lines into two lines based on their slopes.

% Lines with equal or near-equal slopes should be clubbed together.
slopes1 = -lParams1(:,1);
[~,idxS1] = sort((slopes1));
%tan() is positive in the diagonal quadrants, so -ve values will be sorted
%first and the positive values last.
figure(104);
color2 = 'ccmm';
for kk = 1:4
    segment1 = hviData4(idxP==idxS1(kk),:);
    segment1C(kk,:) = mean(segment1,1);
    segment1C(kk,3) = 1;
    text(segment1C(kk,1),segment1C(kk,2),segment1C(kk,3),num2str(kk),'Color',color2(kk),'FontSize',18);
end
line1M = [mean(segment1C(1,1:2),1); mean(segment1C(2,1:2),1)];
line2M = [mean(segment1C(3,1:2),1); mean(segment1C(4,1:2),1)];

% [inside,intPoint] = isIntersecting2(line1M,line2M);
% if (inside == 0)
%     disp('The lines do not intersect');
%     results.FOUND = 0;
% end
% keyboard;


line1 = [hviData4(idxP==idxS1(1),1:2); hviData4(idxP==idxS1(2),1:2);];
line2 = [hviData4(idxP==idxS1(3),1:2); hviData4(idxP==idxS1(4),1:2);];

tp1 = lineFit(line1(:,1),line1(:,2),1); lParamsAB(1,:) = [-tp1(1) 1 -tp1(2)]; clear tp1;
tp2 = lineFit(line2(:,1),line2(:,2),1); lParamsAB(2,:) = [-tp2(1) 1 -tp2(2)]; clear tp2;
AL = lParamsAB(:,1:2);
BL = -lParamsAB(:,3);
XL = AL\BL; % Center in HVI domain
[pCent4(1),pCent4(2),pCent4(3)] = sph2cart(XL(1),XL(2),centerR);
hviCent = [XL' mean(hviData4(:,3))];
figure(104); plot3(hviCent(1), hviCent(2), hviCent(3),'ro','MarkerSize',10);

%keyboard;
%Performa a check to see if the center obtained by intersecting two lines
%in the HV domain is infact at the center  of these. If the logic is not
%right, the selected lines are not in the opposite quadrants always.
inHull = inhull(hviCent, hviData4);

if inHull == 0
    disp('The center is not the intersection of the HV lines');
    %keyboard;
    results.FOUND = 0;
    return;
end





% figure(106); hold on;
% %plot4D(xyzData2,intData2); hold on;
% plot3D(pCent1,'r','x');plot3D(pCent1,'r','o');
% plot3D(pCent2,'b','x');plot3D(pCent2,'b','o');
% plot3D(pCent3,'k','*');plot3D(pCent3,'k','d');
% plot3D(pCent4,'m','.');plot3D(pCent4,'m','s');

resP = planeFit3R(xyzData3,100);
planeEq = resP.Parameters;
[pCent1P] = pointProjOnPlane3(planeEq,pCent1); %Based on 4 planes
[pCent2P] = pointProjOnPlane3(planeEq,pCent2); %Based on 5 planes
[pCent3P] = pointProjOnPlane3(planeEq,pCent3); %Based on 4 lines
[pCent4P] = pointProjOnPlane3(planeEq,pCent4); %Based on 2 lines %We are using this for our 

figure(106); hold on;
plot4D(xyzData2,intData2); hold on;
% plot3D(pCent1P,'r','x');plot3D(pCent1,'r','o');
% plot3D(pCent2P,'b','x');plot3D(pCent2,'b','o');
% plot3D(pCent3P,'k','*');plot3D(pCent3,'k','d');
plot3D(pCent4P,'m','.');plot3D(pCent4P,'m','s');


results.pCent1 = pCent1;
results.pCent1P = pCent1P;
results.pCent2 = pCent2;
results.pCent2P = pCent2P;
results.pCent3 = pCent3;
results.pCent3P = pCent3P;
results.pCent4 = pCent4;
results.pCent4P = pCent4P;
results.finalCenter = pCent4P;
results.xyzData2 = xyzData2;
results.intData2 = intData2;
results.dist1 = 0;
results.FOUND = 1;


results.cent2D = resultsT.cent2D;
results.cent2D_3D = resultsT.cent3D;
results.lines = resultsT.lines;
%keyboard;