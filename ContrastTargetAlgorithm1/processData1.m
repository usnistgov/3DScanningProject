close all; clear all;
%This script processes PTX files of contrast target scans (from a TLS) and
%calculates its center. The standard deviation values of these centers are
%calculated and compared against the output of OEM software. 

%Some setup
a1 = tic(); %To track time
addpath('.\code0730'); %This is the source for the contrast target code
addpath('.\code0730\supportFiles'); %These have some generic functions like planeFit3 etc. 
folder1 = '.\DATA'; %Source of all XYZI data and approx. centers

%Folders under the main data folder
folderCNT = fullfile(folder1,'Centers'); %Centers
folderPTX = fullfile(folder1,'PTX'); %XYZI files
folderPNG = fullfile(folder1,'PNG'); %For output of PNG files
folderMSK = folderCNT; %For the files needed for 2D mask

filesPTX = dir(fullfile(folderPTX,'*.ptx'));
filesCNT = dir(fullfile(folderCNT,'C*.txt'));
filesMSK = dir(fullfile(folderMSK,'Mask1.txt'));
fileMSK = filesMSK(1).name;
fname3 = fullfile(folderMSK,fileMSK);
maskData1 = dlmread(fname3);

%First read all the OEM s/w calculated centers and take their statistics
for kk = 1:length(filesCNT)  
    fileCNT = filesCNT(kk).name;
    fname2 = fullfile(folderCNT,fileCNT);    
    centerDataAll(kk,:,:,:) = dlmread(fname2);
end

for kk = 1:size(centerDataAll,2)
    centerData2 = squeeze(centerDataAll(:,kk,:));
    [az2,el2,rr2]  = cart2sph(centerData2(:,1), centerData2(:,2), centerData2(:,3));
    OEMCenterAVG(kk,:) = mean(centerData2);
    OEMCenterSTD(kk,:) = std(centerData2);
    
    OEMCenterAVGR(kk,:) = mean([az2 el2 rr2]);
    OEMCenterSTDR(kk,:) = std([az2 el2 rr2]);   
end




%For each scan, process the centers using the NIST developed code
for jj = 1:length(filesPTX)
    filePTX = filesPTX(jj).name;
    fname1 = fullfile(folderPTX,filePTX);
    pngFileName1 = fullfile(folderPNG, strcat('PNG1_',filePTX(1:end-4),'.png'));
    pngFileName2 = fullfile(folderPNG, strcat('PNG1_',filePTX(1:end-4),'.png'));
    
    disp('Loading PTX file...');
    tic; [I1,xyzData,intData] = ptx2img3(fname1);toc;
    figure(1); clf; imshow(I1);tic;
    set(gcf,'Position',[12   579   435   383]);
    
    %Obtain the 2D masks and for each 2D center process the PTX file
    for kk = 1:size(maskData1,1)
        maskParams = maskData1(kk,:);
        figure(1); hold on; plot(maskParams(1), maskParams(2),'ro'); 
        text(maskParams(1)+10,maskParams(2),num2str(kk));
        title(sprintf('File %d/%d: Target %d/%d',jj,length(filesPTX),kk,size(maskData1,1)));
        
        pause(1);
        %maskParams(3) = floor(maskParams(3));
        
        %Obtain the results (center of the contrast target)
        results1 = contrastFun9(I1,xyzData,intData,1,maskParams);
        pCent4P = results1.pCent4P;
        
        
        xyzData2 = results1.xyzData2;
        intData2 = results1.intData2;
        
        td = [];
        [taz,tel,trr] = cart2sph(xyzData2(:,1), xyzData2(:,2), xyzData2(:,3));
        [td(:,1),td(:,2),td(:,3)] = sph2cart(taz-mean(taz),tel,trr);
        
        planeRes = planeFit3R(td); 
        %         dist1 = results1.dist1;
        %         allProjDist(jj) = dist1;
        finalCenter1(jj,kk,:) = results1.pCent1P;
        finalCenter2(jj,kk,:) = results1.pCent2P;
        finalCenter3(jj,kk,:) = results1.pCent3P;
        finalCenter4(jj,kk,:) = results1.pCent4P;
        centerFrom2DMethod(jj,kk,:) = results1.cent2D_3D;
        planeAngles(jj,kk,1) = planeRes.yz_angle;
        planeAngles(jj,kk,2) = planeRes.zx_angle;
        planeAngles(jj,kk,3) = planeRes.xy_angle;
        figure(104); print(gcf, '-dpng', pngFileName1 )
        figure(106); print(gcf, '-dpng', pngFileName2 )

    end
    disp('**********************************')
    toc; 
    disp('**********************************')
    pause(2);
    
end

%Calculate the statistics to perform a comparision between the NIST method
%and the OEM method.
for kk = 1:size(finalCenter4,2)
    
    fCenterData3 = squeeze(finalCenter4(:,kk,:));
    fCenterData2 = squeeze(centerFrom2DMethod(:,kk,:));

    [faz3,fel3,frr3]  = cart2sph(fCenterData3(:,1), fCenterData3(:,2), fCenterData3(:,3));
    [faz2,fel2,frr2]  = cart2sph(fCenterData2(:,1), fCenterData2(:,2), fCenterData2(:,3));

    centerAVG3(kk,:) = mean(fCenterData3)
    centerSTD3(kk,:) = std(fCenterData3)
    
    centerAVG3R(kk,:) = mean([faz3 fel3 frr3]);
    centerSTD3R(kk,:) = std([faz3 fel3 frr3]);

    
    centerAVG2(kk,:) = mean(fCenterData2)
    centerSTD2(kk,:) = std(fCenterData2)
    
    centerAVG2R(kk,:) = mean([faz2 fel2 frr2]);
    centerSTD2R(kk,:) = std([faz2 fel2 frr2]);

    
    
end

%Comparing the distances
for kk = 1:size(centerDataAll,1)
c1 = squeeze(centerDataAll(kk,:,:));
d1 = rssq2(bsxfun(@minus,c1,c1(1,:)),2);
d1All(kk,:) = d1';

c2 = squeeze(finalCenter4(kk,:,:));
d2 = rssq2(bsxfun(@minus,c2,c2(1,:)),2);
d2All(kk,:) = d2';

end

RR_OEM = [OEMCenterAVGR(:,3) OEMCenterAVGR(:,3) OEMCenterAVGR(:,3)*0+1];
RR_NEW = [centerAVG3R(:,3) centerAVG3R(:,3) centerAVG3R(:,3)*0+1];

STD_NEW = centerSTD3R.*RR_NEW
STD_OEM = OEMCenterSTDR.*RR_NEW

ratio = centerSTD3./OEMCenterSTD

%Here were are comparing the standard deviations of the centers calculated
%by both the methods. 
ratioSphINST1 = centerSTD3R./OEMCenterSTDR
sum(ratioSphINST1(:))

%Here we are displaying the standard deviations of the centers in spherical
%coordinate system.
disp('STD - NIST method')
round(centerSTD3R*1000,1)

disp('STD - INST1')
round(OEMCenterSTDR*1000,1)

