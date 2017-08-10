function results = findImgCenter9(I1,xyzData,intData,FLAG, maskParameters)
% SUCCESS results.FOUND = 1
% FAILURE results.FOUND = 0
% NOT POSSIBLE TO FIND SOLUTION results.FOUND = -1;

UNIT_MULT = 1000; % xyzData is in mm
%Early on there was an issue with some images with the edges. The code
%below truncates the edges
chopPixels = 0;
I = imcrop(I1,[chopPixels+1 chopPixels+1 size(I1,2)-chopPixels size(I1,1)-chopPixels]);

%Get the rows and columns in the original image.
imgRows = size(I1,1);
imgCols = size(I1,2);


%Flags for plotting images
FLAG = 0; FLAG2 =0;
if(FLAG) figure(1); imshow(I), end

%% Here, generate a mask and apply it to the image

%Get the input center and radius for the mask
xc1 = maskParameters(1); yc1 = maskParameters(2); radius1 = maskParameters(3)*.8;

%Scale the pixesl radius to search
iRowIdx1 = floor(yc1)+1;
iColIdx1 = round(xc1);

%Get the index of the center based on the 2D center
tempCenterIdx = sub2ind(size(I1),floor(yc1),round(xc1));

%We expect about 35 pixels @ 10 m. Scale the radius value in the input
%parameters (maskparameters) according to the distance. More the distance
%less the radius. 
try
    tempPoint1 = xyzData(tempCenterIdx-10:tempCenterIdx+10,:);
catch ME  %It is likely that the 2D center doesn't belong to this scan
    disp(ME.identifier)
    disp(ME.causeException)
    error(ME.message);
end
tempDist1  = rssq2(tempPoint1,2);
tempDist1 = tempDist1(tempDist1>0);
tempPoint2 = diff(tempPoint1,2);
rr3  = rssq2(tempPoint2,2);

%Scale the radius value based on distance; 35 pixels for 10 m distance. New radius needs to be calculated for a different distance
%radius1 = radius1*10*UNIT_MULT/median(tempDist1); 


if isnan(radius1) & sum(tempDist1(:))==0
    warning('The XYZ data corresponding to the selected point are zeros');
    results.FOUND = -1;
    return;
elseif isnan(radius1) & sum(tempDist1(:))>0
    error('It is likely that the target was manually selected wrong');
end
results.newPixelRadius = radius1;

%Perform the masking operation
[rNum, cNum, ~] = size(I1);
[xx, yy] = ndgrid((1:rNum)-yc1, (1:cNum)-xc1); % generate grid with binary mask
mask1 = ( xx.^2 + yy.^2)>radius1^2; % define area outside of cirlce that will be masked
mask2 = (xx.^2 + yy.^2)<(0.05*radius1)^2;
I = I1;
I(mask1) = 0;% color is gray
%I(mask2) = 0;




%%%% The code below is modified from https://en.wikipedia.org/wiki/Chessboard_detection#Corners
USE_DILATION = 0;
if (USE_DILATION == 1)
    % Preprocess image to remove extraneous stuff
    BWs = edge(I,'canny');
    
    %Dilate image
    se90 = strel('line', 3, 90);
    se0 = strel('line', 3, 0);
    BWsdil = imdilate(BWs, [se90 se0]);
    
    
    %Fill interior gaps
    BWdfill = imfill(BWsdil, 'holes');
    
    %Remove Connected Objects on Border
    BWnobord = imclearborder(BWdfill, 4);
    
    
    seD = strel('diamond',1);
    BWFinal = BWnobord;
    for jj = 1:20
        BWFinal = imerode(BWnobord,seD);
    end
    I = BWFinal;
end


% Compute edge image
BW = edge(I,'canny');

% Compute Hough transform
[H theta rho] = hough(BW);

% Find local maxima of Hough transform
numpeaks = 2; %This was 19 in the original code
thresh = ceil(0.05 * max(H(:)));
P = houghpeaks(H,numpeaks,'threshold',thresh);

% Extract image lines
lines = houghlines(BW,theta,rho,P,'FillGap',50,'MinLength',60);

FLAG2 =1;
% If there are no lines, issue an error
if length(lines) >0
    
    
    %--------------------------------------------------------------------------
    % Display results
    %--------------------------------------------------------------------------
    
        if (FLAG2)
            % Original image
            figure(1); imshow(I);
            xlim([xc1-radius1*2 xc1+radius1*2]);
            ylim([yc1-radius1*2 yc1+radius1*2]);
            
            figure(111); clf; subplot(1,4,1);imshow(I1);axis equal;
            figure(111); subplot(1,4,2);imshow(I);axis equal;
            figure(111); subplot(1,4,3);imshow(I);axis equal;
            figure(111); subplot(1,4,4);imshow(I);axis equal;
            xlim([xc1-radius1*2 xc1+radius1*2]);
            ylim([yc1-radius1*2 yc1+radius1*2]);

            
        end
    
    
    
    % Obtain the individual lines
    n = size(I,2);
    color1 = 'gc';
    marker1 = 'do';
    for k = 1:length(lines)
        k2 = mod(k,2)+1;
        % Overlay kth line
        x = [lines(k).point1(1) lines(k).point2(1)] + chopPixels;
        y = [lines(k).point1(2) lines(k).point2(2)] + chopPixels;
        xy = [lines(k).point1; lines(k).point2] + chopPixels;
        params1(k,:) = polyfit(xy(:,1),xy(:,2),1);
        %save(strcat('xy',num2str(k),'.mat'),'xy');
        allLines((k-1)*2+1,:) = lines(k).point1;
        allLines((k-1)*2+2,:) = lines(k).point2;
        if (FLAG2)
            figure(111); subplot(1,4,3); hold on; line(xy(:,1),xy(:,2),'Color',color1(1)); axis equal;
            figure(111); subplot(1,4,4); hold on; line(xy(:,1),xy(:,2),'Color',color1(1)); axis equal;
            figure(1); hold on; line(xy(:,1),xy(:,2),'Color',color1(1));
        end
        

    end
    
    % If there is only one line, display it and then give an error
    if length(lines) < 2
        warning(sprintf('Found %d lines; Need 2 lines',length(lines)));
        results.FOUND = 0;
        return;
    end
    
    if length(lines)>2
        disp('More than 2 lines found');
        pause;
    end

    
    %Find intersection point
    A(:,1) = -params1(:,1);
    A(:,2) = 1;
    B = params1(:,2);
    
    if any(isinf(A(:))) | any(isinf(B(:)))
        A
        B
        disp('Line really vertical - Unlikely scenario');
        results.FOUND = 0;
        return;
    end
    X = A\B;
    cent2 = X; %Add the chopped pixels back
else
    cent2(1) = xc1; cent2(2) = yc1;
    warning(sprintf('Found %d lines; Need 2 lines',length(lines)));
end


if (FLAG) figure(1);plot(cent2(1), cent2(2),'bo'); end
if (FLAG2)
    figure(111); subplot(1,4,4); hold on; 
    plot(cent2(1),cent2(2),'bx');
    
    axis equal;
    %keyboard;
end




% Flip the 'y' back, beacuse the original image was flipped. This is the
% intersection point
centX = cent2(1);
centY = cent2(2); %If we have 6 rows, and 4th row is desired. When flipped, 4th row is not 6-4, but 6-4+1st row


%Now find the index of the 2D center
iRowIdx = floor(centY)+1;
iColIdx = round(centX);
centerIdx = (iColIdx-1)*imgRows + iRowIdx;

%Find the 3D center based on 2D center information.
cent3D = xyzData(centerIdx,:);
%disp(sprintf('%2.6f, ',cent3D));

%If the intensities from the image data and 4D data (XYZI) are not the same
%then there is an issue with the calculation.
if (I1(iRowIdx,iColIdx) ~= intData(centerIdx))
    warning('Intensities not equal');
end

cent2D = [centX,centY];

%Sometimes the radius found by the iterative method is too small, so we
%default to the user given radius for the mask (circle radius over the
%contrast target)
if radius1 < maskParameters(3)/2
    radius1 = maskParameters(3);
end

if results.newPixelRadius < maskParameters(3)/2
    results.newPixelRadius = maskParameters(3);
end


results.cent2D = cent2D;
results.radius1 = radius1;
results.cent3D = cent3D;
results.lines = lines;
results.centerR = 1; %This is not useful
results.maskParams = [cent2D radius1];
results.FOUND = 1;

