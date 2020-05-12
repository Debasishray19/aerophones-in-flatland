clc; 
close all;
clear;

% STEP1: Extract tongue contour and palate contour coordinates from the file data
fileTongueID1 = fopen('vowel_a.txt','r');
fileTongueID2 = fopen('vowel_a.txt','r');

fileUpperPalateID1 = fopen('upperPalate.txt','r');
fileUpperPalateID2 = fopen('upperPalate.txt','r');

formatSpec_x = '%f %*s %f %*s %f %*s %f %*s %f %*s %*[^\n]';
formatSpec_y = '%*s %f %*s %f %*s %f %*s %f %*s %f %*[^\n]';

formatXSpec_palate = '%f %*s %*[^\n]';
formatYSpec_palate = '%*s %f %*[^\n]';

tongueContourSamplesCell_x = textscan(fileTongueID1, formatSpec_x);
tongueContourSamplesCell_y = textscan(fileTongueID2, formatSpec_y);
tongueContourSamplesMat_x = cell2mat(tongueContourSamplesCell_x);
tongueContourSamplesMat_y = cell2mat(tongueContourSamplesCell_y);

upperPalateSampleCell_x = textscan(fileUpperPalateID1, formatXSpec_palate);
upperPalateSampleCell_y = textscan(fileUpperPalateID2, formatYSpec_palate);

upperPalateSampleMat_x = cell2mat(upperPalateSampleCell_x);
upperPalateSampleMat_y = cell2mat(upperPalateSampleCell_y);

% STEP2: Interpolate tongue contour at exactly same coordinate the palate has.
sampleNum = 5; % Out of 5 samples, which one you are considering
tongueContourApproxSampleMat_x = upperPalateSampleMat_x;
tongueContourApproxSampleMat_y = ...
    interp1(tongueContourSamplesMat_x(:,sampleNum), tongueContourSamplesMat_y(:,sampleNum), tongueContourApproxSampleMat_x, 'spline', 'extrap');

% STEP3: Flip the y-coordinate values for tongue and upper palate
tongueContourApproxSampleMat_y = -1.*tongueContourApproxSampleMat_y;
upperPalateSampleMat_y = -1.*upperPalateSampleMat_y;

% STEP4: Find the center line between palate and the tongue
centerLine_y = (tongueContourApproxSampleMat_y+upperPalateSampleMat_y)./2;
centerLine_x = upperPalateSampleMat_x;

% STEP5: Find the perpendiculars lines' equation to the centre line
% centerLine_X length
numPoints = length(centerLine_x-1); 

% Points at which we need to find the centre line
perpendicularPoint_x0 =...
    (centerLine_x(1:numPoints-1)+centerLine_x(2:numPoints))./2;
perpendicularPoint_y0 =...
    (centerLine_y(1:numPoints-1)+centerLine_y(2:numPoints))./2;

slopeCentreLine = (centerLine_y(2:numPoints) - centerLine_y(1:numPoints-1))./...
    (centerLine_x(2:numPoints) - centerLine_x(1:numPoints-1));

slopePerpendicularLine = -1./slopeCentreLine; 

% Open figure window and plot the centre line
figure; hold on;
plot(tongueContourApproxSampleMat_x, tongueContourApproxSampleMat_y);
plot(upperPalateSampleMat_x, upperPalateSampleMat_y);
plot(centerLine_x, centerLine_y);

% STEP6: Draw the line perpendicular to the centre Line
yMax = max(upperPalateSampleMat_y);
yMin = min(tongueContourApproxSampleMat_y);

yMax = yMax+10;
yMin = yMin+20;

% Create a matrix to store the perpendicular line coordinates
lowerBound=5;
upperBoundSubtractor = 20;
stepSize = 10;

rowSize = 100;
columnsize = length(lowerBound:stepSize:numPoints-upperBoundSubtractor);

perpendicularLineYCoordinates = zeros(rowSize, columnsize);
perpendicularLineXCoordinates = zeros(rowSize, columnsize);
matrixCounter=1;

for counter = lowerBound:stepSize:numPoints-upperBoundSubtractor
    perpendicularToCentreLine_y = linspace(yMin, yMax, 100);
    
    y_y0 = perpendicularToCentreLine_y - perpendicularPoint_y0(counter);
    
    perpendicularToCentreLine_x = ((y_y0)./slopePerpendicularLine(counter))...
        +perpendicularPoint_x0(counter);
    
    perpendicularLineXCoordinates(:,matrixCounter) = perpendicularToCentreLine_x;
    perpendicularLineYCoordinates(:,matrixCounter) = perpendicularToCentreLine_y;
     
    matrixCounter=matrixCounter+1;
    plot(perpendicularToCentreLine_x, perpendicularToCentreLine_y,'-k');
end

axis equal;
hold off;

% Close open files
fclose(fileTongueID1);
fclose(fileTongueID2);
fclose(fileUpperPalateID1);
fclose(fileUpperPalateID2);

% Find the upper palate and tongue coordinates through which
% perpendicularLine passes

% For the upper palate
upperPalataLeftCoordinates = zeros(2, columnsize);
upperPalataRightCoordinates = zeros(2, columnsize);
perpendicularLineNum = 1;

for counter = 1: length(upperPalateSampleMat_x)-1
    a1 = (upperPalateSampleMat_x(counter) - perpendicularLineXCoordinates(1, perpendicularLineNum))*...
        (perpendicularLineYCoordinates (2, perpendicularLineNum) - perpendicularLineYCoordinates (1, perpendicularLineNum));
    
    b1 = (upperPalateSampleMat_y(counter) - perpendicularLineYCoordinates(1, perpendicularLineNum))*...
        (perpendicularLineXCoordinates (2, perpendicularLineNum) - perpendicularLineXCoordinates (1, perpendicularLineNum));
    
    d1 = a1-b1; 
     
    a2 = (upperPalateSampleMat_x(counter+1) - perpendicularLineXCoordinates(1, perpendicularLineNum))*...
        (perpendicularLineYCoordinates (2, perpendicularLineNum) - perpendicularLineYCoordinates (1, perpendicularLineNum));
    
    b2 = (upperPalateSampleMat_y(counter+1) - perpendicularLineYCoordinates(1, perpendicularLineNum))*...
        (perpendicularLineXCoordinates (2, perpendicularLineNum) - perpendicularLineXCoordinates (1, perpendicularLineNum));
    
    d2 = a2-b2;
    
    if d1*d2 < 0
        upperPalataLeftCoordinates(:,perpendicularLineNum) = [upperPalateSampleMat_x(counter), upperPalateSampleMat_y(counter)];
        upperPalataRightCoordinates(:,perpendicularLineNum) = [upperPalateSampleMat_x(counter+1), upperPalateSampleMat_y(counter+1)];
        perpendicularLineNum = perpendicularLineNum+1;
        if perpendicularLineNum > columnsize
            break;
        end
    end
end

% For the tongue palate
tongueLeftCoordinates = zeros(2, columnsize);
tongueRightCoordinates = zeros(2, columnsize);
perpendicularLineNum = 1;

for counter = 1: length(upperPalateSampleMat_x)-1
    a1 = (tongueContourApproxSampleMat_x(counter) - perpendicularLineXCoordinates(1, perpendicularLineNum))*...
        (perpendicularLineYCoordinates (2, perpendicularLineNum) - perpendicularLineYCoordinates (1, perpendicularLineNum));
    
    b1 = (tongueContourApproxSampleMat_y(counter) - perpendicularLineYCoordinates(1, perpendicularLineNum))*...
        (perpendicularLineXCoordinates (2, perpendicularLineNum) - perpendicularLineXCoordinates (1, perpendicularLineNum));
    
    d1 = a1-b1; 
     
    a2 = (tongueContourApproxSampleMat_x(counter+1) - perpendicularLineXCoordinates(1, perpendicularLineNum))*...
        (perpendicularLineYCoordinates (2, perpendicularLineNum) - perpendicularLineYCoordinates (1, perpendicularLineNum));
    
    b2 = (tongueContourApproxSampleMat_y(counter+1) - perpendicularLineYCoordinates(1, perpendicularLineNum))*...
        (perpendicularLineXCoordinates (2, perpendicularLineNum) - perpendicularLineXCoordinates (1, perpendicularLineNum));
    
    d2 = a2-b2;
    
    if d1*d2 < 0
        tongueLeftCoordinates(:,perpendicularLineNum) = [tongueContourApproxSampleMat_x(counter), tongueContourApproxSampleMat_y(counter)];
        tongueRightCoordinates(:,perpendicularLineNum) = [tongueContourApproxSampleMat_x(counter+1), tongueContourApproxSampleMat_y(counter+1)];
        perpendicularLineNum = perpendicularLineNum+1;
        if perpendicularLineNum > columnsize
            break;
        end
    end
end

perperndicularFirstTwoCoordinates_x = perpendicularLineXCoordinates(1:2,:);
perperndicularFirstTwoCoordinates_y = perpendicularLineYCoordinates(1:2,:);

% l1: Palate/Tongue l2: Perpendicular Lines
lineSlope_ml2 = (perperndicularFirstTwoCoordinates_y(2,:) - perperndicularFirstTwoCoordinates_y(1,:))./...
            (perperndicularFirstTwoCoordinates_x(2,:) - perperndicularFirstTwoCoordinates_x(1,:));
        
lineSlopeUpperPalalte_ml1 = (upperPalataRightCoordinates(2,:) - upperPalataLeftCoordinates(2,:))./...
            (upperPalataRightCoordinates(1,:) - upperPalataLeftCoordinates(1,:));

lineSlopeTongue_ml1 = (tongueRightCoordinates(2,:) - tongueLeftCoordinates(2,:))./...
            (tongueRightCoordinates(1,:) - tongueLeftCoordinates(1,:));
        
% Intercept coordinates for the upper palate
den1 = upperPalataLeftCoordinates(2,:) - perperndicularFirstTwoCoordinates_y(1,:);

den2 = lineSlope_ml2.*perperndicularFirstTwoCoordinates_x(1,:)-...
       lineSlopeUpperPalalte_ml1.*upperPalataLeftCoordinates(1,:);
   
num1 = lineSlope_ml2 - lineSlopeUpperPalalte_ml1;

x_incepUP = (den1+den2)./num1;
y_incepUP = upperPalataLeftCoordinates(2,:) + lineSlopeUpperPalalte_ml1.*...
            (x_incepUP - upperPalataLeftCoordinates(1,:));
        
% Intercept coordinates for the tongue
den1 = tongueLeftCoordinates(2,:) - perperndicularFirstTwoCoordinates_y(1,:);

den2 = lineSlope_ml2.*perperndicularFirstTwoCoordinates_x(1,:)-...
       lineSlopeTongue_ml1.*tongueLeftCoordinates(1,:);
   
num1 = lineSlope_ml2 - lineSlopeTongue_ml1;

x_incepTongue = (den1+den2)./num1;
y_incepTongue = tongueLeftCoordinates(2,:) + lineSlopeTongue_ml1.*...
            (x_incepTongue - tongueLeftCoordinates(1,:));
        
% Find Diameter
diameter = sqrt((y_incepTongue-y_incepUP).^2 + (x_incepTongue-x_incepUP).^2);
area_function = pi.*(diameter/2).^2;
