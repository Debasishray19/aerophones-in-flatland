% This code implements the following paper to illustrate the tube cross
% sectional area for the vowel sound - /a/, /u/, /i/ : 
% Comparision of Magnetic resonance Imaging-Based vocal tract area function
% obtained from the same speaker in 1994 and 2002.

% We will be using 44 cross sections to generate the vocal tract
% Important: Sectional/segment Length(delta) is different from spatial resolution (ds)

function [listenerX, listenerY, frameH, frameW, depthX, depthY, depthP, PV_N]...
         = crossSectionTube_SymmetricalGeometry(pmlSwitch, ds, pmlLayer, vowelSound, simulation2D, cell_wall, cell_air, cell_excitation, cell_noPressure)
     
    % Define units
    meter = 1;
    centimeter = 1e-2*meter;
    cellHalfLen = ds/2;
	
    % Vocal tract parameters
    numSections = 44;
    diameter_mul = 1;
    depth_mul = 1;
    
    % STEP0: Select the cross-sectional area based on the vowelSound 
    % Selct the cross sectional area based on the vowel sound
    % Tube section area in cm^2 in 3D
    switch vowelSound
        case 1
            sectional_length = 0.00388; % in meter
            tubeSectionArea_incm2_3D = [0.56 0.62 0.66 0.78 ...
                                        0.97 1.16 1.12 0.82 ...
                                        0.55 0.45 0.37 0.29 ...
                                        0.21 0.15 0.16 0.25 ...
                                        0.34 0.43 0.54 0.61 ...
                                        0.67 0.98 1.76 2.75 ...
                                        3.52 4.08 4.74 5.61 ...
                                        6.60 7.61 8.48 9.06 ...
                                        9.29 9.26 9.06 8.64 ...
                                        7.91 6.98 6.02 5.13 ...
                                        4.55 4.52 4.71 4.72];
        case 2
            sectional_length = 0.00445; % in meter
            tubeSectionArea_incm2_3D = [0.54 0.61 0.66 0.75 ...
                                        1.13 1.99 2.83 2.90 ...
                                        2.52 2.40 2.83 3.56 ...
                                        3.99 3.89 3.50 3.04 ...
                                        2.64 2.44 2.31 2.07 ...
                                        1.80 1.52 1.14 0.74 ...
                                        0.42 0.22 0.14 0.20 ...
                                        0.47 0.89 1.15 1.42 ...
                                        2.17 3.04 3.69 4.70 ...
                                        5.74 5.41 3.82 2.34 ...
                                        1.35 0.65 0.29 0.16];
            
        case 3
            sectional_length = 0.00384; % in meter
            tubeSectionArea_incm2_3D = [0.51 0.59 0.62 0.72 ...
                                        1.24 2.30 3.30 3.59 ...
                                        3.22 2.86 3.00 3.61 ...
                                        4.39 4.95 5.17 5.16 ...
                                        5.18 5.26 5.20 5.02 ...
                                        4.71 4.13 3.43 2.83 ...
                                        2.32 1.83 1.46 1.23 ...
                                        1.08 0.94 0.80 0.67 ...
                                        0.55 0.46 0.40 0.36 ...
                                        0.35 0.35 0.38 0.51 ...
                                        0.74 0.92 0.96 0.91];
            
        otherwise
            disp('Something Wrong - Check Your Code');
            return;
    end
                             
    % Tube section area in m^2
    tubeSectionArea_inm2_3D = tubeSectionArea_incm2_3D.*(centimeter*centimeter);    
    
    % STEP1: Calculate the tube section diameter with diameter multiplier
    tubeSectionDiameter_3D = 2*sqrt(tubeSectionArea_inm2_3D./pi).*diameter_mul;
    
    % STEP2: Convert the 3D cross sectional area into 2D using expansion
    % ratio. But we do not need this for 2.5D. Use a variable to switch
    % between 2.5D and 2D.
    
    % If simulation2D is 1 then we will be running 2D simulator otherwise 2.5D
    % will be running 
    if simulation2D==1 % Fpr 2D simulation
        % Tube section diameter in cm^2 in 2D
        % Create the array 
        tubeSectionDiameter_2D = zeros(1, numSections);
        
        % Find out the maxium diameter in 3D and it's index
        [maxDiameter, maxDiameterIdx] = max(tubeSectionDiameter_3D);
        
        % Assign the updated diameter to the tubeSectionDiameter_2D for the same
        % index position: d_2D = d3D(0.5*pi/1.84)   
        tubeSectionDiameter_2D(maxDiameterIdx) = (maxDiameter*0.5*pi)/1.84;
        
        % To Assign other diameters for tubeSectionDiameter_2D implement
        % eqn -(3a) and (3b) from the foloowing paper:
        % Two dimensional vocal tracts with three dimensional behaviour in the
        % numerical generation of vowels.

        % Store the value of m in an array - are expansion ratio
        % The variable 'm' has been named as per the parameters named in the
        % paper to avoid confusion
        m = ones(1, numSections);

        % Loop to define m value
        for counterSection = 1:numSections

            if counterSection < maxDiameterIdx
                m(counterSection) = (tubeSectionDiameter_3D(counterSection)/...
                                    tubeSectionDiameter_3D(counterSection+1))^2;
            end

            if counterSection > maxDiameterIdx
                m(counterSection) = (tubeSectionDiameter_3D(counterSection)/...
                                    tubeSectionDiameter_3D(counterSection-1))^2;
            end
        end

        % Loop to define 2D diameter
        radiiLessMaxIdx = maxDiameterIdx;

        while radiiLessMaxIdx > 1
            tubeSectionDiameter_2D(radiiLessMaxIdx-1) = ...
            tubeSectionDiameter_2D(radiiLessMaxIdx) * m(radiiLessMaxIdx-1);

            radiiLessMaxIdx = radiiLessMaxIdx - 1;
        end

        radiiGreaterMaxIdx = maxDiameterIdx;
        while radiiGreaterMaxIdx < numSections
            tubeSectionDiameter_2D(radiiGreaterMaxIdx+1) = ...
            tubeSectionDiameter_2D(radiiGreaterMaxIdx) * m(radiiGreaterMaxIdx+1);

            radiiGreaterMaxIdx = radiiGreaterMaxIdx + 1;
        end
        
        % Tube section diameter in terms of number of grid cell
        tubeSectionDiameterCells = round(tubeSectionDiameter_2D./ds);
        finalTubeSectionDiameter = tubeSectionDiameter_2D;
        
    else % For 2.5D simulation
        % Tube section diameter in terms of number of grid cell
        tubeSectionDiameterCells = round(tubeSectionDiameter_3D./ds); 
        finalTubeSectionDiameter = tubeSectionDiameter_3D;
    end
           
    % Change the tube diameter to 1 if it contains 0
    tubeSectionDiameterCells(tubeSectionDiameterCells==0)=1;
    
    %STEP3: Choose the best possible odd number from the Diameter array
    for diameterCounter = 1:numSections       
        % Verify if the cellsPerDiameter is odd or not
        if mod(tubeSectionDiameterCells(diameterCounter), 2) == 0
            
            % Find the difference between rounded and actual diameter value
            diff = tubeSectionDiameterCells(diameterCounter) - ...
                   (finalTubeSectionDiameter(diameterCounter)/ds); 
            
            if diff>0
                tubeSectionDiameterCells(diameterCounter) = ...
                    tubeSectionDiameterCells(diameterCounter)-1;
            else
                tubeSectionDiameterCells(diameterCounter)=...
                    tubeSectionDiameterCells(diameterCounter)+1;
            end
        end   
    end
  
    % STEP4: Find the total tube length and calculate the percentage error 
    % in the approximated tube length

    % Number of cells for total tube length
    actualTubeLength = numSections*sectional_length;
    totalTubeLengthinCells = round(actualTubeLength/ds);
    approxTubeLength = totalTubeLengthinCells*ds;
    
    % Percentage error in total tube length
    totalTubeLengthError = abs(actualTubeLength - approxTubeLength)/...
                           (actualTubeLength);
    fprintf('Error percentage in approximated tube length = %.4f \n',totalTubeLengthError*100);
    
    % STEP5: Construct the frame and PV_N array
    offsetW = 8;
    offsetH = 6;
    domainW = totalTubeLengthinCells + offsetW;
    domainH = max(tubeSectionDiameterCells) + offsetH;
    
    % Build frame [Frame = Domain Size + PML Layers]
    if pmlSwitch == 1
        frameW = domainW+2*pmlLayer+2; % 2 is for border/dead cells
        frameH = domainH+2*pmlLayer+2; % 2 is for border/dead cells
    else
        frameW = domainW+2;
        frameH = domainH+2;
    end
    
    % Define PV_N array
    % For pressure and velocity
    % PV_N(:,:,1) = To store cell pressure
    % PV_N(:,:,2) = To store Vx
    % PV_N(:,:,3) = To store Vy
    % PV_N(:,:,4) = To store grid cell type
    PV_N = zeros(frameH, frameW, 4);
    
    % Define cell type and store it in PV_N(,,4)
    % Declare all the cells as air by default
    PV_N(1:frameH, 1:frameW, 4) = cell_air;
    
    % STEP6: Define the depth matrix to store the tube height along z-axis
    if simulation2D == 1
        openSpaceDepth = 1;
    else
        openSpaceDepth = input('Enter the open space depth value: ');
    end
    
    minDepth  = 0.001;
    depthP    = ones(frameH, frameW)*openSpaceDepth*meter;
    depthX    = ones(frameH, frameW)*openSpaceDepth*meter;
    depthY    = ones(frameH, frameW)*openSpaceDepth*meter;
    
    
    % Array to store the tube upper and lower wall positions
    % This array will help us to major the depth/height across the tube
    % Need to calculate depth/height for each cell in a column
    % How many columns = totalTubeLengthinCells
    % ROW1=Upper Wall Idx, ROW2=Lower Wall Idx, ROW3=Radius
    cellDepthProp = zeros(2, totalTubeLengthinCells);
    
    % STEP7: Find the mid point in the frame
    midY = floor(frameH/2);
    midX = floor(frameW/2);
    
    % STEP8: Find the glottal end: Starting point of the tube
    tubeStartX = midX - round(totalTubeLengthinCells/2);
    tubeStartY = midY;
    tubeEndX   = tubeStartX + totalTubeLengthinCells-1;
    tubeEndY   = midY;
    
    % STEP9: Store the cummulative length of each tube section
    tubeCummSectionLength = zeros(1, numSections);
    for sectionCount = 1:numSections
        tubeCummSectionLength(sectionCount) = sectional_length*sectionCount;
    end
    
    % STEP10: Draw the geometrical shape of the tube
    % Set a counter to traverse through tubeCummSectionLength
    sectionCounter = 1;
    
    % To join the walls along the vertical direction if there is gap.
    prevUpperY = 0;
    prevUpperX = 0;
    prevLowerY = 0;
    prevLowerX = 0;
    
    for tubeCellsCount = 1:totalTubeLengthinCells
        
        % Check the current length
        currTubeLength = tubeCellsCount*ds;
        
        % Verify if the currTubeLength is more than the tubeCummSectionLength
        % for the current sectionCounter
        
        % if small or equal then set the tube wall as expected-Normal Case
        if currTubeLength <= tubeCummSectionLength(sectionCounter)
            
            % Get the tube Radius
            % We are subtracting 1 as we'll assume that there is a middle
            % row of cells which will act like a morror.
            getRadius = (tubeSectionDiameterCells(sectionCounter)-1)/2;
            
        else
        % If the currTubeLength is greater than the actual tubeCummSectionLength
            % Find the difference between currTubeLength and tubeCummSectionLength
            diffLength = currTubeLength - tubeCummSectionLength(sectionCounter);
            
            if diffLength>cellHalfLen && sectionCounter~=numSections
                % Increase the section counter
                sectionCounter = sectionCounter+1;
                
                % Get the radius for that cross-section
                getRadius = (tubeSectionDiameterCells(sectionCounter)-1)/2;                
            else
                % Get the radius for that cross-section
                getRadius = (tubeSectionDiameterCells(sectionCounter)-1)/2;
                
                % And then increase the section Counter
                sectionCounter = sectionCounter+1;
            end            
        end  % End of checking currTubeLength and tubeCummSectionLength
        
        % Find the upper and lower wall coordinates
        upperX = tubeStartX + (tubeCellsCount-1);
        upperY = tubeStartY - getRadius - 1;           
        lowerX = tubeStartX + (tubeCellsCount-1);
        lowerY = tubeStartY + getRadius + 1;
        
        % set the upper wall & lower wall
        PV_N(upperY, upperX,4) = cell_wall;
        PV_N(lowerY, lowerX,4) = cell_wall;
        
        % Store in the cellDepthProp to calculate cell depth
        cellDepthProp(1, tubeCellsCount) = upperY;
        cellDepthProp(2, tubeCellsCount) = lowerY;
            
        % Find the difference in prevUpperY and current upperY
        differYwall = upperY - prevUpperY;
        if tubeCellsCount~=1
            if differYwall>1
                % Set upper wall
                PV_N(prevUpperY:upperY, upperX,4) = cell_wall;
                
                % Set lower wall
                PV_N(lowerY:prevLowerY, lowerX,4) = cell_wall;
            elseif differYwall < -1
                % Set upper wall
                PV_N(upperY:prevUpperY, prevUpperX,4) = cell_wall;
                
                % Set lower wall
                PV_N(prevLowerY:lowerY, prevLowerX,4) = cell_wall;
            end
        end
        
        prevUpperY = upperY;
        prevUpperX = upperX;
        prevLowerY = lowerY;
        prevLowerX = lowerX;
    end % End of For loop
    
    % STEP11: Set the depthX, depthY and depthP if simulation2D~=1. 
    % For 2D simulation we do not need depth parameter.
    % Calculate depth/height of the tube along z-axis for depthX and depthY
    % For the boundary set the radius as zero
    if simulation2D~=1
        for tubeCellLengthCounter = 1:totalTubeLengthinCells
        
            % S1: Find the lowerY and upperY
            lowerY = cellDepthProp(2, tubeCellLengthCounter);
            upperY = cellDepthProp(1, tubeCellLengthCounter);

            % S2: Find the xPosition
            posX = tubeStartX + (tubeCellLengthCounter-1);

            % S3: Number of air cells between upper and lower wall cells:
            % Diameter of the current tube segment
            diameterCellsCount = lowerY-upperY-1;
            distanceCounter = ceil(diameterCellsCount/2)-1;

            % S4: Get the square of radius
            sqrRadius = ((diameterCellsCount/2)*ds)^2;

            % S5: Traverse through lower wall to center and set depth
            % as zero for the lower wall. Remember we are storing the depth for the
            % lower left corner of each grid cells
            cellCounter = distanceCounter;

            while cellCounter >= 0
                % Find the depth for the corresponding cell
                % For depthX
                sqrDistanceX = (cellCounter*ds)^2;
                depthSqrX = sqrRadius-sqrDistanceX;
                depthValX = sqrt(depthSqrX);

                % FordepthY
                sqrDistanceYlower = ((cellCounter-0.5)*ds).^2;
                sqrDistanceYupper = ((cellCounter+0.5)*ds).^2;
                depthSqrYlower    = sqrRadius - sqrDistanceYlower;
                depthSqrYupper    = sqrRadius - sqrDistanceYupper;
                
                depthValYlower    = sqrt(depthSqrYlower);
                depthValYupper    = sqrt(depthSqrYupper);
                %Check if depthSqrX or depthSqrY is negative or not
                if depthSqrX<0 || depthSqrYlower<0 || depthSqrYupper<0
                     disp('Error Height');
                     return;
                end

                % Find the upper and lower Y position for depthX and depthY
                posUpperYdepthX = upperY + (distanceCounter - cellCounter) + 1;
                posLowerYdepthX = lowerY - (distanceCounter - cellCounter) - 1;
                posUpperYdepthY = upperY + (distanceCounter - cellCounter) + 1;
                posLowerYdepthY = lowerY - (distanceCounter - cellCounter) - 1;

                % Set the depth for both upper and lower cell position
                depthX(posUpperYdepthX, posX) = 2*depthValX;
                depthX(posLowerYdepthX, posX) = 2*depthValX;
                depthY(posUpperYdepthY, posX) = 2*depthValYupper;
                depthY(posLowerYdepthY, posX) = 2*depthValYlower;
                
                % Decrease the cellCounter
                cellCounter = cellCounter-1;
            end

            depthX(lowerY, posX) = 0;
            depthY(lowerY, posX) = 0;
            depthX(upperY, posX) = 0;
            depthY(upperY, posX) = 0;
        end
        
        % Reset the openSpaceDepth to the central air cell width of the
        % last segment of the tube
        depthX(depthX==openSpaceDepth) = depthX(midY, tubeEndX);
        depthY(depthY==openSpaceDepth) = depthY(midY, tubeEndX);
                
        % To smooth out, average the current cell depthX value with the right
        % most cell.
        depthX(:, 1:frameW-1) = ( depthX(:, 1:frameW-1)+depthX(:, 2:frameW) )/2;
        depthX(:, frameW) = ( depthX(:, frameW)+depthX(midY, tubeEndX) )/2;
        
        depthY(1,:) = (depthY(1,:)+ depthY(midY, tubeEndX))/2;
        depthY(2:frameH,:) = (depthY(2:frameH,:)+ depthY(1:frameH-1,:))/2;
        
        % scale all the depth values
        % depthP = depthP*depth_mul;
        depthX = depthX*depth_mul;
        depthY = depthY*depth_mul;
        
        % Extrapolating depthP from depthX and depthY
        for yCounter = 1:frameH
            for xCounter = 1:frameW
                % Find current depthX and depthY
                currDepthX = depthX(yCounter, xCounter);
                currDepthY = depthY(yCounter, xCounter);

                % Find leftDepthX and downDepthY
                xLeft = xCounter - 1;
                yDown   = yCounter + 1;

                if xLeft < 1
                    leftDepthX = depthX(midY, tubeEndX); % openSpaceDepth;
                else
                    leftDepthX = depthX(yCounter, xLeft); 
                end

                if yDown > frameH
                    upDepthY = depthY(midY, tubeEndX); % openSpaceDepth;
                else
                    upDepthY = depthY(yDown, xCounter);
                end

                % Find the average
                avgDepth = (currDepthX+currDepthY+leftDepthX+upDepthY)/4;

                % Assign the average depth to depthP
                depthP(yCounter, xCounter)= avgDepth;
            end        
        end

        % Assign minDepth if anywhere depth is less than minDepth
        depthP(depthP<minDepth) = minDepth;
        depthX(depthX<minDepth) = minDepth;
        depthY(depthY<minDepth) = minDepth;
    end
           
    % Plot the depth of the tube
     plotTubeDepth(depthP);
     
    % Set the source
    getRadius = (tubeSectionDiameterCells(1)-1)/2;
    
    excitationXstart = tubeStartX-1;
    excitationXend   = excitationXstart;
    
    excitationYstart = tubeStartY - getRadius;
    excitationYend   = tubeStartY + getRadius;
    
    PV_N(excitationYstart:excitationYend, excitationXstart:excitationXend,4) = ...
         cell_excitation;
    
    % Set the cells just above the source as the cell_wall
    PV_N(excitationYstart-1,excitationXstart,4) = cell_wall;
    PV_N(excitationYend+1,excitationXstart,4) = cell_wall;
    
    % Set the noPressure cell
    getRadius = (tubeSectionDiameterCells(numSections)-1)/2;
    
    noPressureXstart  = tubeEndX+1;
    noPressureXend    = noPressureXstart;
    noPressureYstart  = tubeEndY - getRadius -1;
    noPressureYend    = tubeEndY + getRadius +1;
    
    PV_N(noPressureYstart:noPressureYend, noPressureXstart:noPressureXend,4) = ...
         cell_noPressure;
     
    % Set the listener
    listenerX = tubeEndX;
    listenerY = tubeEndY;   
end