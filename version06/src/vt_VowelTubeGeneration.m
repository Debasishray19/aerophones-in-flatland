% This code implements the following paper to illustrate the tube cross
% sectional area for the vowel sound - /a/, /u/, /i/ : 
% Comparision of Magnetic resonance Imaging-Based vocal tract area function
% obtained from the same speaker in 1994 and 2002.

% We will be using 44 cross sections to generate the vocal tract
% Important: Sectional/segment Length(delta) is different from spatial resolution (ds)

function [listenerX, listenerY, glottalDiameter2D, glottalDiameter3D, upperTubeEndY, lowerTubeEndY, tubeEndX, frameH, frameW, mu2D, boundarySegmentType, depthX, depthY, depthP, baffleSwitch, micXposCells, glottalArea, PV_N]...
         = vt_VowelTubeGeneration(pmlSwitch, ds, vowelSound, simulation2D, cell_wall, cell_air, cell_excitation, cell_noPressure, cell_head)

    % Define units
    meter = 1;
    centimeter = 1e-2*meter;
    milimeter = 1e-3*meter;
    cellHalfLen = ds/2;
	
    % Vocal tract parameters
    numSections = 44;
    diameter_mul = 0.85;
    depth_mul = 1;
    micXpos = 3*milimeter;
    minDepth  = 0.001; % Define a minDepth
    
    % Find the diameter at the glottis
    glottalDiameter2D = 1;
    glottalDiameter3D = 1;
    
    % STEP0: Swich ON/OFF the circular baffle walls around the vocal tract 
    % based on the vowelSound. Or consider user input for the circular baffle
    % only in the case of vowelSound /a/, /i/ and /u/.
    if vowelSound ~=0
        baffleSwitch = input('Swicth ON the circular baffle. Press 1:ON 0:OFF = ');
    else
        % Define cell_wall as cell_air if it's an open space simulation [No vowel sound]
        cell_wall = cell_air;
        baffleSwitch=0;
    end
    
    % STEP1: Select the cross-sectional area based on the vowelSound 
    % Selct the cross sectional area based on the vowel sound
    % Tube section area in cm^2 in 3D
    switch vowelSound
        case 0 % For open air-space
            % [Note]: For open air-space, we will create the tube structure
            % first. Then we will switch the cell_wall to cell_air.
            sectional_length = 0.00388; % in meter
            tubeSectionArea_incm2_3D = [0 0 0 0 ...
                                        0 0 0 0 ...
                                        0 0 0 0 ...
                                        0 0 0 0 ...
                                        0 0 0 0 ...
                                        0 0 0 0 ...
                                        0 0 0 0 ...
                                        0 0 0 0 ...
                                        0 0 0 0 ...
                                        0 0 0 0 ...
                                        0 0 0 0];
                                    
        case 1 % For vowel /a/
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
        case 2 % For vowel /u/
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
            
        case 3 % For vowel /o/
            sectional_length = 0.00417; % in meter
            tubeSectionArea_incm2_3D = [0.38 0.45 0.57 0.77 ...
                                        1.31 1.92 1.74 1.11 ...
                                        0.75 0.59 0.57 0.68 ...
                                        0.73 0.67 0.58 0.49 ...
                                        0.44 0.42 0.49 0.53 ...
                                        0.38 0.30 0.45 0.61 ...
                                        0.71 0.79 0.86 1.01 ...
                                        1.41 2.09 3.00 4.10 ...
                                        5.16 6.22 7.34 8.15 ...
                                        8.61 8.37 6.76 4.37 ...
                                        2.30 1.06 0.58 0.47];
        case 4 % For vowel /i/
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
        case 5
            % Section length of each segment
            sectional_length = 0.00230; % in meter
            
            % Microphone Position - in meter
            mic_pos1 = 2*0.0091+0.003; % First microphone position 
            mic_pos2 = mic_pos1+0.01;  % Second microphone position 
            
            % Microphone Position - in cells
            mic_pos1_cells = round(mic_pos1/ds); % First microphone position 
            mic_pos2_cells = round(mic_pos2/ds);  % Second microphone position 
            
            tubeSectionArea_incm2_3D = [0.91 0.91 0.91 0.91 ...
                                        0.91 0.91 0.91 0.91 ...
                                        0.91 0.91 0.91 0.91 ...
                                        0.91 0.91 0.91 0.91 ...
                                        0.91 0.91 0.91 0.91 ...
                                        0.91 0.91 0.91 0.91 ...
                                        0.91 0.91 0.91 0.91 ...
                                        0.91 0.91 0.91 0.91 ...
                                        0.91 0.91 0.91 0.91 ...
                                        0.91 0.91 0.91 0.91 ...
                                        0.91 0.91 0.91 0.91];
            
        otherwise
            disp('Something Wrong - Check Your Code');
            return;
    end
                             
    % Tube section area in m^2
    tubeSectionArea_inm2_3D = tubeSectionArea_incm2_3D.*(centimeter*centimeter);    
    glottalArea = tubeSectionArea_inm2_3D(1);
    
    % STEP2: Calculate the tube section diameter with diameter multiplier
    tubeSectionDiameter_3D = 2*sqrt(tubeSectionArea_inm2_3D./pi).*diameter_mul;
    
    % STEP3: Convert the 3D cross sectional area into 2D using expansion
    % ratio. But we do not need this for 2.5D. Use a variable to switch
    % between 2.5D and 2D.

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
             
    % If simulation2D is 1 then we will be running 2D simulator otherwise 2.5D
    % will be running
    % Tube section diameter in terms of number of grid cell
    if simulation2D ==0
        tubeSectionDiameterCells = round(tubeSectionDiameter_3D./ds); 
        finalTubeSectionDiameter = tubeSectionDiameter_3D;
    else
        tubeSectionDiameterCells = round(tubeSectionDiameter_2D./ds);
        finalTubeSectionDiameter = tubeSectionDiameter_2D;
    end
    
    % Storing the boundary admittacnce value
    if vowelSound == 5
        mu3D = ones(1,numSections).*0.03;
    else
        mu3D = ones(1,numSections).*0.005; %(Just for testing vowel sound) 0.045
    end
    
    areaExpansionFlag = input('Apply area expansion ratio to boundary admittance [1-ON 0-Off]: ');
    if areaExpansionFlag ==1
        glottalDiameter2D = tubeSectionDiameter_2D(1);
        glottalDiameter3D = tubeSectionDiameter_3D(1);
        mu2D = (mu3D*2).*(tubeSectionDiameter_2D./tubeSectionDiameter_3D);
    else
        mu2D = (mu3D.*(2*0.5*pi/1.84));
    end
    

    % Change the tube diameter to 1 if it contains 0
    tubeSectionDiameterCells(tubeSectionDiameterCells==0)=1;
    
    %STEP4: Choose the best possible odd number from the Diameter array
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
  
    % STEP5: Find the total tube length and calculate the percentage error 
    % in the approximated tube length

    % Number of cells for total tube length
    actualTubeLength = numSections*sectional_length;
    totalTubeLengthinCells = round(actualTubeLength/ds);
    approxTubeLength = totalTubeLengthinCells*ds;
    
    % Percentage error in total tube length
    totalTubeLengthError = abs(actualTubeLength - approxTubeLength)/...
                           (actualTubeLength);
    fprintf('Error percentage in approximated tube length = %.4f \n',totalTubeLengthError*100);
    
    % STEP6: Determine circular baffle radius 
    diameterBaffleVar = 0.1800; % User Input to change the baffle diameter
    baffleDiameter = diameterBaffleVar;
    
    % Find the baffle radius (r)
    baffleRadius = baffleDiameter/2;
    
    % Find the height of the baffle at the mouth opening from the mid-row
    % cells of vocal tract (p)
    
    % Get the radius value at the tube end [lips]
    getRadiusVal = (tubeSectionDiameterCells(numSections)-1)/2;
    
    % baffleHeight[at the end of the tube] = radius+boundaryCell+1
    baffleHeight = (getRadiusVal*ds)+(1*ds)+ds; 
    
    % Find the Centre Cell
    baffleBase = sqrt(baffleRadius^2 - baffleHeight^2);
    baseNumCells = round(baffleBase/ds);
    
    % STEP7: Construct the domain, frame and PV_N array
    offsetW = 20;
    offsetH = 10;
    
    if baffleSwitch == 1
        domainW = 2*baseNumCells + offsetW;
        domainH = 2*baseNumCells + offsetH;
    else
        domainW = totalTubeLengthinCells + offsetW;
        if vowelSound == 0 % For vowelSound=0, set increase the domain height
            domainH = max(tubeSectionDiameterCells) + offsetH+50;
        else
            domainH = max(tubeSectionDiameterCells) + offsetH;
        end
    end
    
    % Build frame [Frame = Domain Size + PML Layers]
    if pmlSwitch == 1
        % fOR SRATE_MUL: 15 frameW = 335 frameH = 268 (with circular baffle)
        % fOR SRATE_MUL: 20 frameW = 350 frameH = 350 (with circular baffle)
        frameW = 300; % domainW+2*pmlLayer+2; % 2 is for border/dead cells
        frameH = 268; % domainH+2*pmlLayer+2; % 2 is for border/dead cells
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
    
        
    % Define each boundary segment belongs to which type [Out Of 1-44]
    boundarySegmentType = zeros(frameH, frameW);
    
    % STEP8: Define the depth matrix to store the tube height along z-axis
    if simulation2D == 1
        openSpaceDepth = 1;
    else
        openSpaceDepth = input('Enter the open space depth value: ');
    end
    
    depthP    = ones(frameH, frameW)*openSpaceDepth*meter;
    depthX    = ones(frameH, frameW)*openSpaceDepth*meter;
    depthY    = ones(frameH, frameW)*openSpaceDepth*meter;
       
    % Array to store the tube upper and lower wall positions
    % This array will help us to major the depth/height across the tube
    % Need to calculate depth/height for each cell in a column
    % How many columns = totalTubeLengthinCells
    % ROW1=Upper Wall Idx, ROW2=Lower Wall Idx, ROW3=Radius
    cellDepthProp = zeros(2, totalTubeLengthinCells);
    
    % STEP9: Find the mid point in the frame
    midY = floor(frameH/2);
%   midX = floor(frameW/2);
    
    % Get radius at the end of the tube
    tubeEndRadius = (tubeSectionDiameterCells(numSections)-1)/2;
        
    % STEP10: Find the glottal end: Starting point of the tube
    %  tubeStartX = midX - round(totalTubeLengthinCells/2);
    % Consider the best position for tubeStartX, not the mid of frame.
    % This will reduce the unnecessary domain size.
    % [Note] subtract 10 from tubeStartX, if you are simulating vowel
    % sounds with no PML layer or no circular baffle
    % subtract 35, ifyou are simulating same vowel sound with srate=15, PML
    % switch on and circular baffle switch on.
    tubeStartX = frameW-totalTubeLengthinCells-10;
    tubeStartY = midY;
    tubeEndX   = tubeStartX + totalTubeLengthinCells-1;
    tubeEndY   = midY;
    
    % Find the y-coordinate for the both upper and lower cell at the end of
    % the VT tube and just aboube the tube wall.
    upperTubeEndY = tubeEndY-tubeEndRadius-2;
    lowerTubeEndY = tubeEndY+tubeEndRadius+2;
        
    % STEP11: Set the circular baffle walls aound the vocal tract
    if baffleSwitch==1
        baffleCentreY = midY;
        baffleCentreX = tubeEndX-baseNumCells;

        % Assign cell_head to the cells just aboove and below the end of the VT 
        % tube boundaries   
        PV_N(upperTubeEndY, tubeEndX, 4) = cell_head;
        PV_N(lowerTubeEndY, tubeEndX, 4) = cell_head;

        % Create the circular baffle
        for columnLoop = tubeEndX-1:-1:1        
            % baffle base length for the column
            baffleBaseLen = abs(columnLoop - baffleCentreX)*ds;

            if baffleBaseLen <= baffleRadius
                % Find the baffle height = sqrt(bafflleR^2 - baffleBaseLen^2)
                baffleHeightLen = sqrt(baffleRadius^2 - baffleBaseLen^2);

                % Find the height
                baffleHeightCells = round(baffleHeightLen/ds);
                upperBaffleCellY = baffleCentreY - baffleHeightCells;
                lowerBaffleCellY = baffleCentreY + baffleHeightCells;
                
                % Set the cell as cell head
                PV_N(upperBaffleCellY, columnLoop,4) = cell_head;
                PV_N(lowerBaffleCellY, columnLoop,4) = cell_head;        
            end
        end

        % Connect the disjoint circular baffle
        % From tube end/mouth to centre of the baffle
        for columnLoop = tubeEndX:-1:baffleCentreX-1
            currColumnHeadPositions = find(PV_N(:,columnLoop,4)==cell_head);
            prevColumnHeadPositions = find(PV_N(:,columnLoop-1,4)==cell_head);

            PV_N(prevColumnHeadPositions(1)+1:currColumnHeadPositions(1), columnLoop, 4) = cell_head;
            PV_N(currColumnHeadPositions(2):prevColumnHeadPositions(2)-1, columnLoop, 4) = cell_head;
        end

        % Find the column from where circular buffer starts
        bufferStartX = 1;

        while isempty(find(PV_N(:,bufferStartX,4)==cell_head,1))
            bufferStartX =  bufferStartX+1;
        end

        % Fill the bufferStart column
        bufferStartPositions = find(PV_N(:,bufferStartX,4)==cell_head);

        for columnLoop = bufferStartX:baffleCentreX-1
            currColumnHeadPositions = find(PV_N(:,columnLoop,4)==cell_head);
            nextColumnHeadPositions = find(PV_N(:,columnLoop+1,4)==cell_head);

            PV_N(nextColumnHeadPositions(1)+1:currColumnHeadPositions(1), columnLoop, 4) = cell_head;
            PV_N(currColumnHeadPositions(2):nextColumnHeadPositions(2)-1, columnLoop, 4) = cell_head;
        end
        PV_N(bufferStartPositions(1):bufferStartPositions(2), bufferStartX, 4) = cell_head;
    end
    
    % STEP12: Store the cummulative length of each tube section
    tubeCummSectionLength = zeros(1, numSections);
    for sectionCount = 1:numSections
        tubeCummSectionLength(sectionCount) = sectional_length*sectionCount;
    end
    
    % STEP13: Draw the geometrical shape of the tube
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
        
        % Set the wall segmentation type here
        boundarySegmentType(upperY, upperX)= sectionCounter;
        boundarySegmentType(lowerY, lowerX)= sectionCounter;
        
        if sectionCounter<45
            boundarySegmentType(upperY, upperX)= sectionCounter;
            boundarySegmentType(lowerY, lowerX)= sectionCounter;
        else
            boundarySegmentType(upperY, upperX)= 44;
            boundarySegmentType(lowerY, lowerX)= 44;
        end
        
        % Store in the cellDepthProp to calculate cell depth
        cellDepthProp(1, tubeCellsCount) = upperY;
        cellDepthProp(2, tubeCellsCount) = lowerY;
            
        % Find the difference in prevUpperY and current upperY
        differYwall = upperY - prevUpperY;
        if tubeCellsCount~=1
            if differYwall>1
                % Set upper wall
                PV_N(prevUpperY:upperY, upperX,4) = cell_wall;
                boundarySegmentType(prevUpperY:upperY, upperX)= sectionCounter;
                
                % Set lower wall
                PV_N(lowerY:prevLowerY, lowerX,4) = cell_wall;
                boundarySegmentType(lowerY:prevLowerY, lowerX)= sectionCounter;
                
            elseif differYwall < -1
                % Set upper wall
                PV_N(upperY:prevUpperY, prevUpperX,4) = cell_wall;
                boundarySegmentType(upperY:prevUpperY, prevUpperX)= sectionCounter;
                
                % Set lower wall
                PV_N(prevLowerY:lowerY, prevLowerX,4) = cell_wall;
                boundarySegmentType(prevLowerY:lowerY, prevLowerX)= sectionCounter;
            end
        end
        
        prevUpperY = upperY;
        prevUpperX = upperX;
        prevLowerY = lowerY;
        prevLowerX = lowerX;
    end % End of For loop
    
    % STEP14: Set the depthX, depthY and depthP if simulation2D~=1. 
    % For 2D simulation we do not need depth parameter.
    % Calculate depth/height of the tube along z-axis for depthX and depthY
    % For the boundary set the radius as zero
    if simulation2D~=1 && vowelSound ~=0
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
        % last segment of the tube (if the circular baffle is not on)
        if baffleSwitch==0
            depthX(depthX==openSpaceDepth) = depthX(midY, tubeEndX);
            depthY(depthY==openSpaceDepth) = depthY(midY, tubeEndX);
        end
                
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
    
    % Set the circular baffle cells to minDepth
    depthP(PV_N(:,:,4)==cell_head)=minDepth;
    depthX(PV_N(:,:,4)==cell_head)=minDepth;
    depthY(PV_N(:,:,4)==cell_head)=minDepth;
    
    % Plot the depth of the tube
%      plotTubeDepth(depthP);
     
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
    
    % Define the grid column at the tube end as cell_noPressure to
    % implement Dirchilet Boundary Condition
    if baffleSwitch ~= 1 && vowelSound~=0
        PV_N(noPressureYstart:noPressureYend, noPressureXstart:noPressureXend,4) = ...
                cell_noPressure;
    end
     
    % Find the distance of the mic position from the noPressure cells
    % (Inside the tube)   
    micXposCells = round(micXpos/ds);
    
    % Define the place of Listener
    % We can position our mic based upon the experiment thet we are running
    
    switch baffleSwitch       
        case 0 % Exp1: Simmulationg without a circular baffle             
            listenerX = tubeEndX-micXposCells;
            listenerY = tubeEndY;
            
        case 1 % Exp2: Simulating with circular baffle
            listenerX = tubeStartX + mic_pos1_cells;
            listenerY = tubeEndY;
        otherwise
    end   
end