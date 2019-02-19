% This code implements the following paper to illustrate the tube cross
% sectional area for the vowel sound - /a/ : 
% Comparision of Magnetic resonance Imaging-Based vocal tract area function
% obtained from the same speaker in 1994 and 2002.

% We will be using 44 cross sections to generate the vocal tract
% Important: Sectional Length is different from spatial resolution (ds)

function [listenerX, listenerY, frameH, frameW, PV_N]...
         = aaCrossSectionTube(pmlSwitch, ds, pmlLayer, cell_wall, cell_air, cell_excitation, cell_noPressure)
     
    % Define units
    meter = 1;
    centimeter = 1e-2*meter;
    
    % Vocal tract parameters
    sectional_length = 0.00388; % in meter
    numSections = 44;
    
    %STEP0: Convert the 3D area function vecor to 2D
    % Tube section area in cm^2 in 3D
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
                             
    % Tube section area in m^2
    tubeSectionArea_inm2_3D = tubeSectionArea_incm2_3D.*(centimeter*centimeter);    
    
    % Tube section diameter
    tubeSectionDiameter_3D = 2*sqrt(tubeSectionArea_inm2_3D./pi);
        
    % Tube section diameter in cm^2 in 2D
    % Create the array 
    tubeSectionDiamete_2D = zeros(1, numSections);
    
    % Find out the maxium diameter in 3D and it's index
    [maxDiameter, maxDiameterIdx] = max(tubeSectionDiameter_3D);
    
    % Assign the updated diameter to the tubeSectionDiamete_2D for the same
    % index position: d_2D = d3D(0.5*pi/1.84)   
    tubeSectionDiamete_2D(maxDiameterIdx) = (maxDiameter*0.5*pi)/1.84;
    
    % To Assign other diameters for tubeSectionDiamete_2D implement
    % eqn -(3a) and (3b) from the foloowing paper:
    % Two dimensional vocal tracts with three dimensional behaviour in the
    % numerical generation of vowels.
    
    % Store the value of m in an array - are expansion ratio
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
        tubeSectionDiamete_2D(radiiLessMaxIdx-1) = ...
        tubeSectionDiamete_2D(radiiLessMaxIdx) * m(radiiLessMaxIdx-1);
    
        radiiLessMaxIdx = radiiLessMaxIdx - 1;
    end
    
    radiiGreaterMaxIdx = maxDiameterIdx;
    while radiiGreaterMaxIdx < numSections
        tubeSectionDiamete_2D(radiiGreaterMaxIdx+1) = ...
        tubeSectionDiamete_2D(radiiGreaterMaxIdx) * m(radiiGreaterMaxIdx+1);
    
        radiiGreaterMaxIdx = radiiGreaterMaxIdx + 1;
    end
       
    % Tube section diameter in terms of number of grid cell
    tubeSectionDiameterCells = round(tubeSectionDiamete_2D./ds);
    
    % Change the tube diameter to 1 if it contains 0
    tubeSectionDiameterCells(tubeSectionDiameterCells==0)=1;
    
    %STEP1: Choose the best possible odd number from the Diameter array
    for diameterCounter = 1:numSections       
        % Verify if the cellsPerDiameter is odd or not
        if mod(tubeSectionDiameterCells(diameterCounter), 2) == 0
            
            % Find the difference between rounded and actual diameter value
            diff = tubeSectionDiameterCells(diameterCounter) - ...
                    (tubeSectionDiamete_2D(diameterCounter)/ds);            
            
            if diff>0
                tubeSectionDiameterCells(diameterCounter) = ...
                    tubeSectionDiameterCells(diameterCounter)-1;
            else
                tubeSectionDiameterCells(diameterCounter)=...
                    tubeSectionDiameterCells(diameterCounter)+1;
            end
        end   
    end
       
    % STEP2: Find the total tube length and calculate the percentage error 
    % in the approximated tube length

    % Number of cells for total tube length
    actualTubeLength = 44*sectional_length;
    totalTubeLengthinCells = round(actualTubeLength/ds);
    approxTubeLength = totalTubeLengthinCells*ds;
    
    % Percentage error in total tube length
    totalTubeLengthError = abs(actualTubeLength - approxTubeLength)/...
                           (actualTubeLength);
    fprintf('Error percentage in approximated tube length = %.4f \n',totalTubeLengthError*100);
    
    % STEP2: Construct the frame
    offsetW = 6;
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
        
    % STEP2: Find the mid point in the frame
    midY = floor(frameH/2);
    midX = floor(frameW/2);
    
    % STEP3: Find the glottal end: Starting point of the tube
    tubeStartX = midX - round(totalTubeLengthinCells/2);
    tubeStartY = midY;
    tubeEndX   = tubeStartX + totalTubeLengthinCells-1;
    tubeEndY   = midY;
    
    % STEP4: Store the cummulative length of each tube section
    tubeCummSectionLength = zeros(1, numSections);
    for sectionCount = 1:numSections
        if sectionCount == 1
           tubeCummSectionLength(sectionCount) =  sectional_length;
        else
            tubeCummSectionLength(sectionCount) = sectional_length + tubeCummSectionLength(sectionCount-1);
        end
    end
    
    % STEP5: Draw the geometrical shape of the tube
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
        
        % We are subtracting 1 as we'll assume that there is a middle
        % row of cells which will act like a morror.
        getRadius = (tubeSectionDiameterCells(sectionCounter)-1)/2;
            
        % Find the upper and lower wall coordinates
        upperX = tubeStartX + (tubeCellsCount-1);
        upperY = tubeStartY - getRadius - 1;           
        lowerX = tubeStartX + (tubeCellsCount-1);
        lowerY = tubeStartY + getRadius + 1;
        
        % Verify if the currTubeLength is more than the tubeCummSectionLength
        % for the current section counter
        
        % if small or equal then set the tube wall as expected-Normal Case
        if currTubeLength <= tubeCummSectionLength(sectionCounter)            
            % set the upper wall & lower wall
            PV_N(upperY, upperX,4) = cell_wall;
            PV_N(lowerY, lowerX,4) = cell_wall;
        else
        % If the currTubeLength is greater than the actual tubeCummSectionLength
            % Find the difference between currTubeLength and tubeCummSectionLength
            diffLength = currTubeLength - tubeCummSectionLength(sectionCounter);
            
            if diffLength>0.5 && sectionCounter~=numSections
                % Increase the section counter
                sectionCounter = sectionCounter+1;
                
                % Get the radius for that cross-section
                getRadius = (tubeSectionDiameterCells(sectionCounter)-1)/2;
                
                % Find the upper and lower wall coordinates
                upperX = tubeStartX + (tubeCellsCount-1);
                upperY = tubeStartY - getRadius - 1;           
                lowerX = tubeStartX + (tubeCellsCount-1);
                lowerY = tubeStartY + getRadius + 1;
                
                % Set the upper and lower wall
                PV_N(upperY, upperX,4) = cell_wall;
                PV_N(lowerY, lowerX,4) = cell_wall;
            else
                % Set the upper and lower wall
                PV_N(upperY, upperX,4) = cell_wall;
                PV_N(lowerY, lowerX,4) = cell_wall;
                
                % And then increase the section Counter
                sectionCounter = sectionCounter+1;
            end            
        end  % End of checking currTubeLength and tubeCummSectionLength
                
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
    getRadius = (tubeSectionDiameterCells(44)-1)/2;
    
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