function [refFrameBeta, xSrc, ySrc, xLis, yLis] = buildFrameBeta(domainW, domainH, tubeHorizontalLength,...
                   tubeVerticalLength, tubeWidth, pmlLayer)
            
    % BUILDFRAMEBETA function helps to add the beta value to the grid cells
    % Beta value = 1 for PML layers at the absorbing boundary
    % Beta value = 1 for Air
    % Beta value = 0 for Tube Wall     
    % Beta Value = 1 for listener
     betaTube = 0;

    % Frame Size = Domain Size + PML Layer(Up/Left + Down/Right)
    Nx = domainW + 2*pmlLayer;
    Ny = domainH + 2*pmlLayer;
    
    % Default xSrc and ySrc value
    xSrc = ceil(Nx/2);
    ySrc = ceil(Ny/2);
    % STEP1: Create a single reference frame to define PML, Tube wall , source. 
    % And Initialize every cell of the frame as Air with beta = 1;
    refFrameBeta= ones(Nx, Ny);
    
    % Note: No need to separate PML layers and air space here, because both 
    % the PML layers and air have the same beta value.
    
    % STEP2: Build Tube : Sigma value for the tube wall should be 
    % [1] Check if the tube has any horizontal length, but no vertical length
    % [2] Check if the tube has any vertical length, but not horizontal length
    % [3] Check if the tube has both horizontal and vertical length

    % Begin: 
    
    % [1] : Check if the tube has any horizontal length, but no vertical length   
    if tubeHorizontalLength > 0 && tubeVerticalLength == 0
        
		% Find the mid position
		column_mid = xSrc;
		row_mid = ySrc;
		
		% Tube wall coordinates
		tubeStartRowPos_UpperWall = row_mid;
		tubeStartColumnPos_UpperWall = column_mid;
		
		tubeEndRowPos_UpperWall = tubeStartRowPos_UpperWall;
		tubeEndColumnPos_UpperWall = tubeStartColumnPos_UpperWall ...
                                     + tubeHorizontalLength;
				
		tubeStartRowPos_LowerWall = tubeStartRowPos_UpperWall + tubeWidth;
		tubeStartColumnPos_LowerWall = tubeStartColumnPos_UpperWall;
        
        tubeEndRowPos_LowerWall = tubeStartRowPos_LowerWall;
		tubeEndColumnPos_LowerWall = tubeEndColumnPos_UpperWall;
		
		% Build upper tube wall
		refFrameBeta(tubeStartRowPos_UpperWall:tubeEndRowPos_UpperWall,...
             tubeStartColumnPos_UpperWall:tubeEndColumnPos_UpperWall) = ...
                betaTube;
		
		% Build lower tube wall
        refFrameBeta(tubeStartRowPos_LowerWall:tubeEndRowPos_LowerWall,...
             tubeStartColumnPos_LowerWall:tubeEndColumnPos_LowerWall) = ...
                betaTube;
            
        % Build side tube wall
        refFrameBeta(tubeStartRowPos_UpperWall:tubeStartRowPos_LowerWall,...
             tubeStartColumnPos_UpperWall:tubeStartColumnPos_LowerWall) = ...
             betaTube;
        
        % Find source position
        xSrc = tubeStartRowPos_UpperWall + ceil(tubeWidth/2);
        ySrc = tubeStartColumnPos_UpperWall + 1;
        
        % Find listener position
        xLis = xSrc;
        yLis = tubeEndColumnPos_UpperWall + 2;
    end
    
    % [2] Check if the tube has any vertical length, but not horizontal length
    if tubeHorizontalLength == 0 && tubeVerticalLength > 0
        % Find the mid position
		column_mid = xSrc;
		row_mid = ySrc;
        
        % Tube wall coordinates
        tubeStartRowPos_LeftWall = row_mid;
        tubeStartColumnPos_LeftWall = column_mid;
        
        tubeEndRowPos_LeftWall = tubeStartRowPos_LeftWall +...
                                 tubeVerticalLength;
        tubeEndColumnPos_LeftWall = tubeStartColumnPos_LeftWall;
        
        tubeStartRowPos_RightWall = tubeStartRowPos_LeftWall; 
        tubeStartColumnPos_RightWall = tubeStartColumnPos_LeftWall +...
                                       tubeWidth;
                                   
        tubeEndRowPos_RightWall = tubeEndRowPos_LeftWall;
        tubeEndColumnPos_RightWall = tubeStartColumnPos_RightWall;
        
        % Build Left Wall
        refFrameBeta(tubeStartRowPos_LeftWall:tubeEndRowPos_LeftWall,...
                 tubeStartColumnPos_LeftWall:tubeEndColumnPos_LeftWall)=...
                 betaTube;
                
        % Build Right Wall
        refFrameBeta(tubeStartRowPos_RightWall:tubeEndRowPos_RightWall,...
                 tubeStartColumnPos_RightWall:tubeEndColumnPos_RightWall)=...
                 betaTube;
                
        % Build Lower Wall
        refFrameBeta(tubeStartRowPos_LeftWall:tubeStartRowPos_RightWall,...
                 tubeStartColumnPos_LeftWall:tubeStartColumnPos_RightWall)=...
                 betaTube;
             
       % Find Source Position
       xSrc = tubeStartRowPos_LeftWall +1;
       ySrc = tubeStartColumnPos_LeftWall + ceil(tubeWidth/2);
       
       % Find Listener position
       xLis = tubeEndRowPos_LeftWall + 2;
       yLis = ySrc;                            
    end
    
    % [3] Check if the tube has both horizontal and vertical length
    if tubeHorizontalLength > 0 && tubeVerticalLength >0
        % Find the mid position
        column_mid = xSrc;
        row_mid = ySrc;
        
        % Tube wall coordinates
		tubeStartRowPos_UpperWall = row_mid;
		tubeStartColumnPos_UpperWall = column_mid;
		
		tubeEndRowPos_UpperWall = tubeStartRowPos_UpperWall;
		tubeEndColumnPos_UpperWall = tubeStartColumnPos_UpperWall ...
                                     + tubeHorizontalLength;
				
		tubeStartRowPos_LowerWall = tubeStartRowPos_UpperWall + tubeWidth;
		tubeStartColumnPos_LowerWall = tubeStartColumnPos_UpperWall + tubeWidth;
        
        tubeEndRowPos_LowerWall = tubeStartRowPos_LowerWall;
		tubeEndColumnPos_LowerWall = tubeEndColumnPos_UpperWall;
        
        tubeStartRowPos_LeftWall = tubeStartRowPos_UpperWall;
        tubeStartColumnPos_LeftWall = tubeStartColumnPos_UpperWall;
        
        tubeEndRowPos_LeftWall = tubeStartRowPos_LeftWall +...
                                 tubeVerticalLength;
        tubeEndColumnPos_LeftWall = tubeStartColumnPos_LeftWall;
        
        tubeStartRowPos_RightWall = tubeStartRowPos_LowerWall; 
        tubeStartColumnPos_RightWall = tubeStartColumnPos_LowerWall;
                                   
        tubeEndRowPos_RightWall = tubeEndRowPos_LeftWall;
        tubeEndColumnPos_RightWall = tubeStartColumnPos_RightWall;
        
        % Build upper tube wall
		refFrameBeta(tubeStartRowPos_UpperWall:tubeEndRowPos_UpperWall,...
             tubeStartColumnPos_UpperWall:tubeEndColumnPos_UpperWall) = ...
                betaTube;
		
		% Build lower tube wall
        refFrameBeta(tubeStartRowPos_LowerWall:tubeEndRowPos_LowerWall,...
             tubeStartColumnPos_LowerWall:tubeEndColumnPos_LowerWall) = ...
                betaTube;
            
        % Build Left Wall
        refFrameBeta(tubeStartRowPos_LeftWall:tubeEndRowPos_LeftWall,...
                 tubeStartColumnPos_LeftWall:tubeEndColumnPos_LeftWall)=...
                 betaTube;
                
        % Build Right Wall
        refFrameBeta(tubeStartRowPos_RightWall:tubeEndRowPos_RightWall,...
                 tubeStartColumnPos_RightWall:tubeEndColumnPos_RightWall)=...
                 betaTube;
                
        % Build Lower Wall
        refFrameBeta(tubeEndRowPos_LeftWall:tubeEndRowPos_RightWall,...
                 tubeEndColumnPos_LeftWall:tubeEndColumnPos_RightWall)=...
                 betaTube;
             
       % Find Source Position
       xSrc = tubeEndRowPos_LeftWall -1;
       ySrc = tubeEndColumnPos_LeftWall + ceil(tubeWidth/2);
       
       % Find Listener position
       xLis = tubeStartRowPos_UpperWall + ceil(tubeWidth/2);
       yLis = tubeEndColumnPos_UpperWall + 2;
    end
end