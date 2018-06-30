function [refFrameBeta, xSrc, ySrc] = buildFrameBeta(domainW, domainH, tubeHorizontalLength,...
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
    xSrc = Nx/2;
    ySrc = Ny/2;
    % STEP1: Create a single reference frame to define PML, Tube wall , source. 
    % And Initialize every cell of the frame as Air with beta = 1;
    refFrameBeta= ones(Nx, Ny);
    
    % Note: No need to separate PML layers and air space here, because both 
    % the PML layers and air have the same beta value.
    
    % STEP2: Build Tube : Sigma value for the tube wall should be 
    % [1] Check if tube has any horizontal length, but no vertical length
    % [2] Check if tube has any vertical length, but not horizontal length
    % [3] Check if tube has both horizontal and vertical length

    % Begin: 
    
    % [1] : Check if the tube has any horizontal length, but no vertical length  
    if tubeHorizontalLength > 0 && tubeVerticalLength == 0
        xStartUpperWall = floor(Nx/2) - round(tubeHorizontalLength/2);
        xEndUpperWall   = floor(Nx/2) + (tubeHorizontalLength - round(tubeHorizontalLength/2));
        
        xStartLowerWall = xStartUpperWall;
        xEndLowerWall   = xEndUpperWall;
        
        yStartUpperWall = floor(Ny/2);
        yEndUpperWall = floor(Ny/2);
        
        yStartLowerWall = yStartUpperWall + tubeWidth + 1;
        yEndLowerWall   = yEndUpperWall   + tubeWidth + 1;
        
        % Build upper wall
        refFrameBeta (xStartUpperWall:xEndUpperWall,...
                     yStartUpperWall: yEndUpperWall) = betaTube;
        
        % Build lower wall
        refFrameBeta (xStartLowerWall:xEndLowerWall,...
                     yStartLowerWall: yEndLowerWall) = betaTube;
                 
        % Build left side wall
        refFrameBeta (xStartLowerWall,...
                     yStartUpperWall: yStartLowerWall) = betaTube;
                 
        % Find source coordinate inside the tube
        xSrc = xStartLowerWall+1;
        ySrc = round((yStartUpperWall + yStartLowerWall)/2);
    end
   
    % [2] Check if the tube has any vertical length, but not horizontal length
    if tubeHorizontalLength == 0 && tubeVerticalLength > 0
        yStartLeftWall  = floor(Ny/2) - round(tubeVerticalLength/2);
        yEndLeftWall    = floor(Nx/2) + (tubeVerticalLength - round(tubeVerticalLength/2));
        
        yStartRightWall = yStartLeftWall;
        yEndRightWall = yEndLeftWall;
        
        xStartLeftWall = floor(Nx/2);
        xEndLeftWall   = floor(Nx/2);
        
        xStartRightWall =  xStartLeftWall + tubeWidth + 1;
        xEndRightWall   =  xEndLeftWall + tubeWidth + 1;
        
        % Build left wall
        refFrameBeta(xStartLeftWall:xEndLeftWall,...
                     yStartLeftWall:yEndLeftWall) = betaTube;
                 
        % Build right wall
        refFrameBeta(xStartRightWall:xEndRightWall,...
                     yStartRightWall:yEndRightWall) = betaTube;  
                 
        % Build lower wall
        refFrameBeta(xEndLeftWall:xEndRightWall,...
                     yEndLeftWall) = betaTube;
                 
        % Find source coordinate inside the tube
        xSrc = round((xEndLeftWall + xEndRightWall)/2);
        ySrc = yEndLeftWall-1;
    end
    
    return; 
end