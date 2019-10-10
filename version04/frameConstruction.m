function [PV_N, frameH, frameW, depthX, depthY, depthP] = ...
         frameConstruction(domainH, domainW, pmlSwitch, pmlLayer, simulation2D)
    
    % Define unit
    meter = 1;
    
    % Build frame [Frame = Domain Size + PML Layers]
    if pmlSwitch == 1
        frameW = domainW+2*pmlLayer+2; % 2 is for border/dead cells
        frameH = domainH+2*pmlLayer+2; % 2 is for border/dead cells
    else
        frameW = domainW+2;
        frameH = domainH+2;
    end
    
    % 2D simulation = set openSpaceWidth to 1
    % 2.5D simulation = Take user input for setting the openSpaceWidth
    if simulation2D == 1
        openSpaceDepth = 1;
    else
        openSpaceDepth = input('Enter a default tube depth value: ');
    end
    
    depthP    = ones(frameH, frameW)*openSpaceDepth*meter;
    depthX    = ones(frameH, frameW)*openSpaceDepth*meter;
    depthY    = ones(frameH, frameW)*openSpaceDepth*meter;
    
    % For pressure and velocity
    % PV_N(:,:,1) = To store cell pressure
    % PV_N(:,:,2) = To store Vx
    % PV_N(:,:,3) = To store Vy
    % PV_N(:,:,4) = To store grid cell type
    PV_N = zeros(frameH, frameW, 4); 
end