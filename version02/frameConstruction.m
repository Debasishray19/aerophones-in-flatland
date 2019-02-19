function [PV_N, frameH, frameW] = frameConstruction(domainH, domainW, pmlSwitch, pmlLayer)

    % Build frame [Frame = Domain Size + PML Layers]
    if pmlSwitch == 1
        frameW = domainW+2*pmlLayer+2; % 2 is for border/dead cells
        frameH = domainH+2*pmlLayer+2; % 2 is for border/dead cells
    else
        frameW = domainW+2;
        frameH = domainH+2;
    end
    
    % For pressure and velocity
    % PV_N(:,:,1) = To store cell pressure
    % PV_N(:,:,2) = To store Vx
    % PV_N(:,:,3) = To store Vy
    % PV_N(:,:,4) = To store grid cell type
    PV_N = zeros(frameH, frameW, 4); 
end