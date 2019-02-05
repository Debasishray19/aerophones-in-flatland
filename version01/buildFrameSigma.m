function [refFrameSigma] = buildFrameSigma(domainW, domainH, pmlLayer, maxSigmaVal, dt)
                         
    % BUILDFRAMESIGMA function helps to add the PML layers to the outer space  
    % of the domain  and store sigma value for each cell. The PML layers 
    % are defined on the top of domainW and domainH.
    % Sigma is the attenuation coefficient which must defined to introduce
    % PML layers in the frame. For all other grid cells the sigma value
    % must be zero.
    
    % Frame Size = Domain Size + PML Layer(Up/Left + Down/Right)
    Nx = domainW + 2*pmlLayer;
    Ny = domainH + 2*pmlLayer;
    
    % STEP1: Create a single reference frame to define PML. 
    % And Initialize every cell of the frame as Air with sigma =0.
    refFrameSigma = zeros(Nx, Ny);

    % STEP2: Define PML layers along X-Axis (Along horizontal)and Y-axis (Along vertical)
    for layerCount = 0:pmlLayer-1
        % define the height from the top and bottom
        yTop    = pmlLayer-layerCount; % To add layers from the top along horizontally
        yBottom = (Ny - pmlLayer + 1)+layerCount; % To add layers from the bottom along horizontally

        % define the range of x-Axis
        xLeft = pmlLayer-layerCount; % To add the layers from the left most position
        xRight = (Nx - pmlLayer + 1)+layerCount; % To add the layers from right most position

        refFrameSigma(xLeft:xRight,yTop) = (layerCount/dt)*maxSigmaVal;
        refFrameSigma(xLeft:xRight,yBottom) = (layerCount/dt)*maxSigmaVal;

        refFrameSigma(xLeft,yTop:yBottom) = (layerCount/dt)*maxSigmaVal;
        refFrameSigma(xRight,yTop:yBottom) = (layerCount/dt)*maxSigmaVal;
    end
    
    return;
end