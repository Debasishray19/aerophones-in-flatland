function plotTubeDepth(depthP)
    
    % Size of depthP matrix
    xLen = size(depthP,2);
    yLen = size(depthP,1);
    
    % Index matrix
    [yNew, xNew] = meshgrid(1:xLen, 1:yLen);
    
    % Plot the depthP
    plotObj = surf(xNew,yNew,depthP); % hold on;
    %surf(xNew,yNew,depthPNeg); hold off;
    
    plotObj.EdgeColor = 'none';
    shading interp
end