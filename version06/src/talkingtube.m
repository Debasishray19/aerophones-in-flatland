% Author: Debasish Ray Mohapatra
% Date: 20 March, 2019
% A special thanks to Victor Zappi who helped me to understand and
% implement this code. To visit Victor Zappi's website: http://toomuchidle.com/

% Source Implementation
% To implement Sinusoidal wave or Impulse source, use: excitationV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE MATLAB ENVIRONMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; 
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE UNITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meter = 1;
centimeter  =1e-2 * meter;

second    = 1;
millisecond = 1e-3 * second;
hertz     = 1/second;
kilohertz = 1e3 * hertz;
megahertz = 1e6 * hertz;
gigahertz = 1e9 * hertz;

gram      = 1e-3;
kilogram  = 1e3*gram;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Toggle this variable to switch to 2.5D
% if set2D=1, then the 2D simulator will run otherwise 2.5D
simulation2D = input('Enter 1 for 2D simulation or 0 for 2.5D simulation: ');
rho = 1.140*kilogram/(meter^3);   % Air density  [kg/m^3]
srate_mul = input('Enter the sample rate multiplier: ');    % srate multiplier
c   = 350*meter/second;            % Sound speed in air [m/s]
maxSigmadt = 0.5;                  % Attenuation coefficient at the PML layer
srate = 44100*hertz*srate_mul;           % Sample frequency
pmlLayer = 6;                      % Number of PML layers
baffleSwitch = 0;                  % By default model should not have head/circular baffle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1/srate;                      % Temporal resolution/ sample time period
dx = dt*c*sqrt( 2.0 );             % Spatial resolution along x-direction: CFL Condition
dy = dt*c*sqrt( 2.0 );             % Spatial resolution along x-direction: CFL Condition
AudioTime = 1*second;              % Total audio signal time
kappa = rho*c*c;                   % Bulk modulus
ds = dx;                           % Spatial resolution(ds) = dx = dy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:dt:AudioTime-dt;            % time steps
STEPS = length(t);                % Total time steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE CELL TYPE [Not going to use 'cell_dynamics']
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cell_wall       = 0;
cell_air        = 1;
cell_excitation = 2;
cell_pml0       = 3;
cell_pml1       = 4;
cell_pml2       = 5;
cell_pml3       = 6;
cell_pml4       = 7;
cell_pml5       = 8;
cell_dynamic    = 9;
cell_dead       = 10;
cell_noPressure = 11; % To implement Dirichlet boundary condition
cell_head       = 12;
cell_numTypes   = 13;

vis_Boundary = 2000;
sigmadt = zeros(pmlLayer, 1);

% To store beta(tube wall) and sigmaPrimedt(PML Layers). 
% Beta for air & PMLLayers = 1 and for tubewall = 0
% sigmaPrimedt = sigmaPrime*dt
% sigmaPrime = 1 - Beta + Sigma
% e.g - 
% Sigma=0 for all the non-PML layers. Hence, sigmaPrime = 1 - Beta inside
% the domain.
% WALL -> beta = 0, sigma_prima*dt = (1-0)*dt = 1*dt = dt 
% AIR  -> beta = 1, sigma_prima*dt = (1-1)*dt = 0*dt = 0
% HEAD CELLS -> beta = 0, sigma_prima*dt = (1-0)*dt = 1*dt = dt 
% [NOTE] - We are considering excitation cell as a special wall cell
typeValues = zeros(2, cell_numTypes);
typeValues(:, cell_wall+1) = [0, dt];
typeValues(:, cell_air+1) = [1, 0];         % air
typeValues(:, cell_noPressure+1) = [1, 0];  % air
typeValues(:, cell_excitation+1) = [0, dt]; % excitation
typeValues(:, cell_head+1) = [0, dt];       % head cells

% Define sigma for PMLLayers 
for pmlCounter = 0:pmlLayer-1
    sigmadt(pmlCounter+1) = (pmlCounter/(pmlLayer-1)) * maxSigmadt;
    typeValues(:, cell_pml0+1+pmlCounter) = [1, sigmadt(pmlCounter+1)];
end
typeValues(:, cell_dead+1)    =  [0, 1000000];  % dead cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION TYPES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simulationType = input('Choose Simulation Type[0-Open Space 1-TubeWall 2-VerticalWall 3-BothEnd OpenTube 4-VowelSound 5-Measure RadiationImpedance]: ');
pmlSwitch = input('Swicth ON PML Layers. Press 1:ON 0:OFF = ');

switch simulationType
    case 0 % For open air space simulation
        % Fix vowelSound to zero for open air-space simulation
        vowelSound = 0;
        
        % Create the frame
        [listenerX, listenerY, glottalDiameter2D, glottalDiameter3D, upperTubeEndY, lowerTubeEndY, tubeEndX, frameH, frameW, mu2D, boundarySegmentType, depthX, depthY, depthP, baffleSwitch, micXposCells, glottalArea, PV_N]= ...
             vt_VowelTubeGeneration(pmlSwitch, ds, vowelSound,...
             simulation2D, cell_wall, cell_air, cell_excitation, cell_noPressure, cell_head); 
                
        % Retrieve excitationV position
        [exeY, exeX] = find(PV_N(:,:,4)==cell_excitation);
    
    case 1 % For fixed size tube-wall simulation
        
        % domainW an domainH signifies the size of problem space
        domainW = input('Enter domain width: ');
        domainH = input('Enter domain height: ');
        
        % Create the frame
        [PV_N, frameH, frameW, depthX, depthY, depthP] = ...
         vt_UniformTubeGeneration(domainH, domainW, pmlSwitch, pmlLayer, simulation2D);
        
        % Define cell type and store it in PV_N(,,4)
        % Declare all the cells as air by default
        PV_N(1:frameH, 1:frameW, 4) = cell_air;
        
        % Define source position
        excitationX = floor(frameW/2);  
        excitationY = floor(frameH/2);
        
        % Define source size
        excitationH = 5;
        excitationW = 1;
        
        % Defining source cell type
        cellType = cell_excitation;
        PV_N(excitationY+(0:excitationH-1), excitationX+(0:excitationW-1), 4) = cellType;
        
        % Check tube length
        tubeLength = input('Enter tube length: ');

        % Fix listener position
        listenerX = excitationX + tubeLength-1;
        listenerY = excitationY;
        
        % To implement Dirichlet Boundary Condition define cell_type
        for i=0:excitationH+1
                PV_N(listenerY-1+i, listenerX+1, 4) = cell_noPressure;
        end
        
        %back walls
        for i=0:excitationH+1
            PV_N(excitationY-1+i, excitationX-1, 4) = cell_wall;
        end

        %tube walls
        for j=excitationX-1:listenerX
            PV_N(excitationY-1, j, 4) = cell_wall;
            PV_N(excitationY+excitationH, j, 4) = cell_wall;
        end
        
    case 2 % For vertical wall simulation
        
        % domainW an domainH signifies the size of problem space
        domainW = input('Enter domain width: ');
        domainH = input('Enter domain height: ');
        
        % Create the frame
        [PV_N, frameH, frameW, depthX, depthY, depthP] = ...
         vt_UniformTubeGeneration(domainH, domainW, pmlSwitch, pmlLayer, simulation2D);
        
        % Define cell type and store it in PV_N(,,4)
        % Declare all the cells as air by default
        PV_N(1:frameH, 1:frameW, 4) = cell_air;
        
        % Define source position
        excitationX = floor(frameW/2);  
        excitationY = floor(frameH/2);
        
        % Define source size
        excitationH = 1;
        excitationW = 1;
        
        % Define Listener postion
        listenerX = excitationX;
        listenerY = excitationY;
        
        % Defining source cell type
        cellType = cell_excitation;
        PV_N(excitationY+(0:excitationH-1), excitationX+(0:excitationW-1), 4) = cellType;
        PV_N(excitationY-1:excitationY+20, excitationX+20, 4) = cell_wall;
        
    case 3 % Both end open tube
        
        % domainW an domainH signifies the size of problem space
        domainW = input('Enter domain width: ');
        domainH = input('Enter domain height: ');
        
        % Create the frame
        [PV_N, frameH, frameW] = vt_UniformTubeGeneration(domainH, domainW, pmlSwitch, pmlLayer);
        
        % Define cell type and store it in PV_N(,,4)
        % Declare all the cells as air by default
        PV_N(1:frameH, 1:frameW, 4) = cell_air;
        
        % Define source position
        excitationX = floor(frameW/2);  
        excitationY = floor(frameH/2);
        
        % Defining source cell type
        cellType = cell_excitation;
        PV_N(excitationY+(0:excitationH-1), excitationX+(0:excitationW-1), 4) = cellType;
        
        % Check tube length
        tubeLength = input('Enter tube length: ');
        
        % Fix listener position
        listenerX = excitationX + tubeLength-1;
        listenerY = excitationY;
        
        % To implement Dirichlet Boundary Condition define cell_type
        for i=0:excitationH+1
                PV_N(listenerY-1+i, listenerX+1, 4) = cell_noPressure;
                % we don't need cell_noPressure behind the excitation wall,
                % if the excitation is only along the forward direction.
                %PV_N(excitationY-1+i, excitationX-1, 4) = cell_noPressure;
        end
        
        %tube walls
        for j=excitationX:listenerX
            PV_N(excitationY-1, j, 4) = cell_wall;
            PV_N(excitationY+excitationH, j, 4) = cell_wall;
        end
        
    case 4 % For vowel sound
        vowelSound = input('Choose vowels [1-\a\ 2-\u\ 3-\o\ 4-\i\]: ');
        
        % Keep asking user enter wrong input
        while vowelSound<1 || vowelSound>4
            disp('Give correct input')
            vowelSound = input('Choose vowels [1-\a\ 2-\u\ 3-\o\ 4-\i\]: ');
        end
              
        % Generate the Tube Shape
        [listenerX, listenerY, glottalDiameter2D, glottalDiameter3D, upperTubeEndY, lowerTubeEndY, tubeEndX, frameH, frameW, mu2D, boundarySegmentType, depthX, depthY, depthP, baffleSwitch, micXposCells, glottalArea, PV_N]= ...
             vt_VowelTubeGeneration(pmlSwitch, ds, vowelSound, ...
             simulation2D, cell_wall, cell_air, cell_excitation, cell_noPressure, cell_head);   
         
        % Retrieve excitationV position
        [exeY, exeX] = find(PV_N(:,:,4)==cell_excitation);
        
     case 5 % For vowel sound - Impedance Measurement
        vowelSound = 4;
        % Generate the Tube Shape
        [listenerX, listenerY, glottalDiameter2D, glottalDiameter3D, upperTubeEndY, lowerTubeEndY, tubeEndX, frameH, frameW, mu2D, boundarySegmentType, depthX, depthY, depthP, baffleSwitch, micXposCells, glottalArea, PV_N]= ...
             vt_VowelTubeGeneration(pmlSwitch, ds, vowelSound, ...
             simulation2D, cell_wall, cell_air, cell_excitation, cell_noPressure, cell_head);   
         
        % Retrieve excitationV position
        [exeY, exeX] = find(PV_N(:,:,4)==cell_excitation);
    otherwise
end

% Calculate the Z inverse
% Compute the impedance for the whole frame
z_inv = zeros(frameH, frameW);

for row_idx = 1:frameH
    for col_idx = 1:frameW
        if boundarySegmentType(row_idx, col_idx) ~=0
            mu2DCurr = mu2D(boundarySegmentType(row_idx, col_idx));
            alphaCurr = 1/(0.5+0.25*(mu2DCurr +(1/mu2DCurr)));
            z_inv(row_idx, col_idx) = 1 / (rho*c*( (1+sqrt(1-alphaCurr))/(1-sqrt(1-alphaCurr)) ));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOURCE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sourceModelType = input('Choose the source model type [1-Sine Wave 2-Gaussian Source 3-Impulse Function 4-Vocal Fold Model]: ');

switch sourceModelType
    case 1 % Sine wave source model
        excitationF = 440*hertz;            
        srcAmplitude =25;
        exeT = linspace(1, STEPS, STEPS);
        excitationV = srcAmplitude * sin(2*pi*excitationF*dt*(exeT(:)-1));
        
    case 2 % Gaussian source model
        f0 = 10*kilohertz;
        bellPeakPos = 0.646/f0;
        bellWidth = 0.29*bellPeakPos;
        excitationV3D = exp(-((t-bellPeakPos)./bellWidth).^2);
        excitationV = excitationV3D.*((glottalDiameter3D/glottalDiameter2D)^2);
        
    case 3 % Impulse response function
        excitationV = src_ImpulseSignal(srate, 10000, 2, 22050*hertz);
        
    case 4 % Vocal fold model - Two Mass Model
        % Add your code here - Set the vocal fold parameters
        [airParam, vf_structuralParam, vf_flowParam, vf_matParam] = vf_SetVocalFoldParams();
    otherwise
end
 
% Add zeros to the source signal if the number of steps are greater than
% the source signal length
if sourceModelType~=4 && STEPS>length(excitationV)
    excitationV = [excitationV; zeros(STEPS-length(excitationV),1)];
end

% Define source propagation direction
% srcDirection index mean: 1 = Left  = -1
%                          2 = Down  = -1
%                          3 = Right =  1
%                          4 = Up    =  1

srcDirection = [0 0 1 0]; % For all the direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MIC Position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 0: To store the data from first mic position [Near to tube end]
Pr_Audio = zeros(1, STEPS);
Vx_Vel = zeros(1, STEPS);
Vy_Vel = zeros(1, STEPS);

% Case 1: To store the data from second mic position [Near to PML Layer]
Pr_Audio1 = zeros(6, STEPS);
Vx_Vel1 = zeros(6, STEPS);
Vy_Vel1 = zeros(6, STEPS);

% case 2: To store the data from the second mic [Near to excitation]
Pr_Audio2 = zeros(1, STEPS);
Vx_Vel2 = zeros(1, STEPS);
Vy_Vel2 = zeros(1, STEPS);

% Case 3: To store the data from the third mic [In between excitation and listener]
Pr_Audio3 = zeros(1, STEPS);
Vx_Vel3 = zeros(1, STEPS);
Vy_Vel3 = zeros(1, STEPS);

% Case 4: To store the data from the fourth mic [At the excitation]
Pr_Audio4 = zeros(1, STEPS);
Vx_Vel4 = zeros(1, STEPS);
Vy_Vel4 = zeros(1, STEPS);

% Case 5: To store the data in an array of mic [TubeEnd] to average out
Pr_Audio5 = zeros((lowerTubeEndY-2)-(upperTubeEndY+2)+1, STEPS);
Vx_Vel5 = zeros((lowerTubeEndY-2)-(upperTubeEndY+2)+1, STEPS);
Vy_Vel5 = zeros((lowerTubeEndY-2)-(upperTubeEndY+2)+1, STEPS);

% Case 6: To store the data in an array of mic [TubeEnd+1] to average out
Pr_Audio6 = zeros((lowerTubeEndY-2)-(upperTubeEndY+2)+1, STEPS);
Vx_Vel6   = zeros((lowerTubeEndY-2)-(upperTubeEndY+2)+1, STEPS);
Vy_Vel6   = zeros((lowerTubeEndY-2)-(upperTubeEndY+2)+1, STEPS);

% Find the upper and lower tube wall at the listenerX position
checkWallX = tubeEndX-micXposCells;
checkWallY = listenerY;
radiiCounter=0;

while PV_N(checkWallY, checkWallX, 4)~= cell_wall
    radiiCounter = radiiCounter+1;
    checkWallY = listenerY-radiiCounter;
end

% Case 7: To store the data in an array of mic [TubeEnd+1] to average out
% [Inside the tube]

Pr_Audio7   = zeros((listenerY+radiiCounter-1)-(listenerY-radiiCounter+1)+1, STEPS);
Vx_Vel7     = zeros((listenerY+radiiCounter-1)-(listenerY-radiiCounter+1)+1, STEPS);
Vy_Vel7     = zeros((listenerY+radiiCounter-1)-(listenerY-radiiCounter+1)+1, STEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define PML Layer Cells and Dead Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Then modify as per the requirement
% Define cell_dead to the outer most layer
PV_N(1:frameH,1,4) = cell_dead;
PV_N(1:frameH,frameW,4) = cell_dead;
PV_N(1,1:frameW,4) = cell_dead;
PV_N(frameH,1:frameW,4) = cell_dead;

if pmlSwitch == 1 
    
    % Define horizontal PML layers - Start from the outer layers
    % -----Activate horizontal PML Layers-------
    cellType = cell_pml5;
    yShift = 1;
    xStart = 2;
    xEnd = frameW-1;

    for pmlCount = 1:pmlLayer   
        for hCount = xStart:xEnd
            PV_N(yShift+pmlCount, hCount, 4) = cellType;
            PV_N(frameH-pmlCount, hCount, 4) = cellType;
        end   
        xStart = xStart+1;
        xEnd = xEnd-1;
        cellType = cellType-1;
    end

    % -----Activate vertical PML Layers-------
    cellType = cell_pml5;
    xShift = 1;
    yStart = 2;
    yEnd = frameH-1;
    for pmlCount = 1:pmlLayer
        for vCount = yStart:yEnd
            PV_N(vCount, xShift+pmlCount, 4) = cellType;
            PV_N(vCount, frameW-pmlCount, 4) = cellType;
        end

        yStart = yStart+1;
        yEnd = yEnd-1;
        cellType = cellType-1;
    end
end

% Test the frame
frame = PV_N(:, :,4);
frame(listenerY, listenerX) = -1;
frame(listenerY, listenerX+5) = -1;
imagesc(frame);
title('Frame [Domain+PML]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE MINBETA AND MAXSIGMA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Storing min(Beta) and max(sigmaPrimedt) - We'll use this to calculate velocity
height = frameH-2;
width  = frameW-2;

cellTypes = zeros(height*width, 3);
typeIndex = zeros(height*width, 3);
beta = zeros(height*width, 3);
sigma_prime_dt = zeros(height*width, 3);

% Define a counter
counter = 1;

for row_idx = 2:frameH-1
    for col_idx = 2: frameW-1
        
        % Find the cellTypes
         cellTypes(counter,:) = ...
             [PV_N(row_idx, col_idx, 4), PV_N(row_idx, col_idx+1, 4), PV_N(row_idx-1, col_idx, 4)];
         
         % For typeIndex add 1 to cellTypes
         typeIndex(counter,:) = cellTypes(counter,:)+1;
         
         % Compute beta
         beta(counter,:) = typeValues(1, typeIndex(counter,:));  
         
         % Computer sigma_prime_dt
         sigma_prime_dt(counter,:) = typeValues(2, typeIndex(counter,:));
                 
         % Increase the counter
         counter=counter+1;
    end   
end

min_beta_Vx = min(beta(:,[1,2]),[],2);
min_beta_Vy = min(beta(:,[1,3]),[],2);

max_sigma_prime_dt_Vx = max(sigma_prime_dt(:,[1,2]),[],2);
max_sigma_prime_dt_Vy = max(sigma_prime_dt(:,[1,3]),[],2);

% For storing minBeta and maxSigmaPrimedt
minVxBeta = reshape(min_beta_Vx,[width,height])';
minVyBeta = reshape(min_beta_Vy,[width,height])';

maxVxSigmaPrimedt = reshape(max_sigma_prime_dt_Vx,[width,height])';
maxVySigmaPrimedt = reshape(max_sigma_prime_dt_Vy,[width,height])';

% For storing PressureSigmaPrimedt
PressureSigmaPrimedt = reshape( typeValues(2,typeIndex(:,1)),[width,height])';

betaVxSqr = minVxBeta.*minVxBeta;
betaVxSqr_dt_invRho = (betaVxSqr.*dt)/rho;

betaVySqr = minVyBeta.*minVyBeta;
betaVySqr_dt_invRho = (betaVySqr.*dt)/rho;

rho_sqrC_dt_invds = (kappa*dt)/dx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDTD PARAMETERS - PART II (COMPUTING ALL THE REQUIRED DATA TO INJECT SOURCE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify if the current cell is a part of circular baffle
% z_invCurr = z_inv.*(cellTypes(:,1)~=cell_head);

% Check whether the current cell is an excitation cell 
is_excitation = [cellTypes(:,1) == cell_excitation, cellTypes(:, 2) == cell_excitation, cellTypes(:, 3) == cell_excitation];

% Verify if we are not the excitation cell
are_we_not_excitations = [ (1 - is_excitation(:,1)) .* (1 - is_excitation(:,2)),...
                           (1 - is_excitation(:,1)) .* (1 - is_excitation(:,3))];
                       
% Compute the excitation weight
excitation_weight = [is_excitation(:,1) is_excitation(:,1)].*srcDirection(3:4) + [is_excitation(:,2) is_excitation(:,3)].*srcDirection(1:2);

is_normal_dir = [beta(:,2) ~= cell_air, beta(:,3) ~= cell_air, beta(:,3) == cell_air, beta(:,2) == cell_air];
xor_term = [beta(:,2) .* (1-beta(:,1)) , beta(:,1) .* (1-beta(:,2)), ...
            beta(:,3) .* (1-beta(:,1)), beta(:,1) .* (1-beta(:,3))];
        
N = [0.707106*is_normal_dir(:,3) + (1-is_normal_dir(:,3)), 0.707106*is_normal_dir(:,2) + (1-is_normal_dir(:,2)), ...
     0.707106*is_normal_dir(:,4) + (1-is_normal_dir(:,4)), 0.707106*is_normal_dir(:,1) + (1-is_normal_dir(:,1))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDTD PARAMETERS - PART III (INITIALISING FIELD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV_Nplus1  = zeros(frameH, frameW, 4);
Vx_test    = zeros(frameH, frameW);
Vy_test    = zeros(frameH, frameW);


audio_Vis = zeros(frameH, frameW); % Dispaly this array during simmulation

CxVx    = zeros(frameH-2,frameW-2);
CyVy    = zeros(frameH-2,frameW-2);
CxP     = zeros(frameH-2,frameW-2);
CyP     = zeros(frameH-2,frameW-2);
Pr_next = zeros(frameH-2,frameW-2);
Vx_next = zeros(frameH-2,frameW-2);
Vy_next = zeros(frameH-2,frameW-2);

% Open a new figure window to visualize the simulation
figure;


% Figure out the row and column number from the middle cell of  the row of
% excitation cells
[rowNums, colNums] = find(PV_N(:,:,4)==cell_excitation);
rowExcitation = rowNums (((length(rowNums)-1)/2) +1);
columnExcitation = colNums (((length(colNums)-1)/2) +1);

for T = 1:STEPS
    
    % STEP1: Calculate (del.V) = (dVx/dx + dVy/dy)
    % CxVx = dVx/dx, where Vx = velocity along x direction
    % CyVy = dVy/dy, where Vy = velocity along y direction
    
    CxVx(1:frameH-2, 1:frameW-2)= (PV_N(2:frameH-1, 2:frameW-1,2).*depthX(2:frameH-1, 2:frameW-1)) - ...
                                  (PV_N(2:frameH-1, 1:frameW-2,2).*depthX(2:frameH-1, 1:frameW-2));
                              
    CyVy(1:frameH-2, 1:frameW-2)= (PV_N(2:frameH-1, 2:frameW-1, 3).*depthY(2:frameH-1, 2:frameW-1)) - ...
                                  (PV_N(3:frameH, 2:frameW-1, 3).*depthY(3:frameH, 2:frameW-1));
                                  
    % STEP2: Calculate Pr_next                           
    Pr_next(1:frameH-2, 1:frameW-2) = ((PV_N(2:frameH-1, 2:frameW-1, 1).*depthP(2:frameH-1, 2:frameW-1)) - ((rho_sqrC_dt_invds.*(CxVx+CyVy))))./...
                                     ((1+PressureSigmaPrimedt).*depthP(2:frameH-1, 2:frameW-1));

    % STEP3: Copy Pr_next  to PV_Nplus1
    PV_Nplus1(2:frameH-1, 2:frameW-1,1) = Pr_next(:,:);
        
    % STEP4: Implement Dirichlet Boundary Condition
    if baffleSwitch~=1
        checkPresCond = PV_N(:,:,4)~=cell_noPressure;
        PV_Nplus1(:,:,1) = PV_Nplus1(:,:,1).*checkPresCond;
    end
         
    % STEP5: Calculate Vx & Vy
    % To compute Vx we need calculate CxP = (del.P) = dPx/dx
    % To compute Vy we need calculate CyP = (del.P) = dPy/dy
    
    CxP(1:frameH-2, 1:frameW-2) = (PV_Nplus1(2:frameH-1,3:frameW,1) - PV_Nplus1(2:frameH-1,2:frameW-1,1))/dx;
    CyP(1:frameH-2, 1:frameW-2) = (PV_Nplus1(1:frameH-2,2:frameW-1,1) - PV_Nplus1(2:frameH-1,2:frameW-1,1))/dy;
    
    Vx_next(1:frameH-2, 1:frameW-2) = (minVxBeta.*PV_N(2:frameH-1,2:frameW-1,2)- (betaVxSqr_dt_invRho.*CxP));                                
    Vy_next(1:frameH-2, 1:frameW-2) = (minVyBeta.*PV_N(2:frameH-1,2:frameW-1,3)- (betaVySqr_dt_invRho.*CyP));
    
    PV_Nplus1(2:frameH-1, 2:frameW-1,2) = Vx_next(:,:);
    PV_Nplus1(2:frameH-1, 2:frameW-1,3) = Vy_next(:,:);
    
    Vx_test = PV_Nplus1(:,:,2);
    Vy_test = PV_Nplus1(:,:,3);
    
    % Reset the counter
    counter = 1;
    
    % For Vocal Fold Model
    if sourceModelType==4
        % Retrive supra glottal presure from the middle cell of the row of
        % excitation cells
        vf_flowParam.p1 = PV_Nplus1(rowExcitation, columnExcitation+1);
        
        % Compute the updated vocalfold parameters
        [vf_flowParam] = vf_TwoMassModel(srate, airParam, vf_structuralParam, vf_flowParam);
    end
    
            
    for row_idx = 2:frameH-1
        for col_idx = 2: frameW-1
            
            % Define the excitation velocity based on the source model
            if sourceModelType==4
                exeV = vf_flowParam.ug_curr/glottalArea;
            else
                exeV = excitationV(T);
            end
          
            % Inject the source to the Vx_next and Vy_next = excitationV(T)
            PV_Nplus1(row_idx, col_idx, 2) = PV_Nplus1(row_idx, col_idx, 2) + exeV*excitation_weight(counter,1)*maxVxSigmaPrimedt(row_idx-1, col_idx-1);
            PV_Nplus1(row_idx, col_idx, 3) = PV_Nplus1(row_idx, col_idx, 3) + exeV*excitation_weight(counter,2)*maxVySigmaPrimedt(row_idx-1, col_idx-1);
            
            % Compute vb_alpha
            vb_alpha = [xor_term(counter,2)*PV_Nplus1(row_idx,col_idx,1)*N(counter,2) - xor_term(counter,1)*PV_Nplus1(row_idx,col_idx+1,1)*N(counter,1), ...
                        xor_term(counter,4)*PV_Nplus1(row_idx,col_idx,1)*N(counter,4) - xor_term(counter,3)*PV_Nplus1(row_idx-1,col_idx,1)*N(counter,3)];
            
            vb_alpha = (vb_alpha.* are_we_not_excitations(counter,:))*z_inv(row_idx, col_idx);
            
            % Update Vx and Vy
            PV_Nplus1(row_idx,col_idx,2) = PV_Nplus1(row_idx,col_idx,2) + maxVxSigmaPrimedt(row_idx-1,col_idx-1) * vb_alpha(1);
            PV_Nplus1(row_idx,col_idx,3) = PV_Nplus1(row_idx,col_idx,3) + maxVySigmaPrimedt(row_idx-1,col_idx-1) * vb_alpha(2);
            
            % Increment the counter
            counter = counter+1;
        end
    end
       
    PV_Nplus1(2:frameH-1, 2:frameW-1,2) = PV_Nplus1(2:frameH-1, 2:frameW-1,2)./(minVxBeta+maxVxSigmaPrimedt);
    PV_Nplus1(2:frameH-1, 2:frameW-1,3) = PV_Nplus1(2:frameH-1, 2:frameW-1,3)./(minVyBeta+maxVySigmaPrimedt);
    
    % STEP8: Re-store the grid cell type
    PV_Nplus1(2:frameH-1, 2:frameW-1,4) = PV_N(2:frameH-1, 2:frameW-1,4);
    
    %STEP9: Clear the border cell
    
    PV_Nplus1(:, 1, 1:3) = 0;
    PV_Nplus1(:, 1, 4) = PV_N(:, 1, 4);
    
    PV_Nplus1(:, frameW, 1:3) = 0;
    PV_Nplus1(:, frameW, 4) = PV_N(:, frameW, 4); 

    PV_Nplus1(1, :, 1:3) = 0;
    PV_Nplus1(1, :, 4) = PV_N(1, :, 4);
    
    PV_Nplus1(frameH, :, 1:3) = 0;
    PV_Nplus1(frameH, :, 4) = PV_N(frameH, :, 4); 
    
    audio_Vis = PV_Nplus1(:,:,1);
    % To visualize the boundary
    audio_Vis(PV_Nplus1(:,:,4)==cell_wall) = vis_Boundary; 
    audio_Vis(PV_Nplus1(:,:,4)==cell_head) = vis_Boundary;
    
    % STEP10: Plot wave simulation
    if ~mod(T,1)
        imagesc(audio_Vis,[-1000 4000]);  %colorbar; % Multiplied with twenty to change the color code
        %xlabel('Spatial Resolution along X');
        %ylabel('Spatial Resolution along Y');
        title('Vowel Simulation');
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        %title(['STEP NUMBER: ' num2str(T) ' OUT OF ' num2str(STEPS)]);
        drawnow;     
    end
    
    % STEP11: Copy PV_Nplus1 to PV_N for the next time step
    PV_N = PV_Nplus1;
    
    
   % Case 0: Store data from the first mic position
    Pr_Audio(T) = PV_Nplus1(listenerY, listenerX,1);
    Vx_Vel(T)   = PV_Nplus1(listenerY, listenerX,2);
    Vy_Vel(T)   = PV_Nplus1(listenerY, listenerX,3);
    
    % Case 1: Store data from the second mic position [Near to PML Layer]
%     Pr_Audio1(1,T) = PV_Nplus1(listenerY, listenerX-1,1);
%     Vx_Vel1(1,T)   = PV_Nplus1(listenerY, listenerX-1,2);
%     Vy_Vel1(1,T)   = PV_Nplus1(listenerY, listenerX-1,3);
%     
%     Pr_Audio1(2,T) = PV_Nplus1(listenerY, listenerX+1,1);
%     Vx_Vel1(2,T)   = PV_Nplus1(listenerY, listenerX+1,2);
%     Vy_Vel1(2,T)   = PV_Nplus1(listenerY, listenerX+1,3);
%     
%     Pr_Audio1(3,T) = PV_Nplus1(listenerY, listenerX+2,1);
%     Vx_Vel1(3,T)   = PV_Nplus1(listenerY, listenerX+2,2);
%     Vy_Vel1(3,T)   = PV_Nplus1(listenerY, listenerX+2,3);
%     
%     Pr_Audio1(4,T) = PV_Nplus1(listenerY, listenerX+3,1);
%     Vx_Vel1(4,T)   = PV_Nplus1(listenerY, listenerX+3,2);
%     Vy_Vel1(4,T)   = PV_Nplus1(listenerY, listenerX+3,3);
%     
%     Pr_Audio1(5,T) = PV_Nplus1(listenerY, listenerX+4,1);
%     Vx_Vel1(5,T)   = PV_Nplus1(listenerY, listenerX+4,2);
%     Vy_Vel1(5,T)   = PV_Nplus1(listenerY, listenerX+4,3);
%     
%     Pr_Audio1(6,T) = PV_Nplus1(listenerY, listenerX+5,1);
%     Vx_Vel1(6,T)   = PV_Nplus1(listenerY, listenerX+5,2);
%     Vy_Vel1(6,T)   = PV_Nplus1(listenerY, listenerX+5,3);
  
    % Case 2: Store data from the second mic position [Near to excitation]
%     Pr_Audio2(T) = PV_Nplus1(listenerY, exeX(1)+3,1);
%     Vx_Vel2(T)   = PV_Nplus1(listenerY, exeX(1)+3,2);
%     Vy_Vel2(T)   = PV_Nplus1(listenerY, exeX(1)+3,3);
     
    % Case 3: Store data from the third mic position [Mid position of source and listener]
%     Pr_Audio3(T) = PV_Nplus1(listenerY, exeX(1)+round((listenerX-exeX(1))/2),1);
%     Vx_Vel3(T)   = PV_Nplus1(listenerY, exeX(1)+round((listenerX-exeX(1))/2),2);
%     Vy_Vel3(T)   = PV_Nplus1(listenerY, exeX(1)+round((listenerX-exeX(1))/2),3);
    
    % Case 4: Store data from the fourth mic position [At the excitation]
%     Pr_Audio4(T) = PV_Nplus1(listenerY, exeX(1),1);
%     Vx_Vel4(T)   = PV_Nplus1(listenerY, exeX(1),2);
%     Vy_Vel4(T)   = PV_Nplus1(listenerY, exeX(1),3);
    
    % Case 5: STore the data in an array of mic [TubeEnd] to average out
%     Pr_Audio5(:,T) = PV_Nplus1(upperTubeEndY+2:lowerTubeEndY-2, tubeEndX,1);
%     Vx_Vel5(:,T)   = PV_Nplus1(upperTubeEndY+2:lowerTubeEndY-2, tubeEndX,2);
%     Vy_Vel5(:,T)   = PV_Nplus1(upperTubeEndY+2:lowerTubeEndY-2, tubeEndX,3);
    
    % Case 6: STore the data in an array of mic [TubeEnd] to average out
%     Pr_Audio6(:,T) = PV_Nplus1(upperTubeEndY+2:lowerTubeEndY-2, tubeEndX+1,1);
%     Vx_Vel6(:,T)   = PV_Nplus1(upperTubeEndY+2:lowerTubeEndY-2, tubeEndX+1,2);
%     Vy_Vel6(:,T)   = PV_Nplus1(upperTubeEndY+2:lowerTubeEndY-2, tubeEndX+1,3);
    
%     Pr_Audio7(:,T) = PV_Nplus1(listenerY-radiiCounter+1:listenerY+radiiCounter-1, checkWallX, 1);
%     Vx_Vel7(:,T)   = PV_Nplus1(listenerY-radiiCounter+1:listenerY+radiiCounter-1, checkWallX, 2);
%     Vy_Vel7(:,T)   = PV_Nplus1(listenerY-radiiCounter+1:listenerY+radiiCounter-1, checkWallX, 3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE DATA IN "impedanceData.mat" FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the velocity magnitude
% vRes = sqrt(Vx_Vel.^2 + Vy_Vel.^2);
% vRes1 = sqrt(Vx_Vel1.^2 + Vy_Vel1.^2);
% vRes2 = sqrt(Vx_Vel2.^2 + Vy_Vel2.^2);
% vRes3 = sqrt(Vx_Vel3.^2 + Vy_Vel3.^2);
% vRes4 = sqrt(Vx_Vel4.^2 + Vy_Vel4.^2);

% Save the data
% save('impedanceData.mat','excitationV','Pr_Audio','Vx_Vel','Vy_Vel',...
%       'Pr_Audio1','Vx_Vel1','Vy_Vel1','Pr_Audio2','Vx_Vel2','Vy_Vel2',...
%       'Pr_Audio3','Vx_Vel3','Vy_Vel3','Pr_Audio4','Vx_Vel4','Vy_Vel4',...
%       'vRes', 'vRes1','vRes2', 'vRes3','vRes4', 'Pr_Audio5','Vx_Vel5',...
%       'Vy_Vel5','Pr_Audio6','Vx_Vel6','Vy_Vel6');

% save('impedanceData.mat','excitationV','Pr_Audio5','Vx_Vel5',...
%       'Vy_Vel5','Pr_Audio6','Vx_Vel6','Vy_Vel6','Pr_Audio7',....
%       'Vx_Vel7','Vy_Vel7');

% Create the audio data and play
generateAudio(Pr_Audio, srate, srate_mul);

save('vowel_i25.mat', 'Pr_Audio', 'srate', 'srate_mul');