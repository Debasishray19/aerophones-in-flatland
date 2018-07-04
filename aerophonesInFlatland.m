%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDTD 2D MODELLING FOR ACOUSTIC WAVE ANALYSIS INSIDE A TUBE
%
% This code implments a cylinderical tube model inside a grid to analyze 
% acoustic wave propagation using FDTD engine. To understand the technical 
% details and conceptual theory behind this implementation, follow these papers - 
% [1] "Aerophones in Flatland: Interactive wave simulation of wind instruments"
% by Andrew allen and Nikunj Raghubansi.
% [2] "Towards a real-time two-dimensional wave propagation for articulatory 
% speech synthesis." by Victor Zappi, Arvind and Sidney Fels
% [3] "Acoustic Analysis of the vocal tract during vowel production" by 
% finite-difference time-domain method by Hironori Takemoto and Parham Mokhtari.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wave equation for non-PML (FOR 2D):
% dp/dt = -(rho*c*c*r*del.(v)) % del = del operator
% dv/dt = (-1/rho)*(del(P))
% du/dt = (-1/rho)*(del(P))

% Wave equation for PML (FOR 2D):
% dp/dt + sig*p = -(rho*c*c*r*del.(v)) % sig = stretching field
% dv/dt + sig*v = (-1/rho)*(del(P))
% du/dt + sig*v = (-1/rho)*(del(P))

% Wave equation to implement Tube with PML
% dp/dt + sigPrime*p = -(rho*c*c*r*del.(v))
% dv/dt + sigPrime*v = (-beta^2/rho)*(del(P)) + sigprime*Vb
% sigPrime = 1-beta+sigma
%***********************************************************************

% Initialize MATLAB environment
close all;
clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE UNITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meter = 1;
centimeter  =1e-2;

second    = 1;
hertz     = 1/second;
kilohertz = 1e3 * hertz;
megahertz = 1e6 * hertz;
gigahertz = 1e9 * hertz;

gram      = 1e-3;
kilogram  = 1e3*gram;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = 1.1760*kilogram/(meter^3);   % Air density  [kg/m^3]
c   = 343*meter/second;            % Sound speed in air [m/s]
maxSigmaVal = 0.5;                 % Attenuation coefficient at the PML layer
alpha = 0.004;                     % Reflection coefficient
srate = 44100;                     % Sample frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1/srate;                      % Temporal resolution/ sample time period
dx = dt*c*sqrt( 2.0 );             % Spatial resolution along x-direction: CFL Condition
dy = dt*c*sqrt( 2.0 );             % Spatial resolution along x-direction: CFL Condition
AudioTime = 2*second;              % Total audio signal time
kappa = rho*c*c;                   % Bulk modulus
Zn = ((1+sqrt(1-alpha))/(1-sqrt(1-alpha)))*rho*c; % Acoustic Impedance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRID CELL CONSTRUCTION : DEFINE SIGMA AND BETA VALUE FOR EACH CELL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domainW an domainH signifies the size of problem space
domainW = input('Enter domain width: ');
domainH = input('Enter domain height: ');

% tubeHorizontalLength an tubeVerticalLength signifies the size of tube
tubeHorizontalLength = input('Enter horizontal tube length: ');
tubeVerticalLength   = input('Enter vertical tube length: ');

% tubeWidth signifies the distance between tube walls : Right now we are
% assuming it's uniform across the tube.
tubeWidth = input('Enter tube width: ');

% Number of PML layers
pmlLayer = input('Number of PML layer: ');

% Assign sigma value to grtid cells
refFrameSigma = buildFrameSigma(domainW, domainH, pmlLayer, maxSigmaVal, dt);

% Assign beta value to grid cells
[refFrameBeta, Xsrc, Ysrc, Xlis, Ylis] = buildFrameBeta(domainW, domainH, tubeHorizontalLength,...
                   tubeVerticalLength, tubeWidth, pmlLayer);
               
% Calculate sigmaPrime
refSigmaPrime = 1- refFrameBeta + refFrameSigma;

% Validate tube structure
figure('color','w'); imagesc('refFrameSigma');

% Define Grid/Frame size
[Nx, Ny] = size(refFrameSigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOURCE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq = 1*kilohertz;             % source frequency
t = 0:dt:AudioTime-dt;          % time steps
STEPS = length(t);              % Total time steps
Esrc = 25*sin(2*pi*freq*t);     % sinusoidal source wave

% Source parameters for Gaussian source
tau = 0.5/freq;
t0 = 6*tau;
% Esrc = exp(-(((t-t0)./tau).^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRESSURE CHANGE - AUDIO GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pr_Audio = zeros(1, STEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDTD PARAMETERS - PART I (DEFINE UPDATE COEFFICIENT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mUx0 = (refFrameBeta./dt) + (refSigmaPrime./2);
mUx1 = ((refFrameBeta./dt) - (refSigmaPrime./2))./mUx0;
mUx2 = (-(refFrameBeta.*refFrameBeta)./rho)./mUx0;
mUx3 = refSigmaPrime./mUx0;

mVy0 = (refFrameBeta./dt) + (refSigmaPrime./2);
mVy1 = ((refFrameBeta./dt) - (refSigmaPrime./2))./mUx0;
mVy2 = (-(refFrameBeta.*refFrameBeta)./rho)./mUx0;
mVy3 = refSigmaPrime./mUx0;

mPr0 = 1/dt + refSigmaPrime./2;
mPr1 = (1/dt - refSigmaPrime./2)./mPr0;
mPr2 = -kappa./mPr0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDTD PARAMETERS - PART II (INITIALISING FIELD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CxP = zeros(Nx, Ny);
CyP = zeros(Nx, Ny);
CxU = zeros(Nx, Ny);
CyV = zeros(Nx, Ny);
Ux  = zeros(Nx, Ny);
Vy  = zeros(Nx, Ny);
Ubx = zeros(Nx, Ny);
Vby = zeros(Nx, Ny);
Pr  = zeros(Nx, Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDTD ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for T = 1: STEPS
    
    % STEP1 : CxU & CyV   
    CxU(1,:)  = (Ux(1,:) - Ux(Nx,:))/dx;
    CxU(2:Nx,:)  = (Ux(2:Nx,:) - Ux(1:Nx-1,:))/dx;    
    
    CyV(:,1) = (Vy(:,1) - Vy(:,Ny))/dy;
    CyV(:,2:Ny) = (Vy(:,2:Ny) - Vy(:,1:Ny-1))/dy;   
    
    % STEP2 : Solve Pr
    Pr = Pr.*mPr1 + (mPr2.*(CxU+CyV));
    
    % STEP3 : Solve CxP & CyP
    CxP(1:Nx-1,:) = (Pr(2:Nx,:)-Pr(1:Nx-1,:))/dx;
    CxP(Nx,:) = (Pr(1,:)-Pr(Nx,:))/dx;
    
    CyP(:,1:Ny-1) = (Pr(:,2:Ny) - Pr(:,1:Ny-1))/dy;
    CyP(:,Ny) = (Pr(:,1) - Pr(:,Ny))/dy;
    
    % STEP4 : Solve Ubx and Vby
    Ubx(domainW:Nx-domainW, domainH:Ny-domainH) = Pr(domainW:Nx-domainW, domainH:Ny-domainH)/Zn;
    Vby(domainW:Nx-domainW, domainH:Ny-domainH) = Pr(domainW:Nx-domainW, domainH:Ny-domainH)/Zn;
    
    % STEP5 : Solve Ux & Vy
    Ux = mUx1.*Ux + mUx2.*CxP + mUx3.*Ubx;
    Vy = mVy1.*Vy + mVy2.*CyP + mVy3.*Vby;
    
    % STEP6: Inject source
    Pr(Xsrc,Ysrc) = Pr(Xsrc,Ysrc) + Esrc(T);
    
    % STEP7: Store pressure change for audio generation
    Pr_Audio(T) = Pr(Xlis, Ylis);
    
    % STEP8 : Draw the graphics
    if ~mod(T,500)
        imagesc(Pr'*50, [-1,1]); colorbar; % Multiplied with twenty to change the color code
        xlabel('Spatial Resolution along X');
        ylabel('Spatial Resolution along Y');
        title(['STEP NUMBER: ' num2str(T) ' OUT OF ' num2str(STEPS)]);
        drawnow;
    end   
end