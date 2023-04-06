%% Preliminaries

clc
clear
close all
AddPathsAndCleanUp;

saveoutputs = 1;
runtitle = '2 Sources (Crossing SNRs, Same Amp, Different Bearings, Asymmetric 3)';
savepath = 'Jons Tools/SavedRuns/';

%% Array Geometry

m = 20;
n = 20;
%Array length
deltax = 7.5*m;
deltaz = 7.5*n;
ArrayCoords = GenerateRectangularArray(m,n,deltax,deltaz);
Arraymat = reshape(ArrayCoords,[m*n,3]);

%% Specify the Array's Orientation

%Euler Angles in radians
yaw = deg2rad(0);
pitch = deg2rad(0);
roll = deg2rad(0);


%% Specify the Array's Location

%Range From the Source
ArrayRange = 0;

%Depth
ArrayDepth = 50 + deltaz/2;

%% Specify Source Information

%Two sources moving in opposite directions through broadside
SourceInitial = [-1.4,6,ArrayDepth/1000;
                  1.4,12,ArrayDepth/1000];
SourceFinal = [-1.4,12,ArrayDepth/1000;
               1.4,6,ArrayDepth/1000];

Velocity = 1000*(SourceFinal - SourceInitial);

SourceLocation = 1000*SourceInitial;

SourceAmplitude = [1,1];
SourceOnswitch = [1,1];
SourceAmplitude = SourceOnswitch.*SourceAmplitude;
numdetects = sum(SourceOnswitch);

dt = 0.01;
time = 0:dt:1;

%% Do some conversions

SourceFreq = 100;
HorizSubIndex = 1;
VertSubIndex = 1;

extraeigs = 4;

%Grab the postion vectors corresponding to the subarrays
HorizIndex = (HorizSubIndex-1)*m + 1:HorizSubIndex*m;
VertIndex = VertSubIndex:m:n*m;

%Eigenvector transformation coefficient
computesubarraycorrelation = 1;
ApplyTransformation = 0;
Alpha = 1;
angle = pi/2;


%% Get everything Set Up the Way Kraken Wants it

% If you want to plot the array's orientation
plotoption = 0;
sigma = 3.6526e-08;

phi = 70:0.1:110;
theta = 0;

RunKrakenMVBSim;  


if saveoutputs == 1
    save([savepath,runtitle]);
end

%% Function Dump


function [] = AddPathsAndCleanUp()
    addpath('Matlab/')
    addpath('Matlab/ReadWrite')
    addpath('Jons Tools')
    delete KRAKEN_MAT.env KRAKEN_MAT.mod KRAKEN_MAT.prt
end