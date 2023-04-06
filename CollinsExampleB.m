clc
clear
close all
AddPathsAndCleanUp()


saveoutputs = 1; %1 to save, anything else to not save
runtitle = 'Collins Example B w corrmaxalg 0.99.mat';
savepath = 'Jons Tools/SavedRuns/';


m = 10;
n = 40;

%Array length
deltax = 7.5*(n-1);
deltaz = 7.5*(m-1);
ArrayCoords = GenerateRectangularArray(m,n,deltax,deltaz);
Arraymat = reshape(ArrayCoords,[m*n,3]);

%% Specify the Array's Orientation

%Euler Angles in radians
yaw = deg2rad(0);
pitch = deg2rad(0);
roll = deg2rad(0);

%% Specify the Array's Location

%Depth
ArrayDepth = 50 + deltaz/2;

%% Specify Source Information

%Collins paper example B
SourceInitial = [ 5,15,0.16;
                  16,-12.5,0.11;
                  -18,6,0.5;
                  1,-5,0.15 ];
SourceFinal = [5,15,0.36;
               8,-8.5,0.11;
               -10,-2,0.5;
               4.5,-0.8,0.15];

Velocity = 1000*(SourceFinal - SourceInitial);

SourceLocation = 1000*SourceInitial;

SourceAmplitude = [0.6,1,0.5,0.75];
SourceOnswitch = [1,1,1,1];
SourceAmplitude = SourceOnswitch.*SourceAmplitude;
numdetects = sum(SourceOnswitch);

dt = 0.01;
time = 0:dt:1;

plotoption = 0;
sigma = 3.6526e-08;

phi = 5:1:180;
theta = 0;

%% Do some conversions

SourceFreq = 100;


extraeigs = 1;

%Eigenvector transformation coefficient (equations 9 and 10 from collins
%paper
computesubarraycorrelation = 1;
ApplyTransformation = 1;
CorrelationThreshold = 0.99;
AlphaRange = 0:0.01:1;
AngleRange = 0:pi/120:2*pi;





RunKrakenMVBSim;  


if saveoutputs == 1
    save([savepath,runtitle]);
end


%%
function [] = AddPathsAndCleanUp()
    addpath('Matlab/')
    addpath('Matlab/ReadWrite')
    addpath('Jons Tools')
    delete KRAKEN_MAT.env KRAKEN_MAT.mod KRAKEN_MAT.prt
end