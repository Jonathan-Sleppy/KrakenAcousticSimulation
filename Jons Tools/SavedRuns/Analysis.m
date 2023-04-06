clc
clear
close all

addpath('../')

%%  Case 1

clear

load '2 Sources Crossing (different Amp, different range)';
Analyze

%% Case 2

clear

load '2 Sources Crossing (Same Amp, different Range).mat';
Analyze

%% Case 3 

clear

load '2 Sources Crossing (Different Amp, Same Range).mat';
Analyze

%% Case 4
%This one is interesting

clear

load '2 Sources Crossing (Same Amp, Same Range).mat';
Analyze

%% Case 5

clear

load '2 Sources Crossing (Same Amp, Slightly Different Range).mat';
Analyze

%% Case 6

clear

load '2 Sources Crossing (Slightly Different Amp, Slightly Different Range).mat';
Analyze

%% Case 7

clear

load '2 Sources (Crossing SNRs, Different Bearings, AntiSymmetric).mat';
Analyze

%% Case 8

clear

load '2 Sources (Crossing SNRs, Different Bearings, Asymmetric).mat';
Analyze

%% Case 9

clear

load '2 Sources (Crossing SNRs, Different Bearings, Asymmetric 2).mat';
Analyze

%% Case 10

clear

load '2 Sources (Crossing SNRs, Different Amp, Different Bearings, Asymmetric 2).mat';
Analyze

%% Case 11

clear

load '2 Sources (Crossing SNRs, Same Amp, Different Bearings, Asymmetric 3).mat';
Analyze
% Make SNR's very similar, bearing very different, don't cross paths

% Make SNR's cross 2 dB apart and cross. Don't cross bearing