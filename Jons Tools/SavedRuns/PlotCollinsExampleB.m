clc
clear
close all

addpath('../')

load 'Collins Example B w corrmaxalg 0.99.mat';
%load 'Collins Example B w corrmaxalg 0.95.mat';

DetectionLocations = squeeze(DetectionLocationsFull(:,1,:));
VertMeanCorrelation = squeeze(mean(VertPairCorrelation,1));
HorizMeanCorrelation = squeeze(mean(HorizPairCorrelation,1));
GeometricBearing = ReflectAnglesOnToRange(GeometricBearing);
VertMeanCorrelationTransformed = squeeze(mean(VertPairCorrelationTransformed,1));
CorrelationChanges = VertMeanCorrelationTransformed(1:4,:) - VertMeanCorrelation;
[ambiguity,timestamps] = ParseAmbiguity(ambiguity,time);



for i = 1:length(DetectionLocations(:,1))
    figure(1)
    color = enumeratecolors(i);
    hold on
    plot(time,SNRfull(i,:),'LineWidth',3,'Color',color)
    ylim([-10,10])
    title('SNR')
    ylabel('dB')
    xlabel('time')
    hold off
    grid on
    figure(2)
    hold on
    plot(time,DetectionLocations(i,:),'o','LineWidth',3,'Color',color)
    plot(time,GeometricBearing(i,:),'LineStyle','--','LineWidth',3,'Color',color)
    hold off
    grid on
    ylabel('degrees')
    xlabel('time')
    title('True and Estimated Bearing')
    figure(3)
    hold on
    plot(time,VertMeanCorrelation(i,:),'LineWidth',3,'Color',color)
    hold off
    grid on
    title('Vertical Subarray Correlation')
    ylim([0.5,1])
    figure(4)
    hold on
    plot(time,VertMeanCorrelationTransformed(i,:),'LineWidth',3,'Color',color)
    hold off
    grid on
    title('Transformed Vertical Subarray Correlation')
    ylim([0.5,1])
    figure(5)
    hold on
    plot(time,CorrelationChanges(i,:),'LineWidth',3,'Color',color)
    hold off
    grid on
    title('Change in Correlation from Transformation')
end

%% Correlation Coefficient


mesh(rad2deg(AngleRange),AlphaRange,ambiguity(:,:,3,3))
ylabel('Alpha')
xlabel('Phi (deg)')
zlabel('Correlation')

function [ambiguityout,timestamp] = ParseAmbiguity(ambiguity,time)
    [a,b,~,~] = size(ambiguity);
    x = 0;
        for i = 1:length(time)
            for q = 1:length(ambiguity(1,1,:,1))
                if ambiguity(:,:,q,i) ~= zeros(a,b)
                    x = x+1;
                    ambiguityout(:,:,:,x) = ambiguity(:,:,:,i);
                    timestamp(x) = time(i);
                end
            end
        end
end

function [reflectedangles] = ReflectAnglesOnToRange(Angles)
    greaterthaninds = find(Angles > 180);
    lessthaninds = find(Angles < 0);
    Angles(greaterthaninds) = 180 + Angles(greaterthaninds);
    Angles(lessthaninds) = -Angles(lessthaninds);
    reflectedangles = Angles;
end

%Plot 0-180 with and without correction