function [plothandle] = plotSubarrayCorrelation(subarraycorrelation,time)
    plothandle = figure();
    hold on
    [a,b,c] = size(subarraycorrelation);
    eig1 = squeeze(subarraycorrelation(:,1,:));
    eig2 = squeeze(subarraycorrelation(:,2,:));
    
end