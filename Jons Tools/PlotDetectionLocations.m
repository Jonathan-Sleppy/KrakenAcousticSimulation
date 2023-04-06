function [plothandle] = PlotDetectionLocations(DetectionLocations,TrueBearing,TrueElevation,time,plottitle)
    plothandle = figure();
    for i = 1:length(DetectionLocations(:,1,1))

        color = enumeratecolors(i);
        Bearing = squeeze(DetectionLocations(i,1,:));
        Elevation = squeeze(DetectionLocations(i,2,:));
        subplot(2,2,1)
        hold on
        plot(time,Bearing,'LineWidth',3,'Color',color)
        ylim([75,105])
        xlim([min(time),max(time)])
        hold off
        title([plottitle,' Calculated Bearing'])
        xlabel('Relative Time')
        ylabel('Bearing (deg)')
        subplot(2,2,2)
        hold on
        plot(time,TrueBearing(i,:),'LineWidth',3,'Color',color)
        ylim([75,105])
        xlim([min(time),max(time)])
        title([plottitle,' True Bearing'])
        xlabel('Relative Time')
        ylabel('Bearing (deg)')
        subplot(2,2,3)
        hold on
        plot(time,Elevation,'LineWidth',3,'Color',color)
        hold off
        title([plottitle,' Calculated Elevation'])
        xlabel('Relative Time')
        ylabel('Elevation (deg)')
        subplot(2,2,4)
        hold on
        plot(time,TrueElevation(i,:),'LineWidth',3,'Color',color)
        hold off
        title([plottitle,' True Elevation'])
        xlabel('Relative Time')
        ylabel('Elevation (deg)')
        
    end
end