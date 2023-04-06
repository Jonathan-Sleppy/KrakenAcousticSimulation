function [] = PlotBearing(DetectionLocations,TrueBearing,time,plottitle)
    hold on
    for i = 1:length(DetectionLocations(:,1,1))
        color = enumeratecolors(i);
        Bearing = squeeze(DetectionLocations(i,1,:));        
        plot(time(1:length(TrueBearing(1,:))),Bearing,'LineWidth',3,'Color',color)
        ylim([75,105])
        xlim([min(time),max(time)])
        plot(time(1:length(TrueBearing(1,:))),TrueBearing(i,:),'LineWidth',3,'Color',color,'LineStyle',':')
        ylim([75,105])
        xlim([min(time),max(time)])
        title([plottitle,' True and Estimated Bearings'])
        xlabel('Relative Time')
        ylabel('Bearing (deg)')
    end
    hold off
end
