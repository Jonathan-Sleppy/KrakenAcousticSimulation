function [plothandle] = PlotArrayOrientation(R,plottitle,Translation)
    if ~isempty(Translation)
        v = ones(length(R(:,1)),1);
        R = (R' - Translation'*v')';
    end
    x = R(:,1);
    y = R(:,2);
    z = R(:,3);
    plothandle = figure();
    hold on
    quiver3([0,0,0],[0,0,0],[0,0,0],[1,0,0],[0,1,0],[0,0,1])
    plot3(x,y,z,'o','linewidth',3)
    grid on
    %axis equal
    title(plottitle)
end