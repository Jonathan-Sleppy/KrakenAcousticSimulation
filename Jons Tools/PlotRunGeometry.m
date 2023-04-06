function [plothandle] = PlotRunGeometry(Receiver,SourceDepth)
    plothandle = figure();
    x = Receiver(:,1);
    y = Receiver(:,2);
    z = Receiver(:,3);
    hold on
    plot3(x,y,z,'o')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    plot3(0,mean(y) - 0.01*mean(y),SourceDepth,'*')

end