function [figurehandle] = PlotTransmissionLoss(pressure,range,depth)
    figurehandle = figure();
    pcolor(range,depth,-10*log10(abs(pressure).^2));shading('flat')
    set(gca,'YDir','Reverse')
    colormap(flipud(jet))
    set(gca,'Fontsize',14)
    xlabel('Range (km)','Interpreter','latex')
    ylabel('Depth (m)','Interpreter','latex')
    title('Transmission Loss (dB)','Interpreter','latex')
    caxis([20,80])
    colorbar
end