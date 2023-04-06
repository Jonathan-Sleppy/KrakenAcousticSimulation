function [] = PlotSNR(SNR,time,plottitle)
    hold on
    for i = 1:length(SNR(:,1))
        color = enumeratecolors(i);
        plot(time(1:length(SNR(1,:))),real(SNR(i,:)),'LineWidth',3,'Color',color)
    end
    title([plottitle,' SNR'])
    ylim([-10,10])
    xlim([min(time),max(time)])
    xlabel('Relative Time')
    ylabel('Signal-to-Noise Ratio (dB)')
    hold off
    grid on
end