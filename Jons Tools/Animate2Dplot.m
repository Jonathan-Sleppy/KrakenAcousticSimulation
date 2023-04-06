function [] = Animate2Dplot(BF_Out_All,PairCorrelation,azimuth,time,SourceInitial,SourceFinal,DetectionLocations,TrueBearing,SNRfull,EigSNR,Eigenvalues)
    h = figure();
    h.Units = 'normalized';
    h.Position = [0.268023255813953,0.094444444444444,0.557267441860465,0.824305555555556];
    dt = time(2) - time(1);
    Velocity = SourceFinal - SourceInitial;
    SourceLocation = zeros([size(SourceInitial),length(time)+1]);
    SourceLocation(:,:,1) = SourceInitial;
    MeanCorrelation = squeeze(mean(PairCorrelation,1));
    
    EstBearing = squeeze(DetectionLocations(:,1,:));
    BearingError = abs(TrueBearing - EstBearing);

    for i = 1:length(time)
        subplot(2,5,1)
        cla
        plot(azimuth,BF_Out_All(1,:,i,1),'LineWidth',3)
        xlabel('Azimuth (deg)')
        ylabel('Beamformer Output')
        title('Eigenvector 1')
        grid on
        ylim([0,2*10^-3])
        subplot(2,5,2)
        cla
        plot(azimuth,BF_Out_All(1,:,i,2),'LineWidth',3)
        xlabel('Azimuth (deg)')
        ylabel('Beamformer Output')
        title('Eigenvector 2')
        grid on
        ylim([0,2*10^-3])
        subplot(2,5,6)
        cla
        bar(PairCorrelation(:,1,i))
        title('Eigenvector 1 Subarray Correlation')
        ylim([0,1])
        grid on
        subplot(2,5,7)
        cla
        bar(PairCorrelation(:,2,i))
        title('Eigenvector 2 Subarray Correlation')
        grid on

        subplot(2,5,10)
        cla
        hold on
        for q = 1:length(MeanCorrelation(:,1))
            color = enumeratecolors(q);
            plot(time(1:i),MeanCorrelation(q,1:i),'LineWidth',3,'Color',color)
        end
        xlim([min(time),max(time)])
        ylim([min(min(MeanCorrelation)),max(max(MeanCorrelation))])
        hold off
        title('Average Subarray Correlation')
        grid on
        
        subplot(2,5,8)
        bar(Eigenvalues(:,i))
        title('Eigenvalue Magnitude')
        xlabel('Eigenvalue')
        ylabel('Magnitude')
        ylim([0,max(max(Eigenvalues))])
        grid on
        
        % Plot SNR's
        subplot(2,5,9)
        cla
        hold on
        %PlotSNR(SNRfull(:,1:i),time,'Eigenvalues')
        PlotSNR(SNRfull(:,1:i),time,'Full Array')
        hold off

        % Plot true and Estimated Bearings
        subplot(2,5,4)
        cla
        PlotBearing(DetectionLocations(:,:,1:i),TrueBearing(:,1:i),time,'Full Array')
        grid on

        % Plot Bearing Error
        subplot(2,5,5)
        cla
        hold on
        for q = 1:length(BearingError(:,1))
            color = enumeratecolors(q);
            plot(time(1:i),BearingError(q,1:i),'LineWidth',3,'Color',color)
        end
        hold off
        grid on
        ylim([0,max(max(BearingError))])
        xlim([min(time),max(time)])
        xlabel('Relative Time')
        ylabel('Error (degrees)')
        title('Bearing Estimation Error')

        % Plot Problem Geometry
        subplot(2,5,3)
        cla
        hold on
        for q = 1:length(SourceLocation(:,1,1))
            color = enumeratecolors(q);
            x = squeeze(squeeze(SourceLocation(q,1,1:i)));
            y = squeeze(squeeze(SourceLocation(q,2,1:i)));
            z = squeeze(squeeze(SourceLocation(q,3,1:i)));
            plot3(x,y,z,'LineWidth',3,'Color',color)
        end
        plot3(0,0,z,'*','LineWidth',3,'Color','k')
        hold off
        grid on
        title('Kinematics')
        xlabel('x (km)')
        ylabel('y (km)')
        zlabel('z (m)')
        SourceLocation(:,:,i+1) = SourceLocation(:,:,i) + dt*Velocity;
        pause(0.01)
    end
end