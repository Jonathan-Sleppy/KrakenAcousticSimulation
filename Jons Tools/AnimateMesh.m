function [] = AnimateMesh(BF_Out_All,azimuth,elevation,time)
    figure();
    for i = 1:length(time)
        subplot(2,1,1)
        cla
        mesh(azimuth,elevation,BF_Out_All(:,:,i,1))
        xlabel('Azimuth (degrees)')
        ylabel('Elevation (degrees)')
        zlabel('Beamformer Output')
        title('Eigenvector 1')
        subplot(2,1,2)
        cla
        mesh(azimuth,elevation,BF_Out_All(:,:,i,2))
        xlabel('Azimuth (degrees)')
        ylabel('Elevation (degrees)')
        zlabel('Beamformer Output')
        title('Eigenvector 2')
        pause(0.5)
    end
end