close all
if length(theta) > 1
    AnimateMesh(BF_Out_All,phi,theta,time);
else
    Animate2Dplot(BF_Out_All,PairCorrelation,phi,time,SourceInitial,SourceFinal,DetectionLocationsFull,GeometricBearing,SNRfull,[],Eigenvalues);
end
