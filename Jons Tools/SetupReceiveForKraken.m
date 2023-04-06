function [ReceiverRange,ReceiverDepth,TransformedR,DCM] = SetupReceiveForKraken(yaw,pitch,roll,ArrayRange,ArrayDepth,ArrayCoords,plotoption)
    Translation = [0,ArrayRange,ArrayDepth];
    [TransformedR,DCM] = TranslateAndRotateCoordinates(yaw,pitch,roll,ArrayCoords,Translation);
    TransformedRCylindrical = ConvertToCylindricalCoords(TransformedR);
    ReceiverRange = TransformedRCylindrical(:,1);
    ReceiverDepth = TransformedRCylindrical(:,3);
    if plotoption == 1
        [~] = PlotArrayOrientation(TransformedR,'Receiver Orientation',Translation);
    end
end