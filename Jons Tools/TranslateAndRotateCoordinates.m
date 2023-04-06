function [InertialCoordinates,DCM] = TranslateAndRotateCoordinates(yaw,pitch,roll,BodyFrameCoords,Translation)
%This function takes an array of 'BodyFrameCoords', Rotates to them to a
%desired orientation, given by yaw, pitch, and roll, and then translates
%them to an arbitrary position in space given by 'Translation'.
%BodyFrameCoords is an n by 3 array where the rows correspond to each point
%of interest on the body, and the three columns have the x-y-z body frame coordinates respictively

    DCM = CreateCompositeRotation(yaw,pitch,roll);
    InertialCoordinates = TransformCoordsToInertialFrame(Translation,DCM,BodyFrameCoords);
end