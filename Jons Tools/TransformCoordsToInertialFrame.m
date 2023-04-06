function [NewCoords] = TransformCoordsToInertialFrame(Translation,Rotation,BodyFrameCoords)

    R = BodyFrameCoords';
    rc = Translation';
    v = ones(max(size(R)),1);
    Rprime = rc*v' + Rotation*R;
    NewCoords = Rprime';
end