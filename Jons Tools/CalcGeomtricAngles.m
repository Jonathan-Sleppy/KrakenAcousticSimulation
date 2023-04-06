function [bearing,elevation] = CalcGeomtricAngles(SourceLocation,ArrayDepth)
    x = SourceLocation(:,1);
    y = SourceLocation(:,2);
    z = SourceLocation(:,3) - ArrayDepth;
    bearing = rad2deg(atan2(y,x));
    elevation = rad2deg(atan2(z,sqrt(x.^2+y.^2)));
end