function [Cylindrical] = ConvertToCylindricalCoords(Cartesian)
    Cylin = zeros(size(Cartesian));
    x = Cartesian(:,1);
    y = Cartesian(:,2);
    z = Cartesian(:,3);
    for i = 1:length(x)
        Cylin(i,1) = norm([x(i),y(i)]);
        Cylin(i,2) = atan2(y(i),x(i));
        Cylin(i,3) = z(i);
    end
    Cylindrical = Cylin;
end