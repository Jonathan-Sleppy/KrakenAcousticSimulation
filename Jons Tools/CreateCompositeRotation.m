function [Rotation] = CreateCompositeRotation(yaw,pitch,roll)


    LambdaRoll = [ 1,     0 ,       0;
                   0, cos(roll), -sin(roll);
                   0, sin(roll), cos(roll)];

    LambdaPitch = [cos(pitch), 0, sin(pitch);
                       0,      1,    0;
                  -sin(pitch), 0, cos(pitch)];

    LambdaYaw = [cos(yaw), -sin(yaw), 0;
                 sin(yaw),  cos(yaw), 0;
                    0,         0,     1];

    Rotation = LambdaYaw*LambdaPitch*LambdaRoll;
end