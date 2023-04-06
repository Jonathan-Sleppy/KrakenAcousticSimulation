clc
clear
close all

DummyR = [0 1 1;
          0 -1 1;
          0 1 -1;
          0 -1 -1];

Translation = [1,0,0];

yaw = deg2rad(45);
pitch = deg2rad(15);
roll = deg2rad(10);

h = PlotArrayOrientation(DummyR,'Unmoved');

TransformedR = TranslateAndRotateCoordinates(yaw,pitch,roll,DummyR,Translation);
q = PlotArrayOrientation(TransformedR,'Transformed');











