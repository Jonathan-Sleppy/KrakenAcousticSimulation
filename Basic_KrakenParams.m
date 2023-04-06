kparms.depth_grid = [800 1];         %zmax dz
kparms.ssps.wat =[0 1500;800 1500];
kparms.bottom.props = [1 0 1700 0 1.5 0.5 0;
     1 50  1700 0 1.5 0.5 0;
     1 100 1700 0 1.5 0.5 0;
     1 250 1700 0 1.5 0.5 0;
     1 300 1700 0 1.5 0.5 0;
     1 400 1700 0 1.5 20 0];      %layer number, depth comp-speed,shear-speed, density, comp-attn ,shear-attn  
kparms.min_depth = 800;
kparms.max_depth = 800;
kparms.bottom_type = 'V';
kparms.complex = 'y';
%kparms.rec_depths = [50];    %receiver depths
kparms.KrOpt = 'CVWT';
kparms.phsspeedlims = [1000 1e9];

