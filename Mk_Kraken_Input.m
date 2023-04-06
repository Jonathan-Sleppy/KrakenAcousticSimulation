%PRECOMP_MODES_CM
%       Precomputes mode functions using KRAKEN. Writes a KRAKEN input
%       file calculates the mode functions and reads them into matlab, makes coupling matrices.
%
%USAGE
%     [kparms] = Mk_Kraken_Input(kparms);
%
%OUTPUTS
%      kparms            ==>  Data structure with modes and wavenumbers.
%                             Needed to compute pressure fields
%
%INPUTS
%      kparms.freqs        ==>  Frequency;
%      kparms.depth_grid   ==>  Grid to compute mode coupling matices and max depth
%      kparms.ssp.wat      ==>  Water sound speed profile
%      kparms.bottom.props ==>  Bottom properties 
%      kparms.bottom_type  ==>  Bottom type (Half-space either Acoustic
%                                                half-space or vacuum)
%      kparms.min_depth    ==>  Minimum bathymetry to precompute
%      kparms.max_depth    ==>  Maximum bathymetry to precompute
%      kparms.complex      ==>  Complex (yes or no) or real axis 
%      kparms.rec_depths   ==>  Receiver depths in mode functions to store
%                      (Specify any that might be used in field calculations)
%      kparms.src_depth    ==>  Source depths in mode functions to store
%      (Specify any that might be used in field calculations)
%      kparms.KrOpt        ==>  Kraken options (for top BC and attenuation)
%      kparms.surface_rough==>  RMS surface roughness (m)
%EXAMPLE
%  kparms.freqs = 200;
%  kparms.depth_grid = [250 1];         %zmax dz
%                                     %water sound speed profile
%  kparms.ssps.wat =[
%             0   1500;               %depth ss
%             70  1490;               %depth ss
%             75  1490;
%             150 1490];
%  kparms.bottom.props = [1 0  1550 0 1.5 0.1 0;     %layer number, depth comp-speed,shear-speed, density, comp-attn ,shear-attn
%               1 20 1750 0 1.5 0.1 0;
%               1 30 1800 0 1.5 0.1 0;
%               2 0  1850 0 1.5 0.1 0;
%               2 50 1850 0 1.5 0.1 0;
%               2 60 1860 0 1.5 0.1 0;
%               3 0  1950 0 1.5 0.1 0;
%               3 20 1975 0 1.5 1.0 0
%               4 0  1975 0 1.5 1 0;
%               4 50 1975 0 1.5 10 0];
%  kparms.min_depth = 20;
%  kparms.max_depth = 200;
%%Bottom type below sediment layers
%  kparms.bottom_type = 'A';
%%complex or real modes (only real is working now)
%  kparms.complex = 'n';
%%These are the depths to be used both for array depth locations and target depths
%  kparms.rec_depths = [0:1:250];    %receiver depths
%%This doesn't really matter if just calculating modes, but if you specify
%%something it will make sure mode functions are stored at that depth
%  kparms.src_depth = [30]; %source depth
%%kraken options
%  kparms.KrOpt = 'CVW';
%  [kparms] = precomp_modes_cm(kparms);
%
%%This can be followed by the following to get the field
%
%%can change any of these w/o re-running precomp_modes_cm (but don't go out
%%of bathymetry min max range
%   kparms.src_depth = [20]; %source depth
%   kparms.rec_ranges = [0.025:.025:10];  %rmin delta r rmax
%   kparms.bathy = [0 40; 4 40; 10 150];
%   kparms.surface_rough = 0;
%   [Prk,ranges,depths,kparms] = juggler_cm(kparms);


%
% M. Siderius


function [kparms] = precomp_modes(kparms);

if(ispc == 1)
    !del KRAKEN_MAT.env KRAKEN_MAT.prt
    % !del MODFIL*
else
    !rm -f KRAKEN_MAT.env KRAKEN_MAT.prt
    % !rm -f MODFIL*
end;
%option to store modes as int16's instead of double to save memorty (or no to not do)
if(isfield(kparms,'store_ints') == 0)
    %kparms.store_ints = 'd';
    %kparms.store_ints = 'y';
    kparms.store_ints = 'n';
end
%subsample from depth grid every depth_samp point 
depth_samp = 1;
%delta z for bathymetry step size (m)
if(isfield(kparms,'pts_per_wlngth') == 0)
    dz = (1500/kparms.freqs)./10;
else 
    dz = (1500/kparms.freqs)./kparms.pts_per_wlngth;
end;
if(isfield(kparms,'no_diags') == 0)
    no_diags = 1e99;
else
    no_diags = kparms.no_diags;
end;
if(isfield(kparms,'bottom_type') == 0)
    kparms.bottom_type = 'A';
end;
if(isfield(kparms,'env') == 0)
    kparms.env = [0 1 1;1e9 1 1];
    kparms.env_step = kparms.env;
end;
if(isfield(kparms,'attn_wat') == 0)
    attn_wat = 0;
else
    attn_wat = kparms.attn_wat;
end;

for nssp = 1:length(kparms.ssps)
    for nbott = 1:length(kparms.bottom)
        kparms.ssp = kparms.ssps(nssp).wat;
        kparms.bott = kparms.bottom(nbott).props;
        if(isfield(kparms,'min_depth')  == 0)
            if(isfield(kparms,'bathy') == 0)
                disp('Define min, max depths or bathymetry')
                return
            else
                kparms.min_depth = min(kparms.bathy(:,2));
            end;
        end;
        if(isfield(kparms,'max_depth')  == 0)
            if(isfield(kparms,'bathy') == 0)
                disp('Define min, max depths or bathymetry')
                return
            else
                kparms.max_depth = max(kparms.bathy(:,2));
            end;
        end;
        kparms.modes_into_mem = 'n';
        kparms.options = 'RC';
        kparms.dispbot = 0;
        bottype = kparms.bottom_type;
        NMEDIA = max(kparms.bott(:,1));
        min_dp = kparms.min_depth;
        max_dp = kparms.max_depth;
        lnth = length(min_dp:dz:max_dp);
        %if(isfield(kparms,'bathy') == 0)
        %find min max depths and decide on number of sectors
        %guaranteed to give a point at max depth and min_depth
        bathy = linspace(min_dp,max_dp,lnth);
        %        bathy = logspace(log10(round(1500/kparms.freqs)),log10(max_dp),100);
        %        bathy = logspace(log10(min_dp),log10(max_dp),50);
        %        bathy = linspace(min_dp,max_dp,300);
        kparms.bathy_min_max = bathy;
        %not guaranteed to give a point at max depth
        %   bathy = min_dp:dz:max_dp;
        %else
        % bathy = interp1(1:length(kparms.bathy(:,2)),kparms.bathy(:,2),linspace(1,length(kparms.bathy(:,2)),lnth));
        % kparms.bathy_min_max = bathy;
        % kparms.RPROF = interp1(1:length(kparms.bathy(:,1)),kparms.bathy(:,1),linspace(1,length(kparms.bathy(:,1)),lnth));
        % kparms.MPROF = 1:length(kparms.RPROF);
        %end
        
        
        %edit sound speed so there is a good point at the top and bottom
        SSP = kparms.ssp;
        if(SSP(1,1) ~= 0)
            SSP(1,1) = 0.0;
            SSP(1,2) = SSP(2,2);
        end;
        if(SSP(length(SSP)) < bathy(1))
            %ssps_tmp = interp1(SSP(end-1:end,1),SSP(end-1:end,2),[bathy(1)],'linear','extrap');
            ssps_tmp = SSP(length(SSP),2);
            SSP(length(SSP)+1,1) = bathy(1);
            SSP(length(SSP),2) = ssps_tmp;
        end;
        ssps = SSP;
        
       % disp('  Examine KRAKEN_MAT.ENV to verify model parameters.')
        %get the basic parameters for Kraken input file
        if(isempty(bottype) == 1)
            bottype = 'V';
        end;
        NRDM = kparms.depth_grid;
        depth = [0:NRDM(2):NRDM(1)];
        FREQ = kparms.freqs;
        SRC_DEPTH = depth(2); %doesn't matter what this is set to
        dispbot = kparms.dispbot;
        NPROF = length(bathy);
        kparms.NPROF = NPROF;
        %open .env file and write the input
        fid1 = -1;
        while(fid1 == -1);
            fid1 = fopen('KRAKEN_MAT.env','w');
        end;
        %loop on #input sectors
        
        for iprof = 1:NPROF
            txt = [''' KRAKEN_MAT ''','  ', int2str(iprof)];
            fprintf(fid1,'%c', txt);
            fprintf(fid1,'\n');
            %first all the header information
            fprintf(fid1,'%4.4f\n', FREQ);
            if(size(kparms.bott,1) == 1)
                fprintf(fid1,'%i\n', NMEDIA);
            else
                fprintf(fid1,'%i\n', NMEDIA+1);
            end
            if(isfield(kparms,'KrOpt') == 1)
                options = ['''',kparms.KrOpt,''''];
            else
                options = '''CVW''';
            end
            fprintf(fid1,'%c', options);
            fprintf(fid1,'\n');
            
            %%%%%%%%%%
            %write out the sound speed in water block
            fprintf(fid1,'%i %4.4f %4.4f\n', [0 0 bathy(iprof)]);
            %               fprintf(fid1,'%i %4.4f %4.4f\n', [round(bathy(iprof)./(1500/kparms.freqs/64)) 0 bathy(iprof)]);
            %make sure the sound speed is extended to the top
            bott = kparms.bott(:,1:6,1);
            ssps_tmp = 0;
            if(ssps(1,1) > 0)
                fprintf(fid1,'%4.4f %4.4f %4.4f %4.4f %.4e', [0.0 ssps(1,2) 0.0 1.0 attn_wat]);
                fprintf(fid1,'%c', ' /');
                fprintf(fid1,'\n');
            end
            for n = 1:size(ssps,1);
                if(ssps(n,1) <= bathy(iprof))
                    fprintf(fid1,'%4.4f %4.4f %4.4f %4.4f %.4e', [ssps(n,1) ssps(n,2) 0.0 1.0 attn_wat]);
                    fprintf(fid1,'%c', ' /');
                    fprintf(fid1,'\n');
                    ssps_tmp = ssps(n,2);
                    dpt_tmp = ssps(n,1);
                    last = n;
                end
            end;
            %make sure the sound speed is extended to the bottom
            if((dpt_tmp - bathy(iprof)) < -.5e-4)
                if(length(ssps) > last)
                    W = (ssps(last+1,1) - bathy(iprof))./(ssps(last+1,1)-ssps(last,1));
                    ssps_tmp = W.*ssps(last,2) + (1-W).*ssps(last+1,2);
                    fprintf(fid1,'%4.4f %4.4f %4.4f %4.4f %.4e', [bathy(iprof) ssps_tmp 0.0 1.0 attn_wat]);
                else
                    ssps_tmp = interp1(ssps(end-1:end,1),ssps(end-1:end,2),[bathy(iprof)],'linear','extrap');
                    fprintf(fid1,'%4.4f %4.4f %4.4f %4.4f %.4e', [bathy(iprof) ssps_tmp 0.0 1.0 attn_wat]);
                end;
                fprintf(fid1,'%c', ' /');
                fprintf(fid1,'\n');
            end;
            
            %%%%%%%%%%
            %write out the sound speed in sediment block
            if(size(kparms.bott,1) > 1 )
                thkns(1) = 0;
                %keyboard
                for nm = 1:NMEDIA
                    layindxmax =  max(find(kparms.bott(:,1) == nm));
                    layindxmin =  min(find(kparms.bott(:,1) == nm));
                    thkns(nm+1) = kparms.bott(layindxmax,2);
                    if(nm < NMEDIA)
                        fprintf(fid1,'%i %4.4f %4.4f\n', [0 0.0 round((sum(thkns(1:nm+1))+bathy(iprof))*10e6)./10e6]);
                    else
                        fprintf(fid1,'%i %4.4f %4.4f\n', [0 0.0 max(kparms.depth_grid(1),round((sum(thkns(1:nm+1))+bathy(iprof))*10e6)/10e6)]);
                    end;
                    %check for a point at the top of the sediment layer
                    for n=layindxmin:layindxmax
                        fprintf(fid1,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f', [round((bathy(iprof)+kparms.bott(n,2)+sum(thkns(1:nm)))*10e6)/10e6 kparms.bott(n,3:7)] );
                        fprintf(fid1,'%c', ' /');
                        fprintf(fid1,'\n');
                    end;
                    if((nm == NMEDIA) & (bathy(iprof)+kparms.bott(n,2)+sum(thkns(1:nm))) < kparms.depth_grid(1))
                        fprintf(fid1,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f', [kparms.depth_grid(1) kparms.bott(n,3:7)] );
                        fprintf(fid1,'%c', ' /');
                        fprintf(fid1,'\n');
                    end;
                end;
            end;
            if(bottype == 'V')
                options = '''V'' 0.0';
                fprintf(fid1,'%c', options);
                fprintf(fid1,'%c\n', ' ');
            elseif(bottype == 'R')
                options = '''R'' 0.0';
                fprintf(fid1,'%c', options);
                fprintf(fid1,'%c\n', ' ');
            else
                options = '''A'' 0.0';
                fprintf(fid1,'%c', options);
                fprintf(fid1,'%c\n', ' ');
                if(size(kparms.bott,1) > 1)
                    fprintf(fid1,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f', [max(kparms.depth_grid(1),(bathy(iprof)+sum(thkns(1:nm+1)))) kparms.bott(n,3:7)] );
                else
                    fprintf(fid1,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f', [bathy(iprof) kparms.bott(1,3:7)] );              
                end
                fprintf(fid1,'%c', ' /');
                fprintf(fid1,'\n');
            end;
            if(isfield(kparms,'phsspeedlims') == 0)
                if(bathy(iprof) < 2000)
                    phspeed = [ 1400 (max(kparms.bott(:,3))+10000)];
                else
                    phspeed = [1400 max([ssps(:,2); ssps_tmp ;1550])];
                end
            else
                phspeed = kparms.phsspeedlims;
            end
            fprintf(fid1,'%4.1f %4.1f\n', phspeed);
            fprintf(fid1,'%i\n', 1000);
            fprintf(fid1,'%i\n', [1]);
            fprintf(fid1,'%4.2f', [SRC_DEPTH]);
            fprintf(fid1,'%c', ' /');
            fprintf(fid1,'%c\n', ' ');
            fprintf(fid1,'%i\n', [length(depth)]);
            fprintf(fid1,'%4.2f %4.2f', [0 NRDM(1)]);
            fprintf(fid1,'%c', ' /');
            fprintf(fid1,'%c\n', ' ');
        end;
        fclose(fid1);
%        disp('  Input file construction complete.')
    end
end

