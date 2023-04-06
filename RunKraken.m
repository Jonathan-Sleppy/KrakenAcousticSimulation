function [pr,rr,rd] = RunKraken(R,SD,ZD,FR)
%% Basic_Params
%Run pressure field 
Basic_KrakenParams
% [C,IA,IC] = unique(R);
% C = C + 1000*[0:1:length(C)-1]'.*1e-6;
% kparms.rec_ranges = C/1000;
kparms.freqs = FR;
kparms.rec_ranges = R/1000;
kparms.src_depth = SD;
kparms.rec_depths = ZD;    %receiver depths
[kparms] = Mk_Kraken_Input(kparms);
!/Users/jonathansleppy/Desktop/NEAR/at/bin/kraken.exe KRAKEN_MAT
[ Modes ] = read_modes_bin( 'KRAKEN_MAT.mod', kparms.freqs );%needed to edit this to remove persistance of fid as a variable
%[~,isd] = min(abs(kparms.src_depth - Modes.z));
phi = Modes.phi;
a = interp1(Modes.z,phi,kparms.src_depth);
%a = phi( isd, : );

b = interp1(Modes.z,phi,kparms.rec_depths);
%phi = phi * diag( a, 0 );	% scale modes by a

% ******************************************************
% form pressure field
% ******************************************************
ck = Modes.k;
tmp = zeros(size(b,1),Modes.M);
pr2 = zeros(size(b,1),length(R));
for m = 1:length(R)
    for n = 1:Modes.M
        tmp(:,n) = (1/sqrt(ck(n))).*a(n).*b(:,n).*exp(-1i*ck(n).*R(m))./sqrt(2*pi*R(m));
    end
    pr2(:,m) = (2*pi*sum(tmp,2).*exp(i*3*pi/4));
end
%Using Field program. Problematic because the precision gets lost when
%reading in the range values
% Mk_Field(kparms);
% !rm -f KRAKEN_MAT.shd
% eval(['!/Users/jonathansleppy/Desktop/NEAR/at/bin/field.exe KRAKEN_MAT'])
% [ Title, PlotType, freqVec, freq0, atten, Pos, pressure ] = read_shd_bin( 'KRAKEN_MAT.shd' );
% prf = squeeze(pressure);
% rr = Pos.r.r;
% rd = Pos.r.z;
pr = pr2;
rr = R;
rd = kparms.rec_depths;
% 
% figure(1)
% pcolor(Pos.r.r./1000,Pos.r.z,10*log10(abs(pr).^2));shading('flat')
% set(gca,'YDir','Reverse')
%%