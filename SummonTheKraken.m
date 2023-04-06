function [pressure,range,depth] = SummonTheKraken(Range,SourceDepth,ReceiveDepths,Freq,M,N)

Basic_KParams
kparms.freqs = Freq;
kparms.rec_ranges = Range/1000;
kparms.src_depth = SourceDepth;
kparms.rec_depths = ReceiveDepths;
[kparms] = Mk_Kraken_Input(kparms);
!/Users/jonathansleppy/Desktop/NEAR/at/bin/kraken.exe KRAKEN_MAT
[ Modes ] = read_modes_bin( 'KRAKEN_MAT.mod', kparms.freqs );
phi = Modes.phi;
a = interp1(Modes.z,phi,kparms.src_depth);
b = interp1(Modes.z,phi,kparms.rec_depths);
k = Modes.k;
Pressure = zeros(M,N);

for i = 1:M
    for j = 1:N
        modesum = 0;
        for q = 1:length(k) %loop through modes
            R = norm([Range(i*j),SourceDepth-ReceiveDepths(i*j)]);
            modesum = modesum + (1/sqrt(k(q)))*a(q)*b(i*j,q)*exp(-1i*k(q)*R)/sqrt(2*pi*R);
        end
        Pressure(i,j) = 2*pi*modesum*exp(1i*3*pi/4);
    end
end
pressure = Pressure;
range = Range;
depth = ReceiveDepths;

end