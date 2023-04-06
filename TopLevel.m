clc
clear
close all

%% BellHop set-up options

BHparms.runtitle = 'BellHop: Test Case';            % (1) Run Title
BHparms.freqs = 200;                                % (2) Source Frequency
BHparms.TopOpt = 'CVWT';                             % (4) Top Options
BHparms.src_depth = 50;                             % (7) Source Depth
BHparms.ssp = [0 1500;100 1480;500 1480];                %water column sound speed depth vs speed
cref = 1500;
[vr,rhor,delta,Refl,theta] = SedType_to_ReflCoef('Very Fine Sand',cref);
BHparms.BotOpt = 'A*';                              %a* uses just properties in BHparms.bottom and * for bathymetry file
%BHparms.BotOpt = 'F*';                              %F* uses reflection loss table and * for bathymetry file
BHparms.bottom = [1 0 vr*cref 0 rhor delta*54.58 0];    % (6a) Bottom Halfspace Properties
%BHparms.bottom = [1 0 1500 0 1 0 0];    % (6a) Bottom Halfspace Properties
%BHparms.rec_depths = [0:1:600];                   % (7) Reciever Depth Grid
%BHparms.rec_ranges = [0.01:.02:10];                  % (7) Reciever Range Grid
BHparms.rec_depths = 0:1:600;                   % (7) Reciever Depth Grid
BHparms.rec_ranges = 0.01:.01:10;                  % (7) Reciever Range Grid
BHparms.nbeams = 2000;                              % (9) Number of Beams
BHparms.angles = [-89 89];                          % (9) Beam Angles
%Bathymetry file inputs
BHparms.BtmInterp = 'L';                           %bathymetry table range (km vs depth (m))
BHparms.bathy = [0.00	300;  
10.00	600];


% Run type enumerations: 'a' = arrivals, 'C' = coherent pressure, 'I' =
% incoherent pressure

BHparms.RunType = 'a';

%% Analysis Options

Options.plotray = 0;
Options.ConvolveArrivals = 1;
Options.PlotOutputs = 1;


%% Convolution Options

Options.SourceSignal.type = 'LFM';
Options.SourceSignal.WindowType = 'Hanning';
Options.SourceSignal.fc = 4000;
Options.SourceSignal.PulseLength = 3;
Options.SourceSignal.Bandwidth = 4000;
Options.SourceSignal.IRMeasurementLocation = [1e3,200];  %[Range,Depth]


%% Run the Sim

RunBellhopSim(BHparms);


%% Handle the outputs
 
[~,OutputData] = ParseOutputs(Options,BHparms);

%% Make Some Audio Players

audio = MakeAudioPlayers(OutputData);

%% Try Deconvolution Methods

DirectDeconvData = DirectDeconvolution(OutputData);
MatchedDeconvData = MatchedDeconvolution(OutputData);

[plots,DirectData,MatchedData] = CheckImpulseResponses(DirectDeconvData,MatchedDeconvData,OutputData);
matchedaudio = MakeAudioPlayers(MatchedData);
directaudio = MakeAudioPlayers(DirectData);

%% Function Dump

function [plothandle,DirectData,MatchedData] = CheckImpulseResponses(DirectDeconvData,MatchedDeconvData,OutputData)
    Source = OutputData.SourceSignal;
    MatchedIR = MatchedDeconvData.ImpulseResponse;
    DirectIR = DirectDeconvData.ImpulseResponse;
    MatchedOutput = ConvolveSynthesizedSignal(Source,MatchedIR);
    DirectOutput = ConvolveSynthesizedSignal(Source,DirectIR);
    MatchedData = OutputData;
    MatchedData.ConvolvedSignal = MatchedOutput;
    DirectData = OutputData;
    DirectData.ConvolvedSignal = DirectOutput;
    DirectData.ImpulseResponse = DirectIR;
    MatchedData.ImpulseResponse = MatchedIR;
    plothandle.direct = MakePlots(DirectData,'Direct Deconvolution');
    plothandle.matched = MakePlots(MatchedData,'Matched Filter Deconvolution');
end

function [] = RunBellhopSim(BHparms)

    %Set up inputs and preliminaries
    getpaths
    [~] = Write_BellhopEnvFile(BHparms);    

    %Actually run the executable
    !/Users/jonathansleppy/Desktop/NEAR/at/bin/bellhop.exe BHparms_MAT
    
end

function [plothandle,OutputData] = ParseOutputs(Options,BHparms)
    if Options.plotray == 1
        ray = readmatrix('MunkB_eigenray.ray','FileType','text');
        plothandle = plot2Deigenray(ray);
    end
    if Options.PlotOutputs == 1
        switch BHparms.RunType
            case 'CB'
                [ ~, ~, ~, ~, ~, Pos, pressure ] = read_shd_bin( 'BHparms_MAT.shd' );
                plothandle = PlotPressureField(Pos,pressure);
                OutputData.Pressure = pressure;
            case 'a'
                [ Arr, Pos ] = read_arrivals_bin( 'BHparms_MAT.arr' );
                if Options.ConvolveArrivals == 1
                    SourceSignal = SynthesizeSignal(Options.SourceSignal);
                    ImpulseResponse = ConstructIR(Arr,Pos,Options.SourceSignal.IRMeasurementLocation,SourceSignal.fs);
                    ConvolvedSignal = ConvolveSynthesizedSignal(SourceSignal,ImpulseResponse);
                    OutputData.ConvolvedSignal = ConvolvedSignal;
                    OutputData.SourceSignal = SourceSignal;
                    OutputData.ImpulseResponse = ImpulseResponse;
                    plothandle = MakePlots(OutputData,[]);
                end
        end
    end
end

function getpaths()
    addpath('Matlab/ReadWrite')
end

function [h] = plot2Deigenray(ray)
    x = ray(:,1);
    z = ray(:,2);
    h = figure();
    plot(x,z,'linewidth', 3)
    axis ij
    grid on
    xlabel('Range (m)')
    ylabel('Depth (m)')
    title('Eigenray')
end

function [q] = PlotPressureField(Pos,pressure)
    pr = squeeze(pressure);
    q = figure();
    pcolor(Pos.r.r./1000,Pos.r.z,-10*log10(abs(pr).^2));shading('flat')
    set(gca,'YDir','Reverse')
    colormap(flipud(jet))
    set(gca,'Fontsize',14)
    xlabel('Range (km)','Interpreter','latex')
    ylabel('Depth (m)','Interpreter','latex')
    title('Transmission Loss (dB)','Interpreter','latex')
    clim([20,80])
    colorbar
end

function [q] = PlotArrivals(Arr,Pos)
    rr = Pos.r.r;
    rd = Pos.r.z;
    Amp = Arr(1000,6).A;
    Del = real(Arr(1000,6).delay);
    q = figure();
    stem(Del,real(Amp))
    set(gca,'Fontsize',14)
    xlabel('Delay (s)','Interpreter','latex')
    ylabel('Arrival Amplitude ($\mu {\rm Pa}$)','Interpreter','latex')
    title('Ray Arrivals','Interpreter','latex')
end

function [ImpulseResponse] = ConstructIR(Arr,Pos,TestLocation,fs)
    r = Pos.r.r;
    z = Pos.r.z;
    delta = zeros(length(r),length(z));
    for i = 1:length(r)
        for j = 1:length(z)
            delta(i,j) = norm([r(i),z(j)]' - TestLocation');
        end
    end
    [I,J] = find(delta == min(min(delta))); 
    LocationArr = Arr(I,J);
    Delay = LocationArr.delay;
    Amplitude = LocationArr.A;
    TimeVector = 0:1/fs:max(real(Delay));
    impresp = zeros(size(TimeVector));
    for i = 1:length(Amplitude)
        [~,Q] = min(abs(TimeVector - Delay(i)));
        impresp(Q) = Amplitude(i);
    end
    ImpulseResponse = impresp;
end

function [SourceSignal] = SynthesizeSignal(Options)
    switch Options.type
        case 'LFM'
            SourceSignal = SynthesizeLFM(Options);
        case 'CW'
            SourceSignal = SynthesizeCW(Options);
    end
    switch Options.WindowType
        case 'Hanning'
            window = hanning(length(SourceSignal.TimeSeries));
        case 'Tukey'
            window = tukeywin(length(SourceSignal.TimeSeries));
        case 'None'
            window = ones(length(SourceSignal.TimeSeries),1);
    end
    SourceSignal.TimeSeries = window.*SourceSignal.TimeSeries';

end

function [SourceSignal] = SynthesizeLFM(Options)
    fc = Options.fc;
    Bandwidth = Options.Bandwidth;
    fs = 2*(fc + Bandwidth/2);
    deltat = 1/fs;
    T = Options.PulseLength;
    t = 0:deltat:T-deltat;
    SourceSignal.TimeSeries = exp(1i*(Bandwidth/(2*T)*t.^2 + (fc - Bandwidth/2)*t));
    SourceSignal.fs = fs;
    SourceSignal.time = t;
end

function [AudioPlayers] = MakeAudioPlayers(SignalData)
    fs = SignalData.SourceSignal.fs;
    ImpulseResponse = real(SignalData.ImpulseResponse);
    ImpulseResponse = ImpulseResponse/max(ImpulseResponse);
    SourceSignal = real(SignalData.SourceSignal.TimeSeries);
    SourceSignal = SourceSignal/max(SourceSignal);
    ConvolvedSignal = real(SignalData.ConvolvedSignal.TimeSeries);
    ConvolvedSignal = ConvolvedSignal/max(ConvolvedSignal);
    AudioPlayers.SourceSignal = audioplayer(SourceSignal,fs);
    AudioPlayers.ImpulseResponse = audioplayer(ImpulseResponse,fs);
    AudioPlayers.ConvolvedSignal = audioplayer(ConvolvedSignal,fs);
end

function [ConvolvedSignal] = ConvolveSynthesizedSignal(SourceSignal,ImpulseResponse)
    SourceTime = SourceSignal.time;
    fs = SourceSignal.fs;
    Source = SourceSignal.TimeSeries;
    ConvolvedSignal.Time = [SourceTime,max(SourceTime)+1/fs:1/fs:(length(SourceTime)+length(ImpulseResponse)-2)/fs];
    ConvolvedSignal.TimeSeries = conv(Source,ImpulseResponse);
end

function [h] = MakePlots(OutputData,plottitle)
    SourceData = OutputData.SourceSignal.TimeSeries;
    SourceTime = OutputData.SourceSignal.time;
    IRData = OutputData.ImpulseResponse;
    ConvolvedData = OutputData.ConvolvedSignal.TimeSeries;
    ConvolvedTime = OutputData.ConvolvedSignal.Time;
    h = figure();
    subplot(1,3,1)
    plot(SourceTime,real(SourceData))
    title('Source Signal')
    subplot(1,3,2)
    plot(1:length(IRData),real(IRData),'lineWidth',2)
    title('Impulse Response')
    subplot(1,3,3)
    plot(ConvolvedTime,real(ConvolvedData))
    title('Convolved Signal')
    if ~isempty(plottitle)
        sgtitle(plottitle)
    end
end

function [OutputData] = MatchedDeconvolution(OutputData)
    Received(:) = OutputData.ConvolvedSignal.TimeSeries;
    Transmitted(:) = OutputData.SourceSignal.TimeSeries;
    
    %Make a zero pad
    l = length(Received) - length(Transmitted);
    zeropad = zeros(1,l);
    Transmitted = [zeropad,Transmitted];
    
    %Do the matched filtering
    ReceivedF = fft(Received);
    TransmittedF = fft(Transmitted);
    MatchedF = conj(TransmittedF).*ReceivedF;
    OutputData.ImpulseResponse = ifft(MatchedF);
end

function [OutputData] = DirectDeconvolution(OutputData)
    Received = OutputData.ConvolvedSignal.TimeSeries;
    Transmitted = OutputData.SourceSignal.TimeSeries;
    
    %Make a zero pad
    l = length(Received);
    l = l - length(Transmitted);
    zeropad = zeros(l,1);
    Transmitted = [zeropad;Transmitted];
    
    %Do the direct deconvolution
    OutputData.ImpulseResponse = ifft(fft(Received)./fft(Transmitted));

end