



DetectionLocationsFull = zeros(numdetects,2,length(time));
DetectionLocationsVert = zeros(numdetects,2,length(time));
DetectionLocationsHoriz = zeros(numdetects,2,length(time));
SNRfull = zeros(length(SourceLocation(:,1)),length(time));
SNReig = zeros(length(SourceLocation(:,1)),length(time));
GeometricBearing = zeros(length(SourceLocation(:,1)),length(time));
GeometricElevation = zeros(length(SourceLocation(:,1)),length(time));
HorizPairCorrelation = zeros(m-1,length(SourceLocation(:,1)),length(time));
VertPairCorrelation = zeros(n-1,length(SourceLocation(:,1)),length(time));
VertPairCorrelationTransformed = zeros(n-1,length(SourceLocation(:,1)),length(time));
BF_Out_All = zeros(length(theta),length(phi),length(time),numdetects);
Eigenvalues = zeros(numdetects+extraeigs,length(time));
ambiguity = zeros(length(AlphaRange),length(AngleRange),length(numdetects),length(time));

for q = 0:length(time) - 1 % Loop over time
    
    %Integrate Position
    SourceLocation = SourceLocation + (q>0)*dt*Velocity;

    % Reset K for new time step
    K = sigma*eye(m*n) + zeros(n*m,n*m);
    
    % Do a bunch of coordinate transformations to properly feed the Kraken
    SourceCylin = ConvertToCylindricalCoords(SourceLocation);
    SourceRange = SourceCylin(:,1);
    SourceDepth = SourceCylin(:,3);
    yaw = (pi/2*ones(size(SourceDepth)) - SourceCylin(:,2));
    [GeometricBearing(:,q+1),GeometricElevation(:,q+1)] = CalcGeomtricAngles(SourceLocation,ArrayDepth);

    
    for i = 1:length(SourceRange)   % Loop over sources
        [ReceiveRange,ReceiveDepth,TransformedArray,DCM] = SetupReceiveForKraken(yaw(i),pitch,roll,SourceRange(i),ArrayDepth,Arraymat,0);
        if plotoption == 1
            PlotRunGeometry(TransformedArray,SourceDepth(i));
        end
        %[pressure,rr,rd] = SummonTheKraken(ReceiveRange,SourceDepth(i),ReceiveDepth,SourceFreq,m,n);
        [pressure,rr,rd] = RunKraken(ReceiveRange,SourceDepth(i),ReceiveDepth,SourceFreq);
        
        pressure = SourceAmplitude(i)*diag(pressure);
        
        PressureCov = pressure*pressure';

        K = K + PressureCov;
        SNRfull(i,q+1) = 10*log10(trace(PressureCov + sigma*eye(m*n))/(m*n*sigma));
    end
    
    %% Multi-Valued Bartlett Processing
    
    %Eigen decomposition

    [Vfull,Dfull] = eigs(K,numdetects + extraeigs);

    for i = 1:numdetects

        %Store the eigenvalues
        Eigenvalues(:,q+1) = diag(Dfull);
        
        if computesubarraycorrelation == 1
            for x = 1:m-1
                [~,~,HorizCorInd,AdjHorizCorInd] = GetSubarrayIndex(m,n,x);
                v = Vfull(HorizCorInd,i);
                v = v/norm(v);
                vnext = Vfull(AdjHorizCorInd,i);
                vnext = vnext/norm(vnext);
                HorizPairCorrelation(x,i,q+1) = abs(v'*vnext);
            end
            [VertPairCorrelation(:,i,q+1),MeanVertPairCorrelation] = CalcMeanVertPairCorrelation(Vfull(:,i),m,n);
            VertPairCorrelationTransformed(:,i,q+1) = VertPairCorrelation(:,i,q+1);
        end
    end
    if ApplyTransformation == 1
        for i = 1:numdetects
            [~,MeanVertPairCorrelation] = CalcMeanVertPairCorrelation(Vfull(:,i),m,n);
            if MeanVertPairCorrelation < CorrelationThreshold
                [Vfull(:,i),Vfull(:,i+1),ambiguity(:,:,i,q+1)] = MaximizeCorrelation(Vfull(:,i:i+1),AlphaRange,AngleRange,m,n);
                [VertPairCorrelationTransformed(:,i,q+1),~] = CalcMeanVertPairCorrelation(Vfull(:,i),m,n);
                [VertPairCorrelationTransformed(:,i+1,q+1),~] = CalcMeanVertPairCorrelation(Vfull(:,i+1),m,n);
            end
        end
    end

    for i = 1:numdetects
        %Grab the covariance of each eigenvector
        Kprimefull = Dfull(i,i)*Vfull(:,i)*Vfull(:,i)';
        
        %Compare SNR's calculated from Eigenvalues/vectors
        SNReig(i,q+1) = 10*log10(trace(Kprimefull)/(m*n*sigma));


        %Beamform it
        [BF_Out_full] = Make_Beamformer(Kprimefull,SourceFreq,Arraymat,'B',phi,theta);

        
        BF_Out_All(:,:,q+1,i) = BF_Out_full;
 
        %Find angles of maximum power
        [rfull,cfull] = find(BF_Out_full == max(max(BF_Out_full)));

        %Store the locations
        DetectionLocationsFull(i,:,q+1) = [phi(cfull(1)),theta(rfull(1))];
    end
end

function [Vtransformed,Vnexttransformed,ambiguity] = MaximizeCorrelation(Vpair,alpharange,anglerange,m,n)
    v = Vpair(:,1);
    vnext = Vpair(:,2);
    ambiguity = zeros(length(alpharange),length(anglerange));
    for i = 1:length(alpharange)
        alpha = alpharange(i);
        for j = 1:length(anglerange)
            angle = anglerange(j);
            Vtransformed = alpha*exp(1i*angle)*v + sqrt(1 - alpha^2)*vnext;
            Vnexttransformed = -sqrt(1-alpha^2)*exp(1i*angle)*v + alpha*vnext;
            [~,ambiguity(i,j)] = CalcMeanVertPairCorrelation(Vtransformed,m,n);
        end
    end
    [alphaind,angleind] = find(ambiguity == max(max(ambiguity)));
    angle = anglerange(angleind(1));
    alpha = alpharange(alphaind(1));
    Vtransformed = alpha*exp(1i*angle)*v + sqrt(1 - alpha^2)*vnext;
    Vnexttransformed = -sqrt(1-alpha^2)*exp(1i*angle)*v + alpha*vnext;
end

function [correlationvector,meancorrelation] = CalcMeanVertPairCorrelation(V,m,n)
    correlationvector = zeros(n-1,1);
    for x = 1:n-1
        [VertSubIndex,AdjVertSubIndex] = GetSubarrayIndex(m,n,x);
        v = V(VertSubIndex);
        vnext = V(AdjVertSubIndex);
        v = v/norm(v);
        vnext = vnext/norm(vnext);
        correlationvector(x) = abs(v'*vnext);
    end
    meancorrelation = mean(correlationvector);
end