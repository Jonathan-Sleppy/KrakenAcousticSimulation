%%
clear
sphericalspreading = 0;
BellhopRun = 0; 
KrakenRun = 1;
fr = 100;
c = 1500;
k = 2*pi*fr/c;
lambda = c./fr;
dz = lambda/2;
zd = 50 + [0:1:19]'.*dz;
xd = [0:1:19]'.*dz - 19*dz/2;xd = flipud(xd);
lz = length(zd);
lx = length(xd);
arraycoords = [];
for n = 1:length(xd);
    arraycoords = [arraycoords;  xd(n)*ones(length(zd),1) zeros(length(zd),1) zd ];
end;
Nph = size(arraycoords,1);
s1s = [-750 8000 160];
s1e = [-750 8000 160];
t = linspace(0,1,10);
s1 = [linspace(s1s(1),s1e(1),length(t))' linspace(s1s(2),s1e(2),length(t))' linspace(s1s(3),s1e(3),length(t))'];
s2s = [-250 8000 110];
s2e = [-250 8000 110];
s2 = [linspace(s2s(1),s2e(1),length(t))' linspace(s2s(2),s2e(2),length(t))' linspace(s2s(3),s2e(3),length(t))'];
s3s = [250 8000 500];
s3e = [250 8000 500];
s3 = [linspace(s3s(1),s3e(1),length(t))' linspace(s3s(2),s3e(2),length(t))' linspace(s3s(3),s3e(3),length(t))'];
s4s = [750 8000 150];
s4e = [750 8000 150];
s4 = [linspace(s4s(1),s4e(1),length(t))' linspace(s4s(2),s4e(2),length(t))' linspace(s4s(3),s4e(3),length(t))'];
s5s = [-1200 12000 220];
s5e = [1200 12000 220];
s5 = [linspace(s5s(1),s5e(1),length(t))' linspace(s5s(2),s5e(2),length(t))' linspace(s5s(3),s5e(3),length(t))'];
for m = 1:length(t)
    for n = 1:size(arraycoords,1);
        R1(n,m) = sqrt((s1(m,1) - arraycoords(n,1)).^2 + (s1(m,2) - arraycoords(n,2)).^2 +(s1(m,3) - arraycoords(n,3)).^2);
        R2(n,m) = sqrt((s2(m,1) - arraycoords(n,1)).^2 + (s2(m,2) - arraycoords(n,2)).^2 +(s2(m,3) - arraycoords(n,3)).^2);
        R3(n,m) = sqrt((s3(m,1) - arraycoords(n,1)).^2 + (s3(m,2) - arraycoords(n,2)).^2 +(s3(m,3) - arraycoords(n,3)).^2);
        R4(n,m) = sqrt((s4(m,1) - arraycoords(n,1)).^2 + (s4(m,2) - arraycoords(n,2)).^2 +(s4(m,3) - arraycoords(n,3)).^2);
        R5(n,m) = sqrt((s5(m,1) - arraycoords(n,1)).^2 + (s5(m,2) - arraycoords(n,2)).^2 +(s5(m,3) - arraycoords(n,3)).^2);
        R1c(n,m) = sqrt((s1(m,1) - arraycoords(n,1)).^2 + (s1(m,2) - arraycoords(n,2)).^2);
        R2c(n,m) = sqrt((s2(m,1) - arraycoords(n,1)).^2 + (s2(m,2) - arraycoords(n,2)).^2);
        R3c(n,m) = sqrt((s3(m,1) - arraycoords(n,1)).^2 + (s3(m,2) - arraycoords(n,2)).^2);
        R4c(n,m) = sqrt((s4(m,1) - arraycoords(n,1)).^2 + (s4(m,2) - arraycoords(n,2)).^2);
        R5c(n,m) = sqrt((s5(m,1) - arraycoords(n,1)).^2 + (s5(m,2) - arraycoords(n,2)).^2);
    end;
    b1(m) = atan2(s1(m,2),s1(m,1))*180/pi;
    b2(m) = atan2(s2(m,2),s2(m,1))*180/pi;
    b3(m) = atan2(s3(m,2),s3(m,1))*180/pi;
    b4(m) = atan2(s4(m,2),s4(m,1))*180/pi;
    b5(m) = atan2(s5(m,2),s5(m,1))*180/pi;
end;

sigma1 = 5e-8;
P1 = zeros(Nph,length(t));
P2 = zeros(Nph,length(t));
P3 = zeros(Nph,length(t));
P4 = zeros(Nph,length(t));
P5 = zeros(Nph,length(t));
K = zeros(Nph,Nph,length(t));
A1 = 0.6;
A2 = 1;
A3 = 0.4;
A4 = 0.75;
A5 = .35;
for m = 1:length(t)
    if(sphericalspreading == 1)
        P1(:,m) = A1.*exp(1i*k*R1(:,m))./R1(:,m);
        P2(:,m) = A2.*exp(1i*k*R2(:,m))./R2(:,m);
        P3(:,m) = A3.*exp(1i*k*R3(:,m))./R3(:,m);
        P4(:,m) = A4.*exp(1i*k*R4(:,m))./R4(:,m);
        P5(:,m) = A5.*exp(1i*k*R5(:,m))./R5(:,m);
    elseif(BellhopRun == 1);
        [pr,rr,rd] = RunBellhop(R1c(1:lz:end,m),s1(m,3),zd,fr);
        pr1 = A1.*pr.';
        [pr,rr,rd] = RunBellhop(R2c(1:lz:end,m),s2(m,3),zd,fr);
        pr2 = A2.*pr.';
        [pr,rr,rd] = RunBellhop(R3c(1:lz:end,m),s3(m,3),zd,fr);
        pr3 = A3.*pr.';
        [pr,rr,rd] = RunBellhop(R4c(1:lz:end,m),s4(m,3),zd,fr);
        pr4 = A4.*pr.';
        [pr,rr,rd] = RunBellhop(R5c(1:lz:end,m),s5(m,3),zd,fr);
        pr5 = A5.*pr.';
    elseif(KrakenRun == 1);
        [pr,rr,rd] = RunKraken(R1c(1:lz:end,m),s1(m,3),zd,fr);
        pr1 = A1.*pr.';
        [pr,rr,rd] = RunKraken(R2c(1:lz:end,m),s2(m,3),zd,fr);
        pr2 = A2.*pr.';
        [pr,rr,rd] = RunKraken(R3c(1:lz:end,m),s3(m,3),zd,fr);
        pr3 = A3.*pr.';
        [pr,rr,rd] = RunKraken(R4c(1:lz:end,m),s4(m,3),zd,fr);
        pr4 = A4.*pr.';
        [pr,rr,rd] = RunKraken(R5c(1:lz:end,m),s5(m,3),zd,fr);
        pr5 = A5.*pr.';
    end;
    tmp1 =[];tmp2=[];tmp3=[];tmp4=[];tmp5=[];
    for n = 1:length(xd)%stack into billboard array
        tmp1 = [tmp1;pr1(n,:).'];
        tmp2 = [tmp2;pr2(n,:).'];
        tmp3 = [tmp3;pr3(n,:).'];
        tmp4 = [tmp4;pr4(n,:).'];
        tmp5 = [tmp5;pr5(n,:).'];
    end;
    P1(:,m) = tmp1;
    P2(:,m) = tmp2;
    P3(:,m) = tmp3;
    P4(:,m) = tmp4;
    P5(:,m) = tmp5;
    %add uncorrelated noise
    K(:,:,m) = P1(:,m)*P1(:,m)' + P2(:,m)*P2(:,m)' + P3(:,m)*P3(:,m)' + P4(:,m)*P4(:,m)' + P5(:,m)*P5(:,m)';
    K(:,:,m) = K(:,:,m) + sigma1.*eye(Nph);
    SNR1(m) = 10*log10(trace(P1(:,m)*P1(:,m)')./(Nph*sigma1));
    SNR2(m) = 10*log10(trace(P2(:,m)*P2(:,m)')./(Nph*sigma1));
    SNR3(m) = 10*log10(trace(P3(:,m)*P3(:,m)')./(Nph*sigma1));
    SNR4(m) = 10*log10(trace(P4(:,m)*P4(:,m)')./(Nph*sigma1));
    SNR5(m) = 10*log10(trace(P5(:,m)*P5(:,m)')./(Nph*sigma1));
    m
end;
%Beamform
if(0)
    phi = [0:180];
    theta = [-30:2:30];
    bf_all = zeros(length(t),length(theta),length(phi));
    bf_allM = zeros(length(t),length(theta),length(phi));
    for m = 1:length(t)
        [bf_out] =  Make_Beamformer(squeeze(K(:,:,m)),fr,arraycoords,'M',phi,theta);
        bf_allM(m,:,:) = bf_out;
        [bf_out] =  Make_Beamformer(squeeze(K(:,:,m)),fr,arraycoords,'B',phi,theta);
        bf_all(m,:,:) = bf_out;
        m
    end;
    % clf
    for n = 1:length(theta);
        figure(1)
        set(gcf, 'PaperPosition', [.25, .25, 10, 10]);
        tmp = squeeze(bf_allM(:,n,:)); tmp = tmp./max(max(tmp));
        pcolor(t,phi,10*log10(tmp)');shading('flat')
        caxis([-3,0])
        axis([0,1,75,105])
        title('Example C MVDR')
        xlabel('Relative Time')
        ylabel('Bearing (deg)')
        set(gca,'Fontsize',14)
        colorbar
        drawnow;
        pause
    end;
    for n = 1:length(theta);
        %     print -dpng ExC_Kraken_MVDR_BearingTime_Theta0.png
        figure(2)
        tmp = squeeze(bf_all(:,n,:));tmp = tmp./max(max(tmp));
        set(gcf, 'PaperPosition', [.25, .25, 10, 10]);
        pcolor(t,phi,10*log10(tmp)');shading('flat')
        caxis([-3,0])
        axis([0,1,75,105])
        title('Example C Conv. BF')
        xlabel('Relative Time')
        ylabel('Bearing (deg)')
        set(gca,'Fontsize',14)
        colorbar
        %     print -dpng ExC_Kraken_ConvBF_BearingTime_Theta0.png
        drawnow;
        n
        pause
    end
end;

%Eigen Processing
phi = [-180:.01:180];
theta = [0];
Numeig = 5;
eigbf = zeros(length(t),length(phi),length(Nph-Numeig+1:Nph));
pk = zeros(Numeig,length(t));
for m = 1:length(t)
    [V,D] = eigs(squeeze(K(:,:,m)));
    for n = 1:5;%Nph-4:Nph
        K1 = V(1:lz:end,n)*V(1:lz:end,n)';
        [bf_out] =  Make_Beamformer(K1,fr,arraycoords(1:lz:end,:),'B',phi,theta);
        eigbf(m,:,n) = bf_out;
        [I,J] = max(bf_out);
        pk(n,m) = abs(phi(J));
        %correlation
        tmp = 0;
        for mm = 1:lx-1
            t1 = V((mm - 1)*lz + 1:1:mm*lz,n);
            t2 = V(mm*lz + 1:1:(mm + 1)*lz,n);
            tmp = tmp + abs(t1'*t2./sqrt((t1'*t1)*(t2'*t2)));
        end
        gamma(n,m) = tmp./(lx - 1);
    end;
    m
end
figure(1)
clf
set(gcf, 'PaperPosition', [.25, .25, 10, 10]);
hold on;
plot(t,pk(1,:),'r','Linewidth',2);plot(t,pk(2,:),'g','Linewidth',2);plot(t,pk(3,:),'b','Linewidth',2); ...
    plot(t,pk(4,:),'Color',[0.9290 0.6940 0.1250],'Linewidth',2); ...
    plot(t,pk(5,:),'Color',[0.4940 0.1840 0.5560],'Linewidth',2);
axis([0,1,75,105])
title('Example C MVB')
xlabel('Relative Time')
ylabel('Bearing (deg)')
set(gca,'Fontsize',14)
legend('Eig 1','Eig 2','Eig 3','Eig 4','Eig 5')
hold off
% print -dpng MVB_BillboardArray_400elements_EigsFull_Kraken.png
figure(2)
clf
set(gcf, 'PaperPosition', [.25, .25, 10, 10]);
hold on;
plot(t,SNR1,'r','Linewidth',2);plot(t,SNR2,'g','Linewidth',2);plot(t,SNR3,'b','Linewidth',2); ...
    plot(t,SNR4,'Color',[0.9290 0.6940 0.1250],'Linewidth',2); ...
    plot(t,SNR5,'Color',[0.4940 0.1840 0.5560],'Linewidth',2);
axis([0,1,-10,10])
title('Example C SNR')
xlabel('Relative Time')
ylabel('Signal-to-Noise Ratio (dB)')
set(gca,'Fontsize',14)
legend('Eig 1','Eig 2','Eig 3','Eig 4','Eig 5')
hold off
% print -dpng SNR_ExampleC.png
%correlations
figure(3)
clf
set(gcf, 'PaperPosition', [.25, .25, 10, 10]);
hold on
plot(t,gamma(1,:),'r','Linewidth',2);plot(t,gamma(2,:),'g','Linewidth',2);plot(t,gamma(3,:),'b','Linewidth',2); ...
    plot(t,gamma(4,:),'Color',[0.9290 0.6940 0.1250],'Linewidth',2); ...
    plot(t,gamma(5,:),'Color',[0.4940 0.1840 0.5560],'Linewidth',2);
axis([0,1,.5,1.1])
title('Example C Vertical Stave Correlations')
xlabel('Relative Time')
ylabel('Correlation')
set(gca,'Fontsize',14)
legend('Eig 1','Eig 2','Eig 3','Eig 4','Eig 5')
hold off
figure(4)
clf
set(gcf, 'PaperPosition', [.25, .25, 10, 10]);
hold on;
plot(t,b1(1,:),'r','Linewidth',2);plot(t,b2,'g','Linewidth',2);plot(t,b3,'b','Linewidth',2); ...
    plot(t,b4,'Color',[0.9290 0.6940 0.1250],'Linewidth',2); ...
    plot(t,b5,'Color',[0.4940 0.1840 0.5560],'Linewidth',2);
axis([0,1,75,105])
title('Example C MVB')
xlabel('Relative Time')
ylabel('Bearing (deg)')
set(gca,'Fontsize',14)
legend('Eig 1','Eig 2','Eig 3','Eig 4','Eig 5')
hold off


%%
