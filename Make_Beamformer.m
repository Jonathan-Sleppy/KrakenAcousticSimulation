function [bf_out] =  Make_Beamformer(cov,fr,arraycoords,bftype,phi,theta);
%
%MAKE_BEAMFORMER    Make the beamformed output from data covariance matrix.
%                   
%
%INPUT        
%           cov  ==> Cross spectral density matrix (chan, chan, freqs)
%            fr  ==> Frequencies (Hz)
%   arraycoords  ==> Vector with phones locations (starting at 0) in m.
%        bftype  ==> Beamforming type (1 for CBF, don't recommend using ABF for this)                                                   
%          wind  ==> 'H' for hanning 'T' for Taylorwin (anything else is boxcar)
%
%
%OUTPUT
%         theta  ==> Angles for beamformer 
%         bf_out  ==> Beamformer output
%          
%
%USAGE

%M. Siderius Portland State/NRL 7/26/2022
%
bf_out = zeros(length(theta),length(phi),length(fr));

%bftype = 1 for plane wave beamforming 0 = adaptive;
if(exist('bftype') == 0)
    bftype = 'B';
end;
if(bftype == 'B')
    %disp('Bartlett')
    for ii = 1:length(fr)
        k0 = 2*pi*fr(ii)/1500;
        for kk = 1:length(theta)
            for mm = 1:length(phi)
                steervec = [cosd(theta(kk)).*cosd(phi(mm)); cosd(theta(kk))*sind(phi(mm)); sind(theta(kk))];
                ph_dist_wv = k0.*arraycoords*steervec;
                e = exp(-1i*ph_dist_wv); e = e(:);
                covtmp = cov(:,:,ii);
                bf_out(kk,mm,ii) = e'*covtmp*e;
            end;
        end;
    end;
elseif(bftype == 'M')
    disp('MVDR')
    for ii = 1:length(fr)
        cvinv = inv(cov(:,:,ii));
        k0 = 2*pi*fr(ii)/1500;
         for kk = 1:length(theta)
            for mm = 1:length(phi)
                steervec = [cosd(theta(kk)).*cosd(phi(mm)); cosd(theta(kk))*sind(phi(mm)); sind(theta(kk))];
                ph_dist_wv = k0.*arraycoords*steervec;
                e = exp(-1i*ph_dist_wv); e = e(:);
                covtmp = cov(:,:,ii);
                bf_out(kk,mm,ii) = 1./(eps + e'*cvinv*e);
            end;
        end;
    end;
end;
bf_out = abs(bf_out);

