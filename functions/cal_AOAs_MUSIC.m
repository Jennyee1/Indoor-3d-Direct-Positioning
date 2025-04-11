
% cal_AOAs_MUSIC: Computes MUSIC spectrum for 2D AoA estimation
%   Implements the Multiple Signal Classification (MUSIC) algorithm for
%   joint azimuth-elevation estimation using uniform planar arrays (UPA)
%
%   Inputs:
%   - R [matrix]
%       Sample covariance matrix (Mx*My x Mx*My)
%   - Signal_basic [struct]
%       System configuration parameters containing:
%       Mx/My - UPA dimensions in x/y axes
%       K - Number of signal sources
%       d - Inter-element spacing (m)
%       lambda - Wavelength (m)
%       azi/ele - Angular grids 
%       ......
%
%   Output:
%   - p_music [matrix]
%       MUSIC pseudo-spectrum (azimuth x elevation)
%


function [p_music] = cal_AOAs_MUSIC( R, Signal_basic)
Mx = Signal_basic.Mx;
My = Signal_basic.My;
K = Signal_basic.K;
d = Signal_basic.d;
lambda = Signal_basic.lambda;

azi = Signal_basic.azi;
ele = Signal_basic.ele;

[EV,D] = eig(R);
[EVA,I] = sort( diag(D).' );
EV = fliplr( EV(:,I) );
Un = EV(:,K+1:end);

p_music = zeros(length(azi),1);
for ii = 1:length(azi)
    for jj = 1:length(ele)

        ay = exp( -1j*2*pi*d*(0:My-1)'*sin(azi(ii)*pi/180)*sin(ele(jj)*pi/180) /lambda );
        ax = exp( 1j*2*pi*d*(0:Mx-1)'*cos(azi(ii)*pi/180)*sin(ele(jj)*pi/180) /lambda  );
        
        aa = kron(ay,ax)  ;
        Power = det(aa'*(Un*Un')*aa);
        
        p_music(ii,jj)= 1/abs(Power);
    end
end


end


