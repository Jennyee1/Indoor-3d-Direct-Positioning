
% find_2D_peaks: Detects top K peaks in 2D power spectrum
%   Identifies dominant peaks in 2D spatial spectrum and returns corresponding
%   azimuth-elevation coordinates
%
%   Coded by Yi Jin

%   Inputs:
%   - P [matrix]
%       2D power spectrum matrix
%   - Azi [vector]
%       Azimuth grid vector corresponding to rows of P
%   - Ele [vector]
%       Elevation grid vector corresponding to columns of P
%   - K [integer]
%       Number of dominant peaks to identify

%   Outputs:
%   - azi_est [vector]
%       Estimated azimuth angles (K strongest peaks)
%   - ele_est [vector]
%       Estimated elevation angles (K strongest peaks)

function [azi_est, ele_est] = find_2D_peaks( P,Azi,Ele, K )
PeaksMap = imregionalmax(P);
PeaksValue = PeaksMap .* P;

[row, col] = find( PeaksValue~=0 );
len_row = length(row);
peak_vec = zeros(len_row, 1);
% 
for ri = 1:len_row
    peak_vec(ri) = P( row(ri),col(ri) );
end
[peak_val, in]=sort(peak_vec, 'descend');
if K<length(in)
    ind_est = in(1:K);
    azi_est = Azi( row(ind_est) );
    ele_est = Ele( col(ind_est) );
else
    ind_est = in;
    azi_est = Azi( row(ind_est) );
    ele_est = Ele( col(ind_est) );
end

end
