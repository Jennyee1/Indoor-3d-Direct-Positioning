
% get_Signal_AOAs_by_IQ: Processes IQ data to extract array signals and AOAs
%   Obtain the signal matrices
%   Generates MUSIC spectrum and estimates dominant AoA (Azimuth/Elevation)
%   for multiple data packets using subspace-based MUSIC algorithm
%
%   Note: WAG-DPD framework bypasses explicit AoA estimation step
%
%   Coded by Yi Jin

%   Inputs:
%   - IQ_list [cell array]
%       Raw IQ data packets with complex samples
%   - pack_data_idx [vector]
%       Packet indices for batch processing
%   - signal_basic [struct]
%       System parameters containing:
%       azi/ele - Angular grids 
%       fc - Carrier frequency (Hz)
%       d - Antenna spacing (m)
%       ......

%   Outputs:
%   - IQ_azi_MUSIC_seq [matrix]
%       MUSIC azimuth estimates 
%   - IQ_ele_MUSIC_seq [matrix]
%       MUSIC elevation estimates
%   - rssi_seq [vector]
%       RSSI values per packet
%   - X_sw_seq [cell array]
%       signal matrices

% -----------------------------------------------------
function [IQ_azi_MUSIC_seq, IQ_ele_MUSIC_seq,  rssi_seq, X_sw_seq] = get_Signal_AOAs_by_IQ( IQ_list, pack_data_idx, signal_basic )
K = 1;
c = physconst('LightSpeed');
mx = 4;
my = 4;

len_dataL = length(pack_data_idx);
IQ_azi_MUSIC_seq = zeros(len_dataL,K);
IQ_ele_MUSIC_seq = zeros(len_dataL,K);
rssi_seq = zeros(len_dataL,K);

Power = cell(1, len_dataL);
X_sw_seq = cell(1, len_dataL);
for p_idx = 1:len_dataL
    pack_idx = pack_data_idx(p_idx);
    %% Estimation AOA

    signal_basic.K = K;
    signal_basic.Mx = mx;
    signal_basic.My = my;
    
    signal_basic.fc = 2.44e9;
    signal_basic.d = 0.036;    %3.6cm
    signal_basic.lambda = c/signal_basic.fc;

    Azi = signal_basic.azi;
    Ele = signal_basic.ele;

    [~, X_sw, X_polar, T_all, rssi_seq(p_idx)] = get_Signal_by_IQ( IQ_list, pack_idx );
    
    signal_basic.rssi_seq = rssi_seq;
    signal_basic.T = T_all;
    signal_basic.X = X_polar;
    R_sw = X_sw*X_sw'/T_all;
    Power{p_idx} = zeros(length(Azi), length(Ele));
    X_sw_seq{p_idx} = X_sw;
    
    [P_MUSIC] = cal_AOAs_MUSIC( R_sw, signal_basic);
    [IQ_azi_MUSIC_seq(p_idx,:), IQ_ele_MUSIC_seq(p_idx,:)] = find_2D_peaks( P_MUSIC, Azi, Ele, K );


end
end