
% get_Signal_by_IQ: Preprocesses BLE IQ data for phased array processing
%   Extracts and aligns raw IQ samples from multiple packets, performs phase
%   compensation and frequency offset estimation, then constructs 
%   signal matrices for direction finding applications
%
%   Key Features:
%   - Handles BLE protocol packet structure
%   - Compensates phase discontinuities across packets
%   - Estimates carrier frequency offset via linear regression
%   - Supports dual-polarized uniform planar array (UPA) processing
%
%   Coded by Yi Jin

%   Inputs:
%   - IQ_list [cell array]
%       Raw BLE packets containing complex IQ samples and metadata
%   - pack_idx [vector]
%       Selected packet indices for processing
%
%   Outputs:
%   - X_vert_all [matrix]
%       Vertical polarization signal matrix 
%   - X_hori_all [matrix]
%       Horizontal polarization signal matrix 
%   - X_polar [matrix]
%       Combined polarization signal matrix 
%   - T_allpack [integer]
%       Total observation samples across all packets
%   - RSSI [scalar]
%       Average Received Signal Strength Indicator

function [X_vert_all, X_hori_all, X_polar, T_allpack, RSSI] = get_Signal_by_IQ( IQ_list, pack_idx )

len_onepack = 164;   % Length of samples in one packet
% len_halfpack = 42;   % Original length of phase data
len_halfpack = 82;    % Adjusted length of phase data
num_pack = length(pack_idx);   % Number of packets to process
simple_flag = 1;

% Stack raw IQ samples
IQ_sample = zeros(len_onepack*num_pack ,1);
RSSI_pack = zeros(1, num_pack);
for paci = 1:num_pack
    pidx = pack_idx(paci);
    pack_struct = IQ_list{pidx};
    IQ_sample( (paci-1)*len_onepack+1:paci*len_onepack ) = pack_struct.samples';
    RSSI_pack(paci) = pack_struct.rssi; % RSSI values
end
RSSI = mean(RSSI_pack); 

% Calculate phase from IQ components
I_seq = IQ_sample(1:2:end);
Q_seq = IQ_sample(2:2:end);
IQ_seq = [I_seq,Q_seq];
phase_seq = atan2(-Q_seq, I_seq);

% phase compensation
phase_inc = phase_seq;
for paci = 1:num_pack
    cnt = 0;
    for i = 2: len_halfpack
        idx_seq = (paci-1)*len_halfpack+i;
        if( phase_seq(idx_seq)<phase_seq(idx_seq-1) )
            cnt = cnt + 1;
        end
        phase_inc(idx_seq) = phase_seq(idx_seq) + cnt *2*pi;  
    end
end
% Estimate carrier frequency deviation via least squares fitting
len_ref = 8; % Reference cycle length (typically 8)
f_dev_pack = zeros(num_pack,1);
for paci = 1:num_pack
    phase_x = (paci-1)*len_halfpack+1 : (paci-1)*len_halfpack+8;
    p = polyfit( phase_x, phase_inc(phase_x),  1 );
    f_dev_pack(paci) = p(1)/(2*pi);
    % figure();
    % plot( phase_x, phase_inc(phase_x), 'o--','LineWidth',1.2,'MarkerSize',5 ); hold on;
    % disp( strcat( 'p:', num2str(f_dev_pack(paci))) )
end
% f_dev = mean(f_dev_pack);

% figure( 'color','w');
% plot( phase_inc, 'o-','LineWidth',lw,'MarkerSize',ms );
% hold on;
% 
% xlabel( 'Sample time/($\mu$s)','interpreter','latex' );
% ylabel( 'Phase','interpreter','latex' );
% grid on; xlim([1,82]);

M = 16; % Number of antenna elements
M_DP = 32; % Dual-polarization mode element number
num_cycles = floor( ( len_halfpack - len_ref )/M_DP ); % Cycles per packet
%
T = M_DP*num_cycles;
T_simp = num_cycles;
T_allpack = T*num_pack;
X_all = zeros( M_DP, T_allpack );
X_polar = zeros( M_DP, T_allpack );
% Construct data matrix  
phase_sw_idx = [9:40; 42:73];
%  Separate vertical and horizontal polarization components 
for paci = 1:num_pack
    phase_inc_one = phase_inc( (paci-1)*len_halfpack+1:paci*len_halfpack );
    for cyci = 1:num_cycles
        idx_cyci = phase_sw_idx(cyci,:);
        Phase_inc_sw( 1:M_DP, (cyci-1)*M_DP+1:cyci*M_DP ) = diag( phase_inc_one( idx_cyci ) ); % 原：(cyci-1)*M_DP+start_loc: cyci*M_DP+start_loc-1
        % Cycle block processing
        Amp_sw(:,cyci) = sqrt( IQ_seq( idx_cyci,1).^2 + IQ_seq( idx_cyci,2).^2)';
        for row_i = 1:M_DP
            for col_i = (cyci-1)*M_DP+1:cyci*M_DP
                blc_row = (cyci-1)*M_DP+row_i;
                Phase_inc_sw(row_i,col_i) = (Phase_inc_sw(row_i,blc_row)+  2*2*pi*f_dev_pack(paci)*(col_i-blc_row)); %c_cb *2*pi*f_dev*(col_i-blc_row) );
            end
        end
        X_sw( :, (cyci-1)*M_DP+1:cyci*M_DP ) =  Amp_sw(:,cyci) .* exp( 1j*( Phase_inc_sw(:, (cyci-1)*M_DP+1:cyci*M_DP) ) );
        
    end
    X_all(:, (paci-1)*T+1:paci*T ) = X_sw;
end

% Organize polarization components
X_vert_all = X_all(1:2:M_DP, :);  % Vertical polarization
X_hori_all = X_all(2:2:M_DP, :);   % Horizontal polarization
X_polar = [ X_hori_all; X_vert_all ];

if simple_flag == 1
    T_allpack = T_simp*num_pack;
    X_polar = X_polar(:, [1,33]);
end


end
