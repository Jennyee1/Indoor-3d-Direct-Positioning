
% -----------------------------------------------------------
%   Experiment Setup in break room
% -----------------------------------------------------------
%  Code for paper: "An Enhanced 3-D Direct Positioning Method 
%                  with Adaptive Grid Refinement in Indoor Environments"
%  Experiment configuration:
%  - 500 packets collected
%  - 10 example test points in breakroom scenario
%  - Requires JSON data loading support

clc; 
clear;
addpath('./data/BreakRoom/');
addpath('./settings');
addpath('./functions');

%% ----------------------------------------------------------
%   Get basic information
% -----------------------------------------------------------
PlotFlag = 0;
PlotFlag_loc = 0;
run('TrueValuePlot_3D.m');

%% ----------------------------------------------------------
%   Data Processing Parameters
% -----------------------------------------------------------
pack_data_index = 1:500;    % Selected packet indices for processing
len_dataL = length(pack_data_index);
Azi = -180:2:180;         % Azimuth search grid (degrees)
Ele = 5:2:75;             % Elevation search grid (degrees)

%% ----------------------------------------------------------
%   System Configuration Parameters
% -----------------------------------------------------------
% Physical constants
c = physconst('LightSpeed');  % Speed of light (m/s)

% Array configuration
signal_basic.d = 0.036;       % Element spacing (3.6cm)
carrierFreq = 2.4e9;          % Carrier frequency [Hz]
signal_basic.fc = carrierFreq;  
signal_basic.lambda = c/signal_basic.fc;  % Wavelength calculation

% Algorithm parameters
signal_basic.K = 1;           % Number of signal sources
signal_basic.Mx = 4;          % UPA x-dimension elements
signal_basic.My = 4;          % UPA y-dimension elements
startCellSize = 0.61;         % Initial grid resolution (m)
oversamplingFactor = 6;       % Spectrum oversampling ratio
bandwidth = 1e6;              % Signal bandwidth [Hz]
mapWidth = [6.51, 6.51];      % Spatial search boundaries (m)
x_bound_DPD = 11.01;          % X-axis search limit (m)
y_bound_DPD = 0.01;           % Y-axis search threshold (m)
signal_basic.azi = Azi;
signal_basic.ele = Ele;

%% ----------------------------------------------------------
%   Result Containers Initialization  
% -----------------------------------------------------------
MULocEstim_WDPD = zeros(len_dataL,3);  % 3D position estimates
SelectedPots = 1:10;                   % Selected test points
numPots = length(SelectedPots);

% Error metric containers
Dis_3D_WDPD = zeros(1, numPots*len_dataL);  % 3D distance errors
Dis_2D_WDPD = zeros(1, numPots*len_dataL);  % Horizontal errors  
H_WDPD = zeros(1, numPots*len_dataL);       % Height errors


%% ----------------------------------------------------------
%   Main Processing Loop
% -----------------------------------------------------------

for pot_n = 1:numPots
    potn = SelectedPots(pot_n);
    potstr = ['P', num2str(potn)];      % Test point ID
    Pot = Points_Tearoom(potn,:);       % Ground truth coordinates
    % fprintf(['\n ', potstr,', oF = ', num2str(oversamplingFactor),'，  ', num2str(Pot(1:2)), '\n']);

    % File pattern matching for current test point
    [P_AllFiles, P_stationIndices] = findfiles('./data/BreakRoom', ID_Stations, potstr);

    %% ------------------------------------------------------
    %   Station Configuration
    % -------------------------------------------------------
    P_stationIndices = [12,13,14,15,16];  % Selected base stations
    len_Stations = length(P_stationIndices);
    Station_pos = zeros(len_Stations,3);  % BS positions
    for ssi = 1:len_Stations
        tmp = eval(['Sta',num2str(P_stationIndices(ssi))]);
        Station_pos(ssi, :) = tmp(3,:);
    end

    %% ------------------------------------------------------
    %   Signal Processing Pipeline
    % -------------------------------------------------------
    % Signal Extraction from IQ data and AOA Estimation

    IQ_azi_MUSIC_seq = zeros(len_dataL, max(P_stationIndices));
    IQ_ele_MUSIC_seq = zeros(len_dataL, max(P_stationIndices));

    AntLocs = cell(len_Stations,1);
    receivedSignals = cell(len_Stations,1);
    antLocs = cell(len_Stations,1);
    for n = 1:len_Stations

        stan_temp = P_stationIndices(n);
        station_temp = eval(['Sta', num2str(P_stationIndices(n))]);
        AntLocs{n} = station_temp;

        IQ_list_cell{stan_temp} = loadjson(P_AllFiles{P_stationIndices(n)-11});

        [IQ_azi_MUSIC_seq(:,stan_temp), IQ_ele_MUSIC_seq(:,stan_temp),  RSSI_seq(:,stan_temp), receivedsignals{n}] = get_Signal_AOAs_by_IQ(...
            IQ_list_cell{stan_temp}, pack_data_index, signal_basic );

        RSSI(stan_temp) = mean( RSSI_seq(:,stan_temp) );

    end

    %% ------------------------------------------------------
    %   Positioning Estimation Method
    % -------------------------------------------------------
    % The proposed 3-D WAG-DPD method
    time_WDPD = 0;

    for tidx = 1:len_dataL
        for n = 1:len_Stations
            receivedSignals{n} = (receivedsignals{n}{tidx});
        end
        Tstart = cputime;
        resolution_init = [1.0, 1.0, 0.6];  
        MULocEstim_WDPD_tmp = WAG_DPD_3D( Station_pos, AntLocs, receivedSignals, signal_basic, bandwidth, oversamplingFactor, mapWidth, resolution_init, 3, x_bound_DPD, y_bound_DPD );
        time_WDPD = time_WDPD + cputime - Tstart;
        MULocEstim_WDPD(tidx,:) = MULocEstim_WDPD_tmp;
        Dis_3D_WDPD(1,(pot_n-1)*len_dataL+tidx) = norm( MULocEstim_WDPD(tidx,:)- Pot(1:3) );
        Dis_2D_WDPD(1,(pot_n-1)*len_dataL+tidx) = norm( MULocEstim_WDPD(tidx,1:2)- Pot(1:2) );
        H_WDPD(1,(pot_n-1)*len_dataL+tidx) = norm( MULocEstim_WDPD(tidx,3)- Pot(3) );
        
    end

    Time_WDPD(pot_n) = time_WDPD/len_dataL;  % Average time per estimation

end

%% ----------------------------------------------------------
%   Result Visualization
% -----------------------------------------------------------

figure('Color','w'); hold on;
Method_vec = { 'Dis_2D_WDPD';'H_WDPD'; 'Dis_3D_WDPD'; }; % 
Disp_vec = { 'WAG-DPD (Horizontal Distance)'; 'WAG-DPD (Vertical Distance)';  'WAG-DPD (3-D Distance)'; }; % 
colorvec = [{'#b7d332'};{'#8a1874'};{'#1264a6'};{'#008665'};{'#e69f37'};{'#d15730'}];

% Generate CDF curves
for i = 1:length(Method_vec)
    data =  eval(Method_vec{i}); 
    sorted_data = sort(data);
    cdf_values = (1:length(data)) / length(data);
    plot(sorted_data, cdf_values, 'LineWidth', 2, 'LineStyle', '-', ... 
         'Marker', 'o', 'MarkerSize', 4, ... 
         'MarkerIndices', round(linspace(1, length(data), 10)), ... 
         'color', colorvec{i}, 'DisplayName', Disp_vec{i});
    
    % Calculate 67% line
    x_at_067 = interp1(cdf_values, sorted_data, 0.67, 'linear');    
    disp([' When 67%, the distance error is ', num2str(x_at_067)]);
end
yline(0.67, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);  % 灰色虚线
text(0.85, 0.71, '67%', 'Color', [0.5 0.5 0.5], 'FontSize', 10, 'FontWeight', 'bold');  % 在y=0.67位置上方标注

xlabel('Distance (m)', 'FontSize', 12, 'FontWeight', 'bold');  
ylabel('Cumulative Probability', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0, 1]); 
ylim([0, 1]);  
title('CDF Curves', 'FontSize', 14, 'FontWeight', 'bold');
grid on; box on;
legend('show', 'Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 10, 'LineWidth', 1.2);  




