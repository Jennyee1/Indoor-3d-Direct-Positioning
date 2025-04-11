

% WAG_DPD_3D 3D: Weighted fusion with Adaptive Grid refinement for Direct Positioning Determination
%   Implements an enhanced 3D direct positioning method with adaptive grid refinement
%   for indoor localization scenarios

%   Recommended Parameters:
%   (1) Breakroom: min_resolution = [0.26, 0.28, 0.121]
%   (2) Anteroom: min_resolution = [0.33, 0.28, 0.081]

%   Coded by Yi Jin

%   Inputs:
%   - baseStationsLocs [m]
%       Matrix of size (numBS)x3 containing 3D coordinates of base stations
%   - AntLocs [m]
%       Cell array of antenna locations for each base station
%   - receivedSignals
%       Cell array containing received signals from each base station
%   - Signal_basic [struct]
%       Structure containing basic signal parameters:
%       - Mx: Number of antenna elements in x-direction
%       - My: Number of antenna elements in y-direction
%       - lambda: Wavelength [m]
%       - d: Antenna spacing [m]
%       - K: Number of signal sources
%   - bandwidth [Hz]
%       Signal bandwidth
%   - oversamplingFactor
%       Oversampling ratio for signal processing
%   - diamaterSearchArea [m]
%       Search area dimensions [x, y]
%   - resolution [m]
%       Initial grid resolution [x, y, z]
%   - numIter
%       Number of grid refinement iterations
%   - x_bound [m], y_bound [m]
%       Search area boundaries

%   Output:
%   - mobileUserLocEstimate [m]
%       1x3 vector with estimated 3D coordinates of target


function [ mobileUserLocEstimate ] = WAG_DPD_3D( baseStationsLocs, AntLocs, receivedSignals, Signal_basic, bandwidth, oversamplingFactor, diamaterSearchArea, resolution, numIter, x_bound, y_bound )
% 推荐参数 
% (1) Tearoom
% min_resolution = [0.26, 0.26, 0.121];
% (2) Anterroom
min_resolution = [0.33, 0.28, 0.081];  % 各轴最小分辨率 [x,y,z]  min_resolution = [0.2, 0.2, 0.05];

% -----------------------------------------------------------
%   Initialization
% -----------------------------------------------------------
% Configuration parameters
sigma_base = 0.5;                   % Base standard deviation
prob_threshold = 0.2;               % Probability density threshold

numBS = size(baseStationsLocs,1);
[M, numSamples] = size(receivedSignals{1});
numSensors = M*ones(numBS, 1);
Mx = Signal_basic.Mx;
My = Signal_basic.My;
lambda = Signal_basic.lambda;
d = Signal_basic.d;
K = 1;

%% -----------------------------------------------------------
%   Main Iteration Loop
% -----------------------------------------------------------
prev_estimates = [];
for i=1:numIter
    
    % -------------------------------------------------------
    %   Dynamic Resolution Adjustment
    % -------------------------------------------------------
    if i > 1
        % Covariance-based resolution adaptation
        Sigma = cov(prev_estimates); 
        sigma_axis = sqrt(diag(Sigma));
        sigma_axis(sigma_axis<1e-3) = sigma_axis(sigma_axis<1e-3)+1e-3;
        
        % Adaptive resolution calculation
        % resolution = min([resolution.*(sigma_axis'/sigma_base); min_resolution], [], 2)';
        resolution = min([resolution.*(sigma_axis'/sigma_base); min_resolution]);

        % Non-uniform grid generation
        [gridLocs, x_now, y_now, z_now] = refineGrid3D_enhanced(...
            resolution, prev_estimates, prob_threshold, diamaterSearchArea, x_bound, y_bound );
    else
        % ---------------------------------------------------
        %   Initial Grid Generation (Sparse Z-axis)
        % ---------------------------------------------------
        [x,y,z] = meshgrid(...
            linspace(0, diamaterSearchArea(1), ceil(diamaterSearchArea(1)/resolution(1))),...
            linspace(0, diamaterSearchArea(2), ceil(diamaterSearchArea(2)/resolution(2))),...
            0.2:resolution(3):3);  % Initial z-axis resolution: 0.4m
        x_now = x + x_bound; 
        y_now = y + y_bound; 
        z_now = z;
        gridLocs = [x_now(:), y_now(:), z_now(:)];
    end

    % -------------------------------------------------------
    %   Spatial Likelihood Calculation
    % -------------------------------------------------------
    numLocs = size(gridLocs,1);
    likelihood_BS = zeros(numBS,numLocs);

    % Calculate 3D distances
    temp = repmat(permute(gridLocs,[3,2,1]),[numBS,1,1]) - repmat(baseStationsLocs,[1,1,numLocs]);
    distBS2locs = sqrt(squeeze(sum(temp.^2,2)));
    
    % 
    offsetSamples = ceil(3*oversamplingFactor/pi*sqrt(2*log(2)));

    for l=1:numBS
        % 
        atom = oversamplingFactor^(-1/2)*(pi/log(2))^(1/4)*...
            exp(-((0:numSamples-1).'*ones(1,numLocs)./oversamplingFactor -...
            bandwidth*ones(numSamples,1)*distBS2locs(l,:)./3e8 -...
            offsetSamples/oversamplingFactor).^2.*pi.^2./2./log(2));
        spatialData = receivedSignals{l} *conj(atom);

        A_Steering = zeros(M, numLocs);
        for v = 1:numLocs
            [AziBS2locs_rot, EleBS2locs_rot] = rotAOAs_G2L(AntLocs{l}, gridLocs(v,:));

            ay = exp( -1j*2*pi*d*(0:My-1)'*sin(AziBS2locs_rot*pi/180)*sin(EleBS2locs_rot*pi/180) /lambda );
            ax = exp( 1j*2*pi*d*(0:Mx-1)'*cos(AziBS2locs_rot*pi/180)*sin(EleBS2locs_rot*pi/180) /lambda  );
            aa = kron(ay,ax) ;
            A_Steering(:,v) = aa;

        end
        
        likelihood_BS(l,:) = abs(sum(conj(A_Steering).*spatialData,1)).^2/numSensors(l);
    end

    % -------------------------------------------------------
    %   Probability Density Estimation
    % -------------------------------------------------------
    if i>1
        BS_weight = vecnorm(baseStationsLocs(:,1:2)-mean(prev_estimates(1:10,1:2)),2,2);
        weight = 1./(BS_weight.^2);
        weight = weight./sum(weight);
        likelihood_weight = sum( likelihood_BS.* weight );
        likelihood = likelihood_weight;
    else
        likelihood = sum(likelihood_BS);
    end
    likelihood = likelihood / max(likelihood);  % Normalization and candidate selection
    prob_density = likelihood / sum(likelihood);

    % -------------------------------------------------------
    %   Candidate Point Selection
    % -------------------------------------------------------
    [~, idx] = sort(prob_density, 'descend');
    valid_points = gridLocs(idx(1:ceil(numLocs*prob_threshold)), :);
    prev_estimates = valid_points(1:min(50,end), :);  % Keep top 50 candidates

end

% Final position estimation
mobileUserLocEstimate = mean(prev_estimates(1:3,:));  % Average top 3 candidates
end

%% -----------------------------------------------------------
%   Function: refineGrid3D_enhanced
% -----------------------------------------------------------
% refineGrid3D_enhanced: Generates adaptive 3D search grid
%   Implements dynamic grid refinement based on previous estimates
function [activeLocs, gridx, gridy, gridz] = refineGrid3D_enhanced(resolution, prev_estimates, prob_thresh, diamaterSearchArea, x_bound, y_bound)

    mu = mean(prev_estimates);
    sigma = std(prev_estimates);
    
    regrids = ceil(3*sigma ./ resolution);  
    regrids = max(regrids, [3,3,2]);  % Minimum expansion
    regrids = min(regrids, [15,15,6]);  % Maximum expansion 
    
    x_range = mu(1) + (-regrids(1):regrids(1)) * resolution(1);
    y_range = mu(2) + (-regrids(2):regrids(2)) * resolution(2);
    z_range = mu(3) + (-regrids(3):regrids(3)) * resolution(3);
    % z_range = mu(3) + (0:regrids(3)) * resolution(3);
    
    [X,Y,Z] = meshgrid(x_range, y_range, z_range);
    activeLocs = [X(:), Y(:), Z(:)];
    
    % Apply search area constraints
    x_min = max( (min(prev_estimates(:,1)) - 3*resolution(1)), x_bound-0.25 );
    x_max = min( (max(prev_estimates(:,1)) + 3*resolution(1)), x_bound+diamaterSearchArea(1)+0.25 );
    y_min = max( (min(prev_estimates(:,2)) - 3*resolution(2)), y_bound-0.25 );
    y_max = min( (max(prev_estimates(:,2)) + 3*resolution(2)), y_bound+diamaterSearchArea(2)+0.25 );
    z_min = 0.2; %min(prev_estimates(:,3)) - 2*resolution(2);
    z_max = 4; %max(prev_estimates(:,3)) + 2*resolution(2);

    % Filter valid grid points
    valid_idx = (activeLocs(:,1)>=x_min) & (activeLocs(:,1)<=x_max) & ...
                (activeLocs(:,2)>=y_min) & (activeLocs(:,2)<=y_max) &...
                (activeLocs(:,3)>=z_min) & (activeLocs(:,3)<=z_max);
    activeLocs = activeLocs(valid_idx,:);
    
    gridx = unique(activeLocs(:,1));
    gridy = unique(activeLocs(:,2));
    gridz = unique(activeLocs(:,3));

end
