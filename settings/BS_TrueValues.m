% The known positions of multiple anchors 
% measured by a total station accurately

%   Coded by Yi Jin

% Determine the relationship between the anchor (or base station) index
% and the respective ID
ID_Stations = [
'ble-pd-0C4314F45CD8'; % anchor 3
'ble-pd-0C4314F46D44'; % anchor 4
'ble-pd-E0798D2F3247'; % anchor 5
'ble-pd-0C4314F4686A'; % anchor 6
'ble-pd-0C4314F46B12'; % anchor 12
'ble-pd-0C4314F46C94'; % anchor 13
'ble-pd-0C4314F46BC8'; % anchor 14
'ble-pd-60B64763790D'; % anchor 15
'ble-pd-0C4314F46C77'; % anchor 16
];

%% Anterroom anchors 3-6
% the positions of anchors is acquired by the total station;
% When building the global coordinate system based on the total station, 
% x-axis inverses as a benchmark, and the coordinates need to be reversed symmetrically.
%           X       Y       Z
Ref3 =      [13.589,  10.82, 0];

Sta6_ori =  [15.199,  9.999, 2.740; ...
             15.308,  10.004, 2.740; ...
             15.312,  9.893, 2.739; ...
             15.200,  9.888, 2.738; ...
             ];
Sta6 = Sta6_ori;
Sta6(:,1) = Ref3(1) - (Sta6_ori(:,1)-Ref3(1));
Sta6_center = mean(Sta6);

Sta5_ori =  [13.708,  12.570, 2.738; ...
             13.819,  12.574, 2.738; ...
             13.820,  12.462, 2.736; ...
             13.711,  12.461, 2.737; ...
             ];
Sta5 = Sta5_ori;
Sta5(:,1) = Ref3(1) - (Sta5_ori(:,1)-Ref3(1));
Sta5_center = mean(Sta5);

Sta4_ori =  [11.423,  9.987, 2.740; ...
             11.534,  9.988, 2.741; ...
             11.537,  9.878, 2.739; ...
             11.425,  9.875, 2.738; ...
             ];
Sta4 = Sta4_ori;
Sta4(:,1) = Ref3(1) - (Sta4_ori(:,1)-Ref3(1));
Sta4_center = mean(Sta4);

Sta3_ori =  [9.745,  12.602, 2.736; ...
             9.852,  12.603, 2.738; ...
             9.851,  12.494, 2.737; ...
             9.742,  12.493, 2.737; ...
             ];
Sta3 = Sta3_ori;
Sta3(:,1) = Ref3(1) - (Sta3_ori(:,1)-Ref3(1));
Sta3_center = mean(Sta3);


%% Breakroom anchors 12-16
% y-axis inverses
Ref2 = [18.269-0.016, 2.547, 0];
Sta12_ori =  [16.436-0.016,  1.490,  2.727; ...
             16.328-0.016,   1.488,  2.726; ...
             16.326-0.016,   1.599,  2.725; ...
             16.435-0.016,   1.601,  2.726; ...
             ];
Sta12 = Sta12_ori;
Sta12(:,2) = Ref2(2) - (Sta12_ori(:,2)-Ref2(2));
Sta12_center = mean(Sta12);

Sta14_ori = [15.034-0.016, -0.163, 2.724;...
            14.923-0.016, -0.163, 2.723;...
            14.923-0.016, -0.052, 2.724;...
            15.031-0.016, -0.052, 2.724];
Sta14 = Sta14_ori;

Sta15_ori = [12.135-0.016, -0.191, 2.722;...
            12.025-0.016, -0.188, 2.722;...
            12.028-0.016, -0.078, 2.722;...
            12.136-0.016, -0.081, 2.722]; 
Sta15 = Sta15_ori;


Sta16_ori = [15.054-0.016, 3.348, 2.732;...
            14.942-0.016, 3.347, 2.731;...
            14.941-0.016, 3.458, 2.732;...
            15.054-0.016, 3.459, 2.732]; 
Sta16 = Sta16_ori;
Sta13_ori = [13.678-0.016, 1.462, 2.728;...
            13.570-0.016, 1.459, 2.727;...
            13.566-0.016, 1.571, 2.727;...
            13.677-0.016, 1.572, 2.727]; 
Sta13 = Sta13_ori;

Sta13(:,2) = Ref2(2) - (Sta13_ori(:,2)-Ref2(2));
Sta13_center = mean(Sta13);
Sta14(:,2) = Ref2(2) - (Sta14_ori(:,2)-Ref2(2));
Sta14_center = mean(Sta14);
Sta15(:,2) = Ref2(2) - (Sta15_ori(:,2)-Ref2(2));
Sta15_center = mean(Sta15);
Sta16(:,2) = Ref2(2) - (Sta16_ori(:,2)-Ref2(2));
Sta16_center = mean(Sta16);


%% rotation angle
X_local = zeros(  17,3 );
Y_local = zeros(  17,3 );
Z_local = zeros(  17,3 );
Alpha = zeros(  17,1 );
Beta = zeros( 17,1 );
Gamma = zeros( 17,1 );
for sta = [3:6,12:16]
    [X_local(sta,:), Y_local(sta,:), Z_local(sta,:), Alpha(sta), Beta(sta), Gamma(sta)] = GetRotation( eval(['Sta',num2str(sta)]) );
end

return


function [x_local, y_local, z_local, alpha, beta, gamma] = GetRotation( station )
% x_local, y_local, z_local : 3x1 vector
% the global coordinate of the x/y/z-axis vector in the local coordinate system
x_local = station(4, :) - station(3, :); 
x_local = x_local / norm(x_local); 
y_local = station(3, :) - station(2, :);
y_local = y_local / norm(y_local); 
z_local = cross(x_local, y_local); 
z_local = z_local / norm(z_local); 
origin = station(3, :); % origin point of the local coordinate system

R_local_to_global = [x_local; y_local; z_local]'; % rotation matrix (from local coordinates to global coordinates)
% Euler angles (rotation around the x, y, z axes) using rotm2eul 
% using 'XYZ' rotation order 
% (i.e. first around the x axis, then around the y axis, and finally around the z axis)
eul_angles = rotm2eul(R_local_to_global, 'XYZ');

alpha = rad2deg(eul_angles(1)); % rotation angle about x axis
beta = rad2deg(eul_angles(2));  % rotation angle about y axis
gamma = rad2deg(eul_angles(3)); % rotation angle about z axis

end



