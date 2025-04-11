
%
%   Coded by Yi Jin

function [azi_rot, ele_rot] = rotAOAs_G2L(antLocs, pointLoc)
% azi_rot, ele_rot degree

origin_local = antLocs(3, :); % origin point of local coordinate system

% the local x-axis direction vector (antenna 3 to antenna 4 direction)
x_local = antLocs(4, :) - antLocs(3, :);
x_local = x_local / norm(x_local); % 归一化

% the local y-axis direction vector (antenna 2 to antenna 3 direction)
y_local = antLocs(3, :) - antLocs(2, :);
y_local = y_local / norm(y_local); % 归一化

% the local z-axis direction vector
z_local = cross(x_local, y_local);
z_local = z_local / norm(z_local); % 归一化

% the rotation matrix from the global coordinate system to the local coordinate system
R_global_to_local = [x_local', y_local', z_local'];
point_local = R_global_to_local * (pointLoc' - origin_local');

x_local_coord = point_local(1);
y_local_coord = point_local(2);
z_local_coord = point_local(3);

azi_rot = atan2(y_local_coord, x_local_coord) * (180 / pi);
ele_rot = atan2( sqrt(x_local_coord^2 + y_local_coord^2), z_local_coord ) * (180 / pi);

end