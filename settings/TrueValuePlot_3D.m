
%   Coded by Yi Jin

run('BS_TrueValues.m');
run('Points_TrueValue.m');

% Overall Scene for plot
Scene = [
8.414, 9.44;
8.414, 13.300;
19.9040, 13.300;
19.9040, 1.600;
17.8040, 1.600;
17.8040, 0;
10.414, 0;
10.414, 1.55;
9.6640, 1.55;
9.6640, 6.34;
9.6640+8.2, 6.34;
17.8640, 6.34+3.1;
17.8640-9.45, 9.44;
];
Scene(:,3) = 0;
Mirror = [
17.8640, 6.34;
17.8640, 3;
];
Mirror(:,3) = 0;

if PlotFlag
    figure('Color','w');  hold on;

    % Anterroom anchors
    CalandPlot(Sta3,3,'3D');
    CalandPlot(Sta4,4,'3D');
    CalandPlot(Sta5,5,'3D');
    CalandPlot(Sta6,6,'3D');

    % Breakroom anchors
    CalandPlot(Sta12,12,'3D');
    CalandPlot(Sta13,13,'3D');
    CalandPlot(Sta14,14,'3D');
    CalandPlot(Sta15,15,'3D');
    CalandPlot(Sta16,16,'3D');
    hold on;
    plot3(Points_Anterroom(:, 1), Points_Anterroom(:,2), Points_Anterroom(:,3), 'r*');
    plot3(Points_Tearoom(1:10, 1), Points_Tearoom(1:10,2), Points_Tearoom(1:10,3), 'r*');

    plot3(Scene(:,1), Scene(:,2), Scene(:,3), '-'); 
    
    % xlim([8,21.5]);
    % ylim([-0.5, 14]);
    view(3);
    grid on;
    axis equal;
    box on;
    
end

return
