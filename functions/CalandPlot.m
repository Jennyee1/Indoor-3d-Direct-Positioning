function [] = CalandPlot(vertices, sta, dim)
% % 示例顶点坐标 (顺时针或逆时针顺序)
% vertices = [
%     1, 1;  % 顶点1
%     4, 1;  % 顶点2
%     4, 4;  % 顶点3
%     1, 4   % 顶点4
% ];
if strcmp( dim, '2D' )
    text(vertices(1,1)-0.2, vertices(1,2)+0.2,num2str(sta));
    % 计算各边与x轴正方向的夹角
    % angles = zeros(4, 1);
    % for i = 1:4
    %     p1 = vertices(i, :);
    %     p2 = vertices(mod(i, 4) + 1, :);  % 下一个顶点 (循环)
    %     vector = p2 - p1;
    %     angle = atan2(vector(2), vector(1));  % 计算夹角
    %     angles(i) = rad2deg(angle);  % 转换为度数
    % end
    
    % % 显示各边的夹角
    % disp('各边与x轴正方向的夹角（度数）：');
    % disp( mean(angles([1,3])) );
    
    % fprintf( 'Station %d, 各边与x轴正方向的夹角（度数）:  %f, %f, %f, %f \n\n', sta, angles   );
    
    % 绘制正方形 StaX数组中, 1-Ant16, 2-Ant13, 3-Ant1, 4-Ant4
    hold on;
    for i = 1:4
        p1 = vertices(i, :);
        p2 = vertices(mod(i, 4) + 1, :);  % 下一个顶点 (循环)
        plot([p1(1), p2(1)], [p1(2), p2(2)], 'b-', 'LineWidth', 2);  % 画边
        % if i == 1
        %     text(p1(1), p1(2), '16');
        % elseif i == 2
        %     text(p1(1), p1(2), '13');
        % elseif i == 3
        %     text(p1(1), p1(2), '1');
        % else
        %     text(p1(1), p1(2), '4');
        % end
    end
    % hold off;
    
    % 计算四边形中心点
    center = mean(vertices, 1);
    
    % % 添加小坐标系的 x 和 y 方向箭头
    % quiver(center(1), center(2), 0.5, 0, 'Color', '#e17055', 'LineWidth', 1, 'MaxHeadSize', 0.5);  % x 方向箭头
    % quiver(center(1), center(2), 0, -0.5, 'Color', '#0984e3', 'LineWidth', 1, 'MaxHeadSize', 0.5);  % y 方向箭头
    % % 标注角度值
    % text(center(1) + 0.55, center(2)-0.15, '0°', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', '#e17055','FontSize',10);
    % text(center(1)-0.15, center(2) - 0.8, '90°', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', '#0984e3','FontSize',10);
    % 
    axis equal;
    xlabel('X');
    ylabel('Y');
    % title('Stations');
end
if strcmp( dim, '3D' )
    text(vertices(1,1)-0.2, vertices(1,2)+0.2, vertices(1,3), num2str(sta));
    % 计算各边与x轴正方向的夹角
    % angles = zeros(4, 1);
    % for i = 1:4
    %     p1 = vertices(i, :);
    %     p2 = vertices(mod(i, 4) + 1, :);  % 下一个顶点 (循环)
    %     vector = p2 - p1;
    %     angle = atan2(vector(2), vector(1));  % 计算夹角
    %     angles(i) = rad2deg(angle);  % 转换为度数
    % end
    
    % % 显示各边的夹角
    % disp('各边与x轴正方向的夹角（度数）：');
    % disp( mean(angles([1,3])) );
    
    % fprintf( 'Station %d, 各边与x轴正方向的夹角（度数）:  %f, %f, %f, %f \n\n', sta, angles   );
    
    % 绘制正方形 StaX数组中, 1-Ant16, 2-Ant13, 3-Ant1, 4-Ant4
    hold on;
    for i = 1:4
        p1 = vertices(i, :);
        p2 = vertices(mod(i, 4) + 1, :);  % 下一个顶点 (循环)
        plot3([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], 'b-', 'LineWidth', 2);  % 画边
        % if i == 1
        %     text(p1(1), p1(2), '16');
        % elseif i == 2
        %     text(p1(1), p1(2), '13');
        % elseif i == 3
        %     text(p1(1), p1(2), '1');
        % else
        %     text(p1(1), p1(2), '4');
        % end
    end
    % hold off;
    
    % 计算四边形中心点
    center = mean(vertices, 1);
    
    % 添加小坐标系的 x 和 y 方向箭头
    quiver3(center(1), center(2), center(3), 0.5, 0, 0, 'Color', '#e17055', 'LineWidth', 1, 'MaxHeadSize', 0.5);  % x 方向箭头
    quiver3(center(1), center(2), center(3), 0, -0.5, 0, 'Color', '#0984e3', 'LineWidth', 1, 'MaxHeadSize', 0.5);  % y 方向箭头
    % 标注角度值
    text(center(1) + 0.55, center(2)-0.15, center(3), '0°', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', '#e17055','FontSize',10);
    text(center(1)-0.15, center(2) - 0.8, center(3), '90°', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', '#0984e3','FontSize',10);
    
    axis equal;
    xlabel('X');
    ylabel('Y');
    % title('Stations');
end


end