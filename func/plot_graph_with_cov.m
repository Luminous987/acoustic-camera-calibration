function plot_graph_with_cov(g, iteration, H)

if nargin<2
    iteration = -1;
end

clf;
hold on;
plot3(nan, nan, nan, 'LineStyle','none','Marker','s', 'MarkerSize', 4,'MarkerEdgeColor','g', 'LineWidth',2);
plot3(nan, nan, nan, 'Color','g');
plot3(nan, nan, nan, 'LineStyle','none','Marker','o', 'MarkerSize', 4,'MarkerEdgeColor','r', 'LineWidth',2);

% 去掉了表示声源估计值的浅蓝色点和相应的图例
% plot3(nan, nan, nan, 'LineStyle','none','Marker','s', 'MarkerSize', 4,'MarkerEdgeColor','c'); 

plot3(nan, nan, nan, 'LineStyle','none','Marker','x', 'MarkerSize', 4,'MarkerEdgeColor','b');

% 修改后的图例
legend('Mic. pos. est.','Sigma region of mic. pos. est.','Mic. pos. g. t.','Sound source g. t.');

[p, l] = get_poses_landmarks(g);

if (length(l) > 0)
    landmarkIdxX = l+1;
    landmarkIdxY = l+2;
    landmarkIdxZ = l+3;
    plot3(g.x(landmarkIdxX), g.x(landmarkIdxY), g.x(landmarkIdxZ), 'LineStyle','none','Marker','s', 'MarkerSize', 4,'MarkerEdgeColor','g', 'LineWidth',2);
    plot3(g.x_gt(landmarkIdxX), g.x_gt(landmarkIdxY), g.x_gt(landmarkIdxZ), 'LineStyle','none','Marker','o', 'MarkerSize', 4,'MarkerEdgeColor','r', 'LineWidth',2);
end

if (length(p) > 0)
    pIdxX = p+1;
    pIdxY = p+2;
    pIdxZ = p+3;
    % plot3(g.x(pIdxX), g.x(pIdxY), g.x(pIdxZ), 'LineStyle','none','Marker','s', 'MarkerSize', 4,'MarkerEdgeColor','c');
    % 去掉了浅蓝色的声源估计值
    % plot3(g.x_gt(pIdxX), g.x_gt(pIdxY), g.x_gt(pIdxZ),'Color','b', 'LineStyle','-','Marker','x', 'MarkerSize', 4,'MarkerEdgeColor','b');
    plot3(g.x_gt(pIdxX), g.x_gt(pIdxY), g.x_gt(pIdxZ),'Color','b', 'LineStyle','-', 'Marker','x', 'MarkerSize', 4, 'MarkerEdgeColor','b');
end

if iteration>1
    P = inv(H);
    for m=2:g.M
        mic_m = g.x(3*(m-1)+1:3*(m-1)+3);
        mic_C = full(P(3*(m-1)+1:3*(m-1)+3,3*(m-1)+1:3*(m-1)+3));
        
        if m==g.M_x
            mic_C(2,2) = 0;
            mic_C(2,1) = 0;
            mic_C(1,2) = 0;
            mic_C(3,3) = 0;
            mic_C(3,1:2) = zeros(1,2);
            mic_C(1:2,3) = zeros(2,1);
        end
        
        if m==g.M_x*(g.M_y-1)+1
            mic_C(3,3) = 0;
            mic_C(3,1:2) = zeros(1,2);
            mic_C(1:2,3) = zeros(2,1);
        end
        
        drawprobellipse_3d(mic_m,mic_C,0.683,'g'); % 0.683 0.954 0.997
    end
end

grid on;
view(30,15);
hold off;

gcf;
grid on; axis equal;
xlabel('X (m)');ylabel('Y (m)');zlabel('Z (m)');

% 更新后的图例
legend('Mic. pos. est.','Sigma region of mic. pos. est.','Mic. pos. g. t.','Sound source g. t.');
drawnow;

end
