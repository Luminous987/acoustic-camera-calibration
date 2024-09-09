function g = modify_and_save_g_real(output_filename,u,measurement,outlier)
%%
g = struct();
%% 对麦克风阵列的位置进行修改，对应g.x_gt的前40行
    g.M = 8;
    % 麦克风间距
    % mic_dis = 0.5;
    %     % 正方体顶点的相对坐标
    % mic_positions = [
    %     0.0, -0.0, 0; % 第一个麦克风在原点
    %     mic_dis , 0, 0;
    %     0, mic_dis , 0;
    %     0, 0, mic_dis ;
    %     mic_dis , mic_dis , 0;
    %     mic_dis , 0, mic_dis ;
    %     0 , mic_dis , mic_dis ;
    %     mic_dis , mic_dis , mic_dis 
    % ];
    % 
    % % % 正方体顶点的相对坐标
    % mic_positions = [
    %     0, 0, 0; % 第一个麦克风在原点
    %     0.37, 0, 0;
    %     0, 0.39, 0;
    %     -0.01, -0.01, 0.39;
    %     0.36, 0.39, 0;
    %     0.38, 0, 0.39;
    %     -0.02, 0.39, 0.39;
    %     0.38, 0.38, 0.39
    % ];

    % 正方体顶点在相机坐标系下的坐标
    mic_positions = [
        -0.10, 0, 1.16; % 第一个麦克风在原点
        0.27, 0, 1.16;
        -0.10, 0, 1.55;
        -0.10, -0.39, 1.16;
        0.27, 0, 1.55;
        0.27, -0.39, 1.16;
        -0.10, -0.39, 1.55;
        0.27, -0.39, 1.55
    ];
    % 设置麦克风的位置信息
    for i = 1:g.M
        base_row = (i-1)*3 + 1;
        g.x_gt(base_row:base_row+2, 1) = mic_positions(i, :)';
    end
% +[-0.18;-0.28;-0.24]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%导入声源位置真实数据,并删去离群点
for i = 1:length(outlier)
    % 计算出要删除的索引范围
    idx_start = (outlier(i) - 1) * 6 + 1;  % 每组有 6 个点
    idx_end = idx_start + 5;  
    
    measurement(idx_start:idx_end, :) = [];
    
    u_start = (outlier(i) - 1) * 18 + 1;  % 每组6个点，对应18个坐标
    u_end = u_start + 17;  
    u(u_start:u_end) = [];
end
num_points = (length(u)/3);
g.x_gt(1+3*g.M:3*g.M+length(u),1) = u;

% figure;
for i = 1:7   
    measured_values =  measurement(:,i);
    % plot(measured_values, 'b-o', 'DisplayName', '测量值');
    % hold on
end

%% 下面生成TDOA测量量
    c = 340;              % 声速 (m/s)
    TDOA_std_dev_noise = 0.0666e-3; % 噪声的标准偏差,对应文章的R
    % TDOA_std_dev_noise = 1e-3 * c * 0.01;
    % TDOA_std_dev_noise = 0.007; %单位m
    % 声源位置存储在g.x_gt的41到280行中，每三个行对应一个(x, y, z)坐标
    % 计算TDOA测量值并添加到g.edge.measurement中

    % TDOA_theoretical_value_sum = 0;
    % outlier = [];
    for k = 1:num_points
        % 获取当前声源的位置
        source_idx = 3 * g.M + (k - 1) * 3;
        source_pos = g.x_gt(source_idx + 1:source_idx + 3);
    
        % 计算第一个麦克风的位置
        first_mic_pos = g.x_gt(1:3);
    
        % 计算声源到第一个麦克风的距离
        dist_first_mic = norm(source_pos - first_mic_pos);
    
        % 初始化TDOAs和POSEs数组
        TDOAs = zeros(g.M-1, 1);
    
        % 遍历其余的麦克风，计算TDOA
        for i = 2:g.M
            mic_idx = (i - 1) * 3;
            mic_pos = g.x_gt(mic_idx + 1:mic_idx + 3);
            dist_mic = norm(source_pos - mic_pos);
            TDOA_theoretical_value = (dist_mic - dist_first_mic) / c; 
            % TDOAs(i - 1) = (dist_mic - dist_first_mic) / c ;
            % % TDOAs(i - 1) = (dist_mic - dist_first_mic);
            
            TDOAs(i - 1) = measurement(k, i-1);
            % if abs(measurement(k, i-1) - TDOA_theoretical_value) > 2e-04 
            %     outlier = [outlier; k];
            % end
            % TDOA_theoretical_value_sum = TDOA_theoretical_value_sum + measurement(k, i-1) - TDOA_theoretical_value;
        end
        % 加入测量噪声
        % TDOAs_noise = TDOA_std_dev_noise * randn(g.M-1, 1);
        % TDOAs = TDOAs + TDOAs_noise * 1;
        
        % 计算信息矩阵
        TDOA_var_noise = TDOA_std_dev_noise^2;      
        % 信息矩阵是方差倒数的对角矩阵
        TDOA_information = diag(repmat(1 / TDOA_var_noise, g.M-1, 1));      

        % 保存TDOA测量值到g.edge的measurement字段
        g.edges(k).measurement = TDOAs;
        g.edges(k).information = TDOA_information;
        g.edges(k).fromIdx = 1;
        g.edges(k).toIdx = 3 * g.M + 3 * k - 2;
    end

% TDOA_theoretical_value_sum = TDOA_theoretical_value_sum / (num_points * 7);
% disp('理论TDOA差值');
% disp(TDOA_theoretical_value_sum);
% disp('离群点')
% disp(outlier);
%% 给g.x一个初值

    % 初始化gx
    init_std_dev = 0.5;
    
    % 初始化g.x与g.x_gt的维度相同
    g.x = g.x_gt;
    
    % 为g.x赋值，是g.x_gt的值加上高斯噪声
    g.x(1:24,1) = g.x_gt(1:24,1) + init_std_dev * randn(24,1) ;
    %% 填充idLookup    
    % 填充麦克风的idLookup
    for i = 1:g.M
        g.idLookup(i).offset = (i - 1) * 3;
        g.idLookup(i).dimension = 3;
    end
    
    % 填充声源的idLookup
    for i = 1:num_points
        g.idLookup(g.M + i).offset = 3 * g.M + (i - 1) * 3 ;
        g.idLookup(g.M + i).dimension = 3;
    end
    g.M_x = 2;
    g.M_y = 2;

    % 保存修改后的结构体
    save(output_filename, 'g');

end
