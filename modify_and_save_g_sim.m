function g = modify_and_save_g_sim(output_filename)
%%
g = struct();
%% 对麦克风阵列的位置进行修改，对应g.x_gt的前40行
    g.M = 8;
    % 麦克风间距
    mic_dis = 0.5;
        % 正方体顶点的相对坐标
    mic_positions = [
        0, 0, 0; % 第一个麦克风在原点
        mic_dis , 0, 0;
        0, mic_dis , 0;
        0, 0, mic_dis ;
        mic_dis , mic_dis , 0;
        mic_dis , 0, mic_dis ;
        0 , mic_dis , mic_dis ;
        mic_dis , mic_dis , mic_dis 
    ];

    % % 正方体顶点的相对坐标
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
    % 设置麦克风的位置信息
    for i = 1:g.M
        base_row = (i-1)*3 + 1;
        g.x_gt(base_row:base_row+2, 1) = mic_positions(i, :)';
    end

% 下面生成声源位置
    x = -0.7;
    y = -1.8;
    z = 2.3;
    delta_x = 0.5;
    delta_y = 0.5;
    delta_z = 0.0;
    bias_z = 1.35;
    % x = -5;
    % y = -20;
    % z = 20;
    % delta_x = 5;
    % delta_y = 5;
    % delta_z = 0.5;
    % 菱形顶点的相对坐标
    mic_positions = [
        x, y, 0 - bias_z; 
        x + delta_x, y, 0 - bias_z;
        x + delta_x * 2, y + delta_y, 0 - bias_z;
        x + delta_x, y + delta_y, 0 - bias_z;
        x, y, 0 - bias_z;
        x, y, z - bias_z;
        x + delta_x, y, z - bias_z;
        x + delta_x * 2, y + delta_y, z - bias_z;
        x + delta_x, y + delta_y, z - bias_z;
        x, y, z - bias_z;
        x + delta_x * 2, y + delta_y, z - bias_z;
        x + delta_x * 2, y + delta_y, 0 - bias_z;
        x, y, 0 - bias_z;

        x, y + 6 * delta_y, 0 + delta_z - bias_z; 
        x + delta_x, y + 6 * delta_y, 0 + delta_z * 2 - bias_z;
        x + delta_x * 2, y + delta_y + 6 * delta_y, 0 + delta_z * 3 - bias_z;
        x + delta_x, y + delta_y + 6 * delta_y, 0 + delta_z * 4 - bias_z;
        x, y + 6 * delta_y, 0 + delta_z * 5 - bias_z;
        x, y + 6 * delta_y, z + delta_z * 6 - bias_z;
        x + delta_x, y + 6 * delta_y, z + delta_z * 7 - bias_z;
        x + delta_x * 2, y + delta_y + 6 * delta_y, z + delta_z * 8 - bias_z;
        x + delta_x, y + delta_y + 6 * delta_y, z + delta_z * 9 - bias_z;
        x, y + 6 * delta_y, z + delta_z * 10 - bias_z;
        x + delta_x * 2, y + delta_y + 6 * delta_y, z + delta_z * 11 - bias_z;
        x + delta_x * 2, y + delta_y + 6 * delta_y, 0 + delta_z * 12 - bias_z;
        x, y + 6 * delta_y, 0 + delta_z * 13 - bias_z;

        
        x + 2, y + 2, 0 - bias_z; 
        x + 2, y + 2 + delta_x, 0 - bias_z;
        x + 2 + delta_y, y + 2 + delta_x * 2, 0 - bias_z;
        x + 2 + delta_y, y + 2 + delta_x, 0 - bias_z;
        x + 2, y + 2, 0 - bias_z;
        x + 2, y + 2, z - bias_z;
        x + 2, y + 2 + delta_x, z - bias_z;
        x + 2 + delta_y, y + 2 + delta_x * 2, z - bias_z;
        x + 2 + delta_y, y + 2 + delta_x, z - bias_z;
        x + 2, y + 2, z - bias_z;
        x + 2 + delta_y, y + 2 + delta_x * 2, z - bias_z;
        x + 2 + delta_y, y + 2 + delta_x * 2, 0 - bias_z;
        x + 2, y + 2, 0 - bias_z;

        % x + 2 - 6 * delta_y, y + 2, 0 + delta_z - bias_z; 
        % x + 2 - 6 * delta_y, y + 2 + delta_x, 0 + delta_z * 2 - bias_z;
        % x + 2 - 6 * delta_y - delta_y, y + 2 + delta_x * 2, 0 + delta_z * 3 - bias_z;
        % x + 2 - 6 * delta_y - delta_y, y + 2 + delta_x, 0 + delta_z * 4 - bias_z;
        % x + 2 - 6 * delta_y, y + 2, 0 + delta_z * 5 - bias_z;
        % x + 2 - 6 * delta_y, y + 2, z + delta_z * 6 - bias_z;
        % x + 2 - 6 * delta_y, y + 2 + delta_x, z + delta_z * 7 - bias_z;
        % x + 2 - 6 * delta_y - delta_y, y + 2 + delta_x * 2, z + delta_z * 8 - bias_z;
        % x + 2 - 6 * delta_y - delta_y, y + 2 + delta_x, z + delta_z * 9 - bias_z;
        % x + 2 - 6 * delta_y, y + 2, z + delta_z * 10 - bias_z;
        % x + 2 - 6 * delta_y - delta_y, y + 2 + delta_x * 2, z + delta_z * 11 - bias_z;
        % x + 2 - 6 * delta_y - delta_y, y + 2 + delta_x * 2, 0 + delta_z * 12 - bias_z;
        % x + 2 - 6 * delta_y, y + 2, 0 + delta_z * 13 - bias_z;
    ];

    num_samples_per_segment = 6;
    num_points = (length(mic_positions)-1) * num_samples_per_segment; % 采样点个数

    % 初始化声源位置存储变量
    g.x_gt(3*g.M+1:3*g.M+1+num_points*3-1) = 0; % 乘以3是因为每个时间点有三个坐标值
    current_index = 24;

    % 计算每个采样点的位置
    for i = 1:(length(mic_positions)-1) 
        for j = 1:num_samples_per_segment %每段num_samples_per_segment采样点
            t = (j - 1) / (num_samples_per_segment-1); 
            x_coord = mic_positions(i,1) + t * (mic_positions(i+1,1) - mic_positions(i,1));
            y_coord = mic_positions(i,2) + t * (mic_positions(i+1,2) - mic_positions(i,2));
            z_coord = mic_positions(i,3) + t * (mic_positions(i+1,3) - mic_positions(i,3));

            % 赋值到g.x_gt
            g.x_gt(current_index + 1) = x_coord;
            g.x_gt(current_index + 2) = y_coord;
            g.x_gt(current_index + 3) = z_coord;

            % 更新索引，为下一个采样点准备
            current_index = current_index + 3;
        end
    end


%     % 参数设定
%     r = 2;                    % 螺旋线半径
%     delta_z = 0.1;           % z轴上的步长
%     delta_theta = 2 * pi / 10; % 角速度
% 
%     time_interval = 0.125;    % 每个采样点的时间间隔，单位秒
% 
%     num_points = 120; % 采样点个数
%     threshold1 = 40;
%     threshold2 = 80;
%     threshold3 = 100;
% 
% %     total_time = num_points * time_interval;
%     % 初始化声源位置存储变量
%     g.x_gt(41:41+num_points*3-1) = 0; % 乘以3是因为每个时间点有三个坐标值
% 
%     % 计算每个采样点在螺旋线上的位置
%     for i = 1:threshold1
%         t = (i - 1) * time_interval; % 当前采样点的时间
%         index = 40 + (i - 1) * 3;
%         g.x_gt(index + 1) = r * cos(delta_theta * i) + 0.5; % x坐标
%         g.x_gt(index + 2) = r * sin(delta_theta * i) + 0.5; % y坐标
%         g.x_gt(index + 3) = delta_z * i; % z坐标，从0开始每个间隔上升delta_z
%     end
% 
%     % z轴直线运动
%     for i = (threshold1 + 1):threshold2
%         % 声源的x和y坐标保持最后一个螺旋线坐标不变
%         x_pos = r * cos(delta_theta * threshold1) + 0.5;
%         y_pos = r * sin(delta_theta *threshold1) + 0.5;
%         
%         % z坐标逐渐减小
%         z_pos = delta_z * threshold1 - delta_z * (i - threshold1) * 2;
%         
%         index = 40 + (i - 1) * 3;
%         g.x_gt(index + 1) = x_pos; % x坐标
%         g.x_gt(index + 2) = y_pos; % y坐标
%         g.x_gt(index + 3) = z_pos; % z坐标
%     end
%        % 沿x移动
%     for i = (threshold2 + 1):threshold3
%         x_pos = r * cos(delta_theta * threshold1) + 0.5 + (i -threshold2) * 0.1;
%         y_pos = r * sin(delta_theta *threshold1) + 0.5;
%        
%         z_pos = delta_z * threshold1 - delta_z * (threshold2 - threshold1) * 2;
%         
%         index = 40 + (i - 1) * 3;
%         g.x_gt(index + 1) = x_pos; % x坐标
%         g.x_gt(index + 2) = y_pos; % y坐标
%         g.x_gt(index + 3) = z_pos; % z坐标
%     end
% 
%     % 沿y移动
%     for i = (threshold3 + 1):num_points
%         x_pos = r * cos(delta_theta * threshold1) + 0.5 + (threshold3 -threshold2) * 0.1;
%         y_pos = r * sin(delta_theta *threshold1) + 0.5 + (i - threshold3)*0.1;
%         
%         z_pos = delta_z * threshold1 - delta_z * (threshold2 - threshold1) * 3;
%         
%         index = 40 + (i - 1) * 3;
%         g.x_gt(index + 1) = x_pos; % x坐标
%         g.x_gt(index + 2) = y_pos; % y坐标
%         g.x_gt(index + 3) = z_pos; % z坐标
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %导入真实数据
% %     num_points = 75;
%     addpath('npy-matlab-master/npy-matlab-master/npy-matlab')
%     % 读取.npy文件，使用相对路径
% %     sound_positions = readNPY('data/pose_02.npy');
%     sound_positions = readNPY('data/final_data/pose_2.npy');
%     num_points = size(sound_positions, 1);
%     true_sound_positions = readNPY('data/final_data/pose_true2.npy');
% 
% %     tra = readNPY('data/final_data/traj_1.npy');
% %     % 显示数据
% %     disp(sound_positions);
%     ref_x = -0.20;
%     ref_y = -0.03;
%     ref_z = 0.05;
%     for i = 1: num_points
%         index = 5 * g.M + (i-1) * 3;
%         g.x_gt(index + 1) = true_sound_positions(i,1) + ref_x; % x坐标
%         g.x_gt(index + 2) = true_sound_positions(i,2) + ref_y; % y坐标
%         g.x_gt(index + 3) = true_sound_positions(i,3) + ref_z; % z坐标
%     end
% 
%     g.x = zeros(size(g.x_gt));
%     for i = 1: num_points
%         index = 5 * g.M + (i-1) * 3;
%         g.x(index + 1) = sound_positions(i,1) + ref_x; % x坐标
%         g.x(index + 2) = sound_positions(i,2) + ref_y; % y坐标
%         g.x(index + 3) = sound_positions(i,3) + ref_z; % z坐标
%     end
    
%% 下面生成TDOA测量量
    c = 346;              % 声速 (m/s)
    TDOA_std_dev_noise = 0.025e-3  ; % 噪声的标准偏差,对应文章的R
    % TDOA_std_dev_noise = 1e-3 * c * 0.01;
    % TDOA_std_dev_noise = 0.007; %单位m
    % 声源位置存储在g.x_gt的41到280行中，每三个行对应一个(x, y, z)坐标
    % 计算TDOA测量值并添加到g.edge.measurement中

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

            TDOAs(i - 1) = (dist_mic - dist_first_mic) / c ;
            % TDOAs(i - 1) = (dist_mic - dist_first_mic);
            
            % TDOAs(i - 1) = delay_mean(k, i-1);
        end
        % 加入测量噪声
        TDOAs_noise = TDOA_std_dev_noise * randn(g.M-1, 1);
        TDOAs = TDOAs + TDOAs_noise * 1;
        
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
%% 给g.x一个初值

    % 初始化gx
    init_std_dev = 0.2;
    
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
