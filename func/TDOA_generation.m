function measurement = TDOA_generation()
clc
k = 1;
s = 0.5;
t = 1:43;
measurement = zeros(length(t)*6,7);
%%
for p = 1:length(t)
    
    % 使用 sprintf 生成两位数格式的文件编号
    f = sprintf('%02d', t(p));

    [SIG, fs] = audioread("./9.3audio/"+f+".wav");

    % total_samples = size(SIG, 1);  % 获取音频信号的总样本数
    % time_vector = (0:total_samples-1) / fs;  % 转换为时间轴，以秒为单位
    % 
    % figure;
    % 
    % % 循环绘制每个声道的波形在不同的子图中
    % for j = 1:8
    %     subplot(8, 1, j);  % 创建8行1列的子图，当前位于第j个子图
    %     plot(time_vector, SIG(:, j));  % 绘制第j个声道的信号
    %     xlabel('Time (s)');
    %     ylabel('Amplitude');
    %     title(['Mic ', num2str(j)]);
    %     grid on;
    % end
    
    % for i = 1:6
    %     TDOA=[];
    %     a = fs*(s+1*(i-1))+1;
    %     b = fs*(s+1*i);
    %     for j = 1:7
    %         refsig = SIG(a:b,j);
    %         sig = SIG(a:b,j+1:8);
    %         [tau, cc] = gccphat(sig, refsig, fs);
    %         TDOA = [TDOA tau];
    %     end
    %     measurement(k,:) = TDOA;
    %     k = k+1;
    % end
    for i = 1:6
        TDOA = [];
        a = fs*(s+1*(i-1))+1;
        b = fs*(s+1*i);

        refsig = SIG(a:b, 1);  % 以第1个麦克风作为参考信号

        for j = 2:8  % 从第2个麦克风到第8个麦克风
            sig = SIG(a:b, j);
            [tau, cc] = gccphat(sig, refsig, fs);
            TDOA = [TDOA tau];  % 将结果存储在TDOA向量中
        end

        measurement(k, :) = TDOA;  % 将每次计算的结果存储在measurement矩阵的相应行
        k = k + 1;
    end

end
save('9.3measurement.mat','measurement');
end
