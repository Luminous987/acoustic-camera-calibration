function mic_position = linearize_and_solve_v2(g)

N = length(g.edges); 
M = g.M;

H = zeros((N - 1) * (M - 1), 1);
G = zeros((N - 1) * (M - 1), 3 * (M - 1) + N);
kesai = zeros(3 * (M - 1) + N, 1);
% W = zeros(N * (M - 1), N * (M - 1));
W = 10 * eye((N - 1) * (M - 1), (N - 1) * (M - 1));
epsilon = zeros((N - 1) * (M -1), 1);

%填充G
for j = 1:M-1
    uN = g.x_gt(3*M+3*(N-1)+1 : 3*M+3*(N-1)+3,1);
    %填充左边分块矩阵
    for i = 1:N-1
        ui = g.x_gt(3*M + 3*(i-1) + 1:3*M + 3*(i-1) + 3,1);
        G((j-1)*(N-1)+i, 3*(j-1)+1 : 3*(j-1)+3) = -2 * (ui - uN);
    end
    %填充右边分块矩阵
    for i = 1:N-1
        edge = g.edges(i);
        di = edge.measurement *340;
        dij_1 = di(j);
        G((j-1)*(N-1)+i , 3*(M-1)+i) = -2 * dij_1;
    end
    %填充最右边一列
    edge_N = g.edges(N);
    dN = edge_N.measurement * 340;
    dNj_1 = dN(j);
    G((j-1)*(N-1)+1:(j-1)*(N-1)+(N-1) , 3*(M-1)+N) = 2*dNj_1;
end

%填充H
for j = 1:M-1
    edge_N = g.edges(N);
    dN = edge_N.measurement * 340;
    dNj_1 = dN(j);
    for i = 1:N-1
        edge = g.edges(i);
        di = edge.measurement * 340;
        dij_1 = di(j);
        H((j-1)*(N-1)+i,1) = dij_1^2 - dNj_1^2;
    end
end

% %这里直接将真值带入测试
% B_test = zeros((N - 1) * (M - 1), (N - 1) * (M - 1));
% for j = 1:M-1
%     sj = g.x_gt(j*3+1:j*3+3);
%     for i = 1:N-1
%         ui = g.x_gt(3*M+3*(i-1)+1 : 3*M+3*(i-1)+3,1);
%         B_test((j-1)*(N-1)+i,(j-1)*(N-1)+i) = 2*(norm(uN-sj)-norm(ui-sj));
%     end
% end
% lamda = 10^-1;
% W = inv(B_test * W * B_test + eye((N - 1) * (M - 1))*lamda);
kesai = (G' * W * G) \ (G' * W * H);

%从kesai中还原s,我们只需要还原s1
A = zeros(N-1,3);
b = zeros(N-1,1);
X = zeros(3,1);

uN = g.x_gt(3*M+3*(N-1)+1 : 3*M+3*(N-1)+3,1);
diN = kesai(3*(M-1)+N,1);
for i = 1:N-1
    ui = g.x_gt(3*M + 3*(i-1) + 1:3*M + 3*(i-1) + 3,1);
    di1 = kesai(3*(M-1)+i,1);
    A(i,:) = 2 * (ui - uN)';
    b(i,:) = norm(ui)^2 - norm(uN)^2 + diN^2 - di1^2;
end
% initial estimation
X1 = (A' * A) \ A' * b
S = zeros(21, 1);
ERROR = zeros(24, 1);
ERROR(1:3,1) = abs(X1 - g.x_gt(1:3,1));
for i = 1:7    
    start_idx = (i-1)*3 + 1;
    end_idx = i*3;
    S(start_idx:end_idx, 1) = X1 + kesai(start_idx:end_idx, 1);
    ERROR(start_idx+3:end_idx+3,1) = abs(S(start_idx:end_idx, 1)-g.x_gt(start_idx+3:end_idx+3, 1));
end
% disp(S);
disp(sum(ERROR)/24);
%% 第二步 解算B，代入重新计算
S_init = zeros(3*M,1);
S_init(1:3,1) = X1;
for j = 2:M
    S_init((j-1)*3+1:(j-1)*3+3 , 1) = kesai((j-2)*3+1:(j-2)*3+3 , 1) + X;
end

B = zeros((N - 1) * (M - 1), (N - 1) * (M - 1));
uN = g.x_gt(3*M+3*(N-1)+1 : 3*M+3*(N-1)+3,1);
for j = 1:M-1
    sj = S_init((j-1)*3+1:(j-1)*3+3);
    for i = 1:N-1
        ui = g.x_gt(3*M + 3*(i-1) + 1:3*M + 3*(i-1) + 3,1);
        B((j-1)*(N-1)+i,(j-1)*(N-1)+i) = 2*(norm(uN-sj)-norm(ui-sj));
    end
end
lamda = 10^-1;
W_final = inv(B * W * B + eye((N - 1) * (M - 1))*lamda);
kesai_final = (G' * W_final * G) \ (G' * W_final * H);

%% 再提取S1
A = zeros(N-1,3);
b = zeros(N-1,1);
X1_final = zeros(3,1);

uN = g.x_gt(3*M+3*(N-1)+1 : 3*M+3*(N-1)+3,1);
diN = kesai_final(3*(M-1)+N,1);
for i = 1:N-1
    ui = g.x_gt(3*M + 3*(i-1) + 1:3*M + 3*(i-1) + 3,1);
    di1 = kesai_final(3*(M-1)+i,1);
    A(i,:) = 2 * (ui - uN)';
    b(i,:) = norm(ui)^2 - norm(uN)^2 + diN^2 - di1^2;
end
X1_final = (A' * A) \ A' * b

S = zeros(21, 1);
ERROR_final = zeros(24, 1);
ERROR_final(1:3,1) = abs(X1_final - g.x_gt(1:3,1));
for i = 1:7    
    start_idx = (i-1)*3 + 1;
    end_idx = i*3;
    S(start_idx:end_idx, 1) = X1_final + kesai_final(start_idx:end_idx, 1);
    ERROR_final(start_idx+3:end_idx+3,1) = abs(S(start_idx:end_idx, 1)-g.x_gt(start_idx+3:end_idx+3, 1));
end
% disp(S);
disp('初始ERROR');
disp(sum(ERROR_final)/24);
mic_position = [X1_final; S];
end