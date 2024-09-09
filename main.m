%% refresh
clear;
close all;
clc;

% rng(0);
%% add path for including some tool functions
addpath('func');
sum = 0;
closed_form_sum = 0;
for i = 1:100
% g = modify_and_save_g_sim(".\data\test_sim.mat");
% fig4a.graph_file = './data/test_sim.mat';

u = source_position_generation();
% measurement = TDOA_generation();

measurement = load("9.3measurement.mat");
outlier = finderror(measurement.measurement);
g = modify_and_save_g_real(".\data\9.3test_real.mat", u, measurement.measurement, outlier);
% 
fig4a.graph_file = './data/9.3test_real.mat';
% fig4a.graph_file = './data/test_sim.mat';

%% params

% Fig.4(a)
fig4a.eps = 1e-2;
fig4a.fig.title = 'Fig.4(a)';
fig4a.fig.view_a = 30; fig4a.fig.view_e = 15;
disp('随机初始值');
error = calib_func(fig4a);
sum = sum + error;
disp('闭式解赋初值');
closed_form_error = closed_form_calib_func(fig4a);
closed_form_sum = closed_form_sum + closed_form_error;
end
disp(sum)
disp(closed_form_sum)

%%%%%%%%%%%%%%%%%%%%仿真实验结果记录
%0.2误差时，总error分别为26.4331,7.7711
%0.1误差时，总error分别为14.5188，7.7084
%0.05误差时，总error分别为7.1249，7.9052
%0.5误差时，总error分别为73.3027,7.8191

%%%%%%%%%%%%%%%%%%%%真实实验结果记录
%0.2误差时，总error分别为34.8909 26.9814
%0.1误差时，总error分别为28.1551 26.9814
%0.05误差时，总error分别为27.0480，26.9814
%0.5误差时，总error分别为72.9131 26.9814