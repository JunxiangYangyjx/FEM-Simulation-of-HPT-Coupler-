% =========================================================================
% 阶段一 (最终黄金版 v15.2 - 精确点加权): 最终分析与可视化
%
% 描述:
%   此版本为最终分析脚本。它执行以下所有任务：
%   1. 运行初始模型仿真作为基准。
%   2. (v15.2 修改) 使用对谐振点及其邻近三点进行精确加权的方式进行参数提取。
%   3. 对提取的参数进行后处理平滑。
%   4. 可视化基线电阻 Rij_base 与最终修正电阻 Rij_total 的对比。
%   5. 运行最终验证仿真。
%   6. 计算初始模型和最终拟合模型相对于实验数据的 NRMSE 值。
%   7. 绘制包含实验、初始模型、最终拟合模型的三曲线对比图。
%   8. 在 Cs1 和 Cs2 支路中加入了可设置的串联电阻 Rs1 和 Rs2。
% =========================================================================
clear; clc; close all;
%% 1. 加载实验数据和定义固定/初始参数
% -------------------------------------------------------------------------
disp('正在加载实验数据和固定参数...');
try
    data = readmatrix('doublepcboutinshort1.5m3.5m.csv');
    freq_exp = data(:,1);
    Zmag_exp = data(:,2);
    Zang_exp = data(:,3);
    Z_complex_exp = Zmag_exp .* exp(1j * Zang_exp * pi/180);
    disp('实验数据加载成功。');
catch ME
    error('无法加载CSV文件! 错误: %s', ME.message);
end

% --- 固定参数定义 ---
fixed_params.R_L=1e-6; 
fixed_params.Cs1=28.4e-12; 
fixed_params.Cs2=28.4e-12;
fixed_params.L11=39.99e-6; 
fixed_params.L12=fixed_params.L11;
fixed_params.L21=34.07e-6; 
fixed_params.L22=fixed_params.L21;
fixed_params.Rs1 = 0; % Ohm
fixed_params.Rs2 = fixed_params.Rs1;
k1=0.2422; k2=0.2422; k3=0.65; k4=0.62; k5=0.234; k6=0.234;
fixed_params.M1=-k1*sqrt(fixed_params.L11*fixed_params.L21);
fixed_params.M2=-k2*sqrt(fixed_params.L12*fixed_params.L22);
fixed_params.M3=k3*sqrt(fixed_params.L11*fixed_params.L12);
fixed_params.M4=k4*sqrt(fixed_params.L21*fixed_params.L22);
fixed_params.M5=-k5*sqrt(fixed_params.L11*fixed_params.L22);
fixed_params.M6=-k6*sqrt(fixed_params.L21*fixed_params.L12);
fixed_params.V_S=20;
fixed_params.R_p_stabilizer=1e8;
fixed_params.RM_fixed = 0.0;
initial_params_scaled = [0, 0, 5.0]; % [ΔR11, ΔR21, CM_pF]
disp('固定参数定义完毕。');


%% 2. 运行初始模型仿真 (ΔR=0)
% -------------------------------------------------------------------------
disp('正在运行使用初始参数的基线模型仿真 (ΔR=0)...');
Z_in_initial = zeros(length(freq_exp), 1, 'like', 1j);
for i = 1:length(freq_exp)
    Z_in_initial(i) = run_single_freq_sim_3D(freq_exp(i), initial_params_scaled, fixed_params);
end
disp('基线模型仿真完成。');


%% 3. 定位加权点
% -------------------------------------------------------------------------
disp('正在定位串联谐振点及其邻近点以进行加权...');
[~, min_idx] = min(Zmag_exp);
f_series = freq_exp(min_idx);
fprintf('已定位串联谐振点在: %.4f MHz (数据点索引: %d)\n', f_series / 1e6, min_idx);
weight_amplitude = 1000; % 定义权重幅值
disp('权重将应用于此索引以及前后各一个点。');


%% 4. 逐点优化
% -------------------------------------------------------------------------
num_freq_points = length(freq_exp);
extracted_params_raw = zeros(num_freq_points, 3);
lower_bounds = [-20, -20, 0.01]; 
upper_bounds = [20, 20, 25];
current_guess = initial_params_scaled;
lambda_dR11 = 5e-2; lambda_dR21 = 5e-2; lambda_C = 1e-6;
lambda_smooth = 0.05;
options = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'sqp');
disp('开始进行三维参数的逐点反演...');
tic;
p_previous = current_guess;
for i = 1:num_freq_points
    f = freq_exp(i);
    Z_target = Z_complex_exp(i);
    
    % 【修改】计算频率权重：对谐振点及其邻近的±1个点进行加权
    if abs(i - min_idx) <= 1
        weight = 1 + weight_amplitude; % 应用高权重
    else
        weight = 1; % 其他所有点使用基础权重
    end
    
    % 定义目标函数 (包含加权、正则化和平滑项)
    objective_fun = @(p) weight * abs(run_single_freq_sim_3D(f, p, fixed_params) - Z_target)^2 ...
        + lambda_dR11 * p(1)^2 + lambda_dR21 * p(2)^2 + lambda_C * p(3)^2 ...
        + lambda_smooth * sum((p - p_previous).^2);
        
    [p_optimal, ~] = fmincon(objective_fun, current_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
    
    extracted_params_raw(i, :) = p_optimal;
    current_guess = p_optimal;
    p_previous = p_optimal;
end
elapsed_time = toc;
fprintf('参数提取完成! 总耗时: %.2f 秒。\n', elapsed_time);


%% 5. 后处理平滑
% (此部分及之后的所有代码均无需修改)
% -------------------------------------------------------------------------
disp('正在对提取的参数进行后处理平滑...');
window_length = 51; 
poly_order = 2;
extracted_params_smooth = zeros(size(extracted_params_raw));
extracted_params_smooth(:,1) = smoothdata(extracted_params_raw(:,1), 'sgolay', window_length, 'Degree', poly_order);
extracted_params_smooth(:,2) = smoothdata(extracted_params_raw(:,2), 'sgolay', window_length, 'Degree', poly_order);
extracted_params_smooth(:,3) = smoothdata(extracted_params_raw(:,3), 'sgolay', window_length, 'Degree', poly_order);
disp('正在保存平滑后的参数数据以供阶段二使用...');
save('extracted_params_3D_smooth.mat', 'freq_exp', 'extracted_params_smooth');
disp('已保存为 extracted_params_3D_smooth.mat');

%% 导出提取的 C_M 曲线到 Excel
disp('正在导出 C_M(f) 曲线到 Excel 文件...');

% 构建表格
T_CM = table( ...
    freq_exp/1e6, ...                                % MHz
    extracted_params_raw(:,3), ...                   % 原始 C_M
    extracted_params_smooth(:,3), ...                % 平滑 C_M
    'VariableNames', {'Frequency_MHz', 'C_M_raw_pF', 'C_M_smooth_pF'} ...
);

% 导出到 Excel 文件
output_filename = 'Extracted_CM_Curve.xlsx';
writetable(T_CM, output_filename);

disp(['C_M 曲线已成功导出到 ', output_filename]);

%% 6. (图 1) 可视化原始提取参数 vs 平滑后参数 (*** 最终专业格式化版 ***)
% -------------------------------------------------------------------------
% --- 选择要绘制的离散点 ---
num_points_to_plot = 50; % 定义要绘制的离散点数量
indices_to_plot = round(linspace(1, length(freq_exp), num_points_to_plot));
% --- 修改结束 ---

figure('Name', 'Raw vs. Smoothed Extracted Parameters', 'Position', [100, 100, 900, 700]);

% --- Subplot 1: ΔR11 ---
subplot(3,1,1);
% 原始数据: 蓝色透明圆圈, 线宽1.5
plot(freq_exp(indices_to_plot)/1e6, extracted_params_raw(indices_to_plot,1), 'o', 'MarkerEdgeColor', 'b', 'LineWidth', 1.5, 'DisplayName', 'Raw Extracted $\Delta R_{11}$');
hold on;
% 平滑数据: 红色实心方块
plot(freq_exp(indices_to_plot)/1e6, extracted_params_smooth(indices_to_plot,1), 's', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'DisplayName', 'Smoothed $\Delta R_{11}$');
xlabel('Frequency (MHz)', 'Interpreter', 'latex');
ylabel('$\Delta R_{11}$ (Ohms)', 'Interpreter', 'latex'); % <-- 修改了此处的Y轴标签
grid on;
legend('Location','northwest', 'Interpreter', 'latex', 'FontSize', 10);
set(gca, 'FontSize', 12, 'XTickLabel', []); % 设置字体大小, 并移除顶部图的X轴刻度
hold off;

% --- Subplot 2: ΔR21 ---
subplot(3,1,2);
% 原始数据: 蓝色透明圆圈, 线宽1.5
plot(freq_exp(indices_to_plot)/1e6, extracted_params_raw(indices_to_plot,2), 'o', 'MarkerEdgeColor', 'b', 'LineWidth', 1.5, 'DisplayName', 'Raw Extracted $\Delta R_{21}$');
hold on;
% 平滑数据: 红色实心方块
plot(freq_exp(indices_to_plot)/1e6, extracted_params_smooth(indices_to_plot,2), 's', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'DisplayName', 'Smoothed $\Delta R_{21}$');
xlabel('Frequency (MHz)', 'Interpreter', 'latex');
ylabel('$\Delta R_{21}$ (Ohms)', 'Interpreter', 'latex'); % <-- Y轴标签保持为通用的 ΔR
grid on;
legend('Location','northwest', 'Interpreter', 'latex', 'FontSize', 10);
set(gca, 'FontSize', 12, 'XTickLabel', []); % 设置字体大小, 并移除中间图的X轴刻度
hold off;

% --- Subplot 3: CM ---
subplot(3,1,3);
% 原始数据: 蓝色透明圆圈, 线宽1.5
plot(freq_exp(indices_to_plot)/1e6, extracted_params_raw(indices_to_plot,3), 'o', 'MarkerEdgeColor', 'b', 'LineWidth', 1.5, 'DisplayName', 'Raw Extracted $C_M$');
hold on;
% 平滑数据: 红色实心方块
plot(freq_exp(indices_to_plot)/1e6, extracted_params_smooth(indices_to_plot,3), 's', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'DisplayName', 'Smoothed $C_M$');
xlabel('Frequency (MHz)', 'Interpreter', 'latex');
ylabel('$C_\mathrm{M}$ (pF)', 'Interpreter', 'latex');
grid on;
legend('Location','northwest', 'Interpreter', 'latex', 'FontSize', 10);
set(gca, 'FontSize', 12); % 设置字体大小
hold off;

%% 7. (图 2) 可视化基线电阻 vs 修正后电阻
% -------------------------------------------------------------------------
freq_MHz = freq_exp / 1e6;
R11_base = 8.1969e-9*freq_exp.^1.318 + 2.00;
R21_base = 5.9395e-8*freq_exp.^1.181 + 1.66;
R11_total = R11_base + extracted_params_smooth(:,1);
R21_total = R21_base + extracted_params_smooth(:,2);
figure('Name', 'Resistance Model Comparison', 'Position', [150, 150, 900, 600]);
subplot(2,1,1);
plot(freq_MHz, R11_base, 'r:', 'LineWidth', 2, 'DisplayName', 'R_{11,base}(f) (Single Coil Fit)');
hold on;
plot(freq_MHz, R11_total, 'b-', 'LineWidth', 2, 'DisplayName', 'R_{11,total}(f) = R_{base} + \DeltaR_{fit}');
hold off;
title('Comparison of Base and Final Corrected R_{11}');
xlabel('Frequency (MHz)'); ylabel('Resistance (Ohms)'); grid on; legend('Location', 'best');
subplot(2,1,2);
plot(freq_MHz, R21_base, 'r:', 'LineWidth', 2, 'DisplayName', 'R_{21,base}(f) (Single Coil Fit)');
hold on;
plot(freq_MHz, R21_total, 'b-', 'LineWidth', 2, 'DisplayName', 'R_{21,total}(f) = R_{base} + \DeltaR_{fit}');
hold off;
title('Comparison of Base and Final Corrected R_{21}');
xlabel('Frequency (MHz)'); ylabel('Resistance (Ohms)'); grid on; legend('Location', 'best');

%% 8. 使用最终平滑的参数进行验证仿真
% -------------------------------------------------------------------------
disp('正在使用最终平滑的参数进行验证仿真...');
Z_in_fitted = zeros(num_freq_points, 1, 'like', 1j);
for i = 1:num_freq_points
    Z_in_fitted(i) = run_single_freq_sim_3D(freq_exp(i), extracted_params_smooth(i,:), fixed_params);
end

%% 9. 计算 NRMSE 值
% -------------------------------------------------------------------------
disp('正在计算 NRMSE 误差...');
rmse_mag_initial = sqrt(mean((abs(Z_in_initial) - Zmag_exp).^2));
rmse_mag_fitted  = sqrt(mean((abs(Z_in_fitted) - Zmag_exp).^2));
norm_factor_mag  = max(Zmag_exp) - min(Zmag_exp);
nrmse_mag_initial = rmse_mag_initial / norm_factor_mag;
nrmse_mag_fitted  = rmse_mag_fitted / norm_factor_mag;
rmse_ang_initial = sqrt(mean((angle(Z_in_initial)*180/pi - Zang_exp).^2));
rmse_ang_fitted  = sqrt(mean((angle(Z_in_fitted)*180/pi - Zang_exp).^2));
norm_factor_ang  = max(Zang_exp) - min(Zang_exp);
nrmse_ang_initial = rmse_ang_initial / norm_factor_ang;
nrmse_ang_fitted  = rmse_ang_fitted / norm_factor_ang;
fprintf('\n--- 性能评估 (NRMSE, 越小越好) ---\n');
fprintf('NRMSE (Magnitude) - Initial Model (ΔR=0): %.2f%%\n', nrmse_mag_initial * 100);
fprintf('NRMSE (Magnitude) - Final Fitted Model:     %.2f%%\n', nrmse_mag_fitted * 100);
fprintf('----------------------------------------\n');
fprintf('NRMSE (Phase)     - Initial Model (ΔR=0): %.2f%%\n', nrmse_ang_initial * 100);
fprintf('NRMSE (Phase)     - Final Fitted Model:     %.2f%%\n', nrmse_ang_fitted * 100);
fprintf('----------------------------------------\n\n');

%% 10. (图 3 & 4) 最终验证绘图 (三条曲线对比)
% -------------------------------------------------------------------------
figure('Name', 'Final Validation: Magnitude Comparison', 'Position', [200, 200, 900, 500]);
semilogy(freq_MHz, Zmag_exp, 'k--', 'LineWidth', 2.5, 'DisplayName', 'Experimental Data');
hold on;
% semilogy(freq_MHz, abs(Z_in_initial), 'r:', 'LineWidth', 2, 'DisplayName', 'Initial Model');
semilogy(freq_MHz, abs(Z_in_fitted), 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Model');
hold off;
title('Magnitude Response Comparison');
xlabel('Frequency (MHz)'); ylabel('|Z_{in}| (Ohms)');
grid on; legend('Location', 'best'); set(gca, 'YScale', 'log'); ylim([10, 1e6]);
figure('Name', 'Final Validation: Phase Comparison', 'Position', [250, 150, 900, 500]);
plot(freq_MHz, Zang_exp, 'k--', 'LineWidth', 2.5, 'DisplayName', 'Experimental Data');
hold on;
% plot(freq_MHz, angle(Z_in_initial)*180/pi, 'r:', 'LineWidth', 2, 'DisplayName', 'Initial Model');
plot(freq_MHz, angle(Z_in_fitted)*180/pi, 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Model');
hold off;
title('Phase Response Comparison');
xlabel('Frequency (MHz)'); ylabel('Phase (Degrees)');
grid on; legend('Location', 'best'); ylim([-100, 100]);
disp('脚本运行完毕。');

% =========================================================================
%                        辅助函数 (与之前 v15.1 完全相同)
% =========================================================================
function Z_out = run_single_freq_sim_3D(f, params_scaled, fixed_params)
delta_R11=params_scaled(1); delta_R21=params_scaled(2);
CM=params_scaled(3)*1e-12; omega=2*pi*f;
R_L=fixed_params.R_L; RM=fixed_params.RM_fixed;
Cs1=fixed_params.Cs1; Cs2=fixed_params.Cs2;
L11=fixed_params.L11; L12=fixed_params.L12; L21=fixed_params.L21;
L22=fixed_params.L22; M1=fixed_params.M1; M2=fixed_params.M2;
M3=fixed_params.M3; M4=fixed_params.M4; M5=fixed_params.M5;
M6=fixed_params.M6; V_S=fixed_params.V_S; R_p=fixed_params.R_p_stabilizer;
Rs1 = fixed_params.Rs1; Rs2 = fixed_params.Rs2;
R11_base=8.1969e-9*f^1.318+2.00;
R21_base=5.9395e-8*f^1.181+1.66;
R11_f=R11_base+delta_R11; R12_f=R11_f;
R21_f=R21_base+delta_R21; R22_f=R21_f;
Y_full=zeros(8,8,'like',1j); Y_RL=1/R_L;
Z_series1=Rs1+1/(1j*omega*Cs1); Y_Cs1=1/Z_series1+1/R_p;
Z_series2=Rs2+1/(1j*omega*Cs2); Y_Cs2=1/Z_series2+1/R_p;
Y_CM1=1/(RM+1/(1j*omega*CM)); Y_CM2=Y_CM1;
n1=2;n2=3;n3=4;n4=5;n5=6;n6=7;
Y_full(n3,n3)=Y_full(n3,n3)+Y_RL;Y_full(n4,n4)=Y_full(n4,n4)+Y_RL;Y_full(n3,n4)=Y_full(n3,n4)-Y_RL;Y_full(n4,n3)=Y_full(n4,n3)-Y_RL;
Y_full(n1,n1)=Y_full(n1,n1)+Y_Cs1;Y_full(n6,n6)=Y_full(n6,n6)+Y_Cs1;Y_full(n1,n6)=Y_full(n1,n6)-Y_Cs1;Y_full(n6,n1)=Y_full(n6,n1)-Y_Cs1;
Y_full(n2,n2)=Y_full(n2,n2)+Y_Cs2;Y_full(n5,n5)=Y_full(n5,n5)+Y_Cs2;Y_full(n2,n5)=Y_full(n2,n5)-Y_Cs2;Y_full(n5,n2)=Y_full(n5,n2)-Y_Cs2;
Y_full(n1,n1)=Y_full(n1,n1)+Y_CM1;Y_full(n2,n2)=Y_full(n2,n2)+Y_CM1;Y_full(n1,n2)=Y_full(n1,n2)-Y_CM1;Y_full(n2,n1)=Y_full(n2,n1)-Y_CM1;
Y_full(n6,n6)=Y_full(n6,n6)+Y_CM2;Y_full(n5,n5)=Y_full(n5,n5)+Y_CM2;Y_full(n6,n5)=Y_full(n6,n5)-Y_CM2;Y_full(n5,n6)=Y_full(n5,n6)-Y_CM2;
R_diag=diag([R11_f,R21_f,R12_f,R22_f]);
L_matrix=[L11,M1,M3,M5;M1,L21,M6,M4;M3,M6,L12,M2;M5,M4,M2,L22];
Z_branch=R_diag+1j*omega*L_matrix;
Y_branch=Z_branch\eye(size(Z_branch));
C=zeros(8,4); C(1,1)=1;C(2,1)=-1;C(8,2)=1;C(7,2)=-1;
C(3,3)=1;C(4,3)=-1;C(6,4)=1;C(5,4)=-1;
Y_coupled = C * Y_branch * C';
Y_full=Y_full+Y_coupled;
k_idx=[2,3,4,5,6,7];s_idx=[1,8];
Y_kk=Y_full(k_idx,k_idx);Y_ks=Y_full(k_idx,s_idx);
Y_sk=Y_full(s_idx,k_idx);Y_ss=Y_full(s_idx,s_idx);
V_s_vec=[V_S;0];
b = -Y_ks * V_s_vec;
V_k=Y_kk\b;
I_s = Y_sk * V_k + Y_ss * V_s_vec;
i_in=I_s(1);
Z_out=V_S/i_in;
end

%% 11. 计算 Figure 1 和 Figure 3 的 RMSE
% =========================================================================
disp(' '); % 添加一个空行以美化输出
disp('--- 正在计算 Figure 1 和 Figure 3 的 RMSE 值 ---');

% --- Figure 1: 原始提取参数 vs 平滑后参数的 RMSE ---
% 这个指标衡量了Savitzky-Golay平滑滤波器对原始反演数据的"平滑"或"修改"程度。
% 值越小，说明平滑曲线与原始数据点越贴近。

rmse_dR11_smoothing = sqrt(mean((extracted_params_smooth(:,1) - extracted_params_raw(:,1)).^2));
rmse_dR21_smoothing = sqrt(mean((extracted_params_smooth(:,2) - extracted_params_raw(:,2)).^2));
rmse_CM_smoothing   = sqrt(mean((extracted_params_smooth(:,3) - extracted_params_raw(:,3)).^2));

fprintf('\n--- Figure 1: 平滑操作引入的RMSE ---\n');
fprintf('RMSE (ΔR11, Raw vs. Smoothed): %.4f Ohms\n', rmse_dR11_smoothing);
fprintf('RMSE (ΔR21, Raw vs. Smoothed): %.4f Ohms\n', rmse_dR21_smoothing);
fprintf('RMSE (CM,   Raw vs. Smoothed): %.4f pF\n', rmse_CM_smoothing);
fprintf('----------------------------------------\n');


%% 12. 计算 Figure 1 的 NRMSE (平滑效果量化)
% =========================================================================
disp('--- 正在计算 Figure 1 的 NRMSE 值 (量化平滑效果) ---');

% --- NRMSE for ΔR11 ---
rmse_dR11_smoothing = sqrt(mean((extracted_params_smooth(:,1) - extracted_params_raw(:,1)).^2));
norm_factor_dR11  = max(extracted_params_raw(:,1)) - min(extracted_params_raw(:,1));
nrmse_dR11_smoothing = rmse_dR11_smoothing / norm_factor_dR11;

% --- NRMSE for ΔR21 ---
rmse_dR21_smoothing = sqrt(mean((extracted_params_smooth(:,2) - extracted_params_raw(:,2)).^2));
norm_factor_dR21  = max(extracted_params_raw(:,2)) - min(extracted_params_raw(:,2));
nrmse_dR21_smoothing = rmse_dR21_smoothing / norm_factor_dR21;

% --- NRMSE for CM ---
rmse_CM_smoothing   = sqrt(mean((extracted_params_smooth(:,3) - extracted_params_raw(:,3)).^2));
norm_factor_CM    = max(extracted_params_raw(:,3)) - min(extracted_params_raw(:,3));
nrmse_CM_smoothing  = rmse_CM_smoothing / norm_factor_CM;

fprintf('\n--- Figure 1: 平滑操作引入的 NRMSE ---\n');
fprintf('说明: 该值表示平滑曲线与原始提取数据点的偏离程度，\n');
fprintf('      是原始数据动态范围的一个百分比。\n');
fprintf('-----------------------------------------------------\n');
fprintf('NRMSE (ΔR11, Raw vs. Smoothed): %.2f%%\n', nrmse_dR11_smoothing * 100);
fprintf('NRMSE (ΔR21, Raw vs. Smoothed): %.2f%%\n', nrmse_dR21_smoothing * 100);
fprintf('NRMSE (CM,   Raw vs. Smoothed): %.2f%%\n', nrmse_CM_smoothing * 100);
fprintf('-----------------------------------------------------\n\n');