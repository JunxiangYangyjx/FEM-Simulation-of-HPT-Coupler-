clear; clc; close all;

%% ========================================================================
%% 1) 用户参数（基波模型 + SS 一次侧磁芯损耗）
%% ========================================================================
V_INV_PEAK = 80;                     % [V] 逆变器方波峰值
V1_rms = 4*V_INV_PEAK/(pi*sqrt(2));  % 基波 RMS 电压

R_DC = 500;                          % [Ohm] DC 负载
R_AC_eq = 8/pi^2 * R_DC;             % 基波等效 AC 负载

k_core = 0;                        % [Ohm / A^2] SS 一次侧磁芯损耗系数

%% ========================================================================
%% 2) 载入参数查找表
%% ========================================================================
disp('正在加载参数查找表...');
load('extracted_params_3D_smooth.mat', ...
     'freq_exp','extracted_params_smooth');
freq_sim   = freq_exp;
params_tbl = extracted_params_smooth;
disp('参数加载完成。');

%% ========================================================================
%% 3) Coupler 参数（无铁氧体）
%% ========================================================================
fixed.R_L     = R_AC_eq;
fixed.L_load  = 0;
fixed.Cs1     = 28.4e-12;
fixed.Cs2     = 28.4e-12;
fixed.Rs1     = 0;
fixed.Rs2     = 0;
fixed.L11     = 39.99e-6;
fixed.L12     = 39.99e-6;
fixed.L21     = 34.07e-6;
fixed.L22     = 34.07e-6;

k1=0.2422; k2=0.2422; k3=0.73;
k4=0.713;  k5=0.238;  k6=0.238;

fixed.M1 = -k1*sqrt(fixed.L11*fixed.L21);
fixed.M2 = -k2*sqrt(fixed.L12*fixed.L22);
fixed.M3 =  k3*sqrt(fixed.L11*fixed.L12);
fixed.M4 =  k4*sqrt(fixed.L21*fixed.L22);
fixed.M5 = -k5*sqrt(fixed.L11*fixed.L22);
fixed.M6 = -k6*sqrt(fixed.L21*fixed.L12);

fixed.V_S      = 0;
fixed.RM_fixed = 0.0;

%% ========================================================================
%% 4) SS 前端参数（铁氧体在一次侧）
%% ========================================================================
ss.Lp = 7.54e-6;
ss.Cp = 800e-12;
ss.Ls = 30.4e-6;
ss.Cs = 172e-12;
ss.k  = 0.54;
ss.Rp = 0.05;        % 一次侧铜损
ss.Rs = 0.05;        % 二次侧铜损
ss.M  = ss.k * sqrt(ss.Lp * ss.Ls);

%% ========================================================================
%% 5) 扫频（仅基波）
%% ========================================================================
Nf = numel(freq_sim);

Pin_AC   = zeros(Nf,1);
Pout_AC  = zeros(Nf,1);
eta_AC   = zeros(Nf,1);
Ip_vec   = zeros(Nf,1);
Rcore_p  = zeros(Nf,1);

disp('开始基波扫频（SS 一次侧磁芯损耗）...');
for i = 1:Nf
    f = freq_sim(i);
    w = 2*pi*f;

    % ---- Coupler 参数插值（基波）----
    prm = interp1(freq_sim, params_tbl, f, 'linear', 'extrap');
    prm(:,3) = max(prm(:,3),1e-3);

    [Zin,~,~,~] = run_single_freq_sim_3D_detail_thevenin( ...
        f, prm, fixed, 1, 0);

    % ---- 第一次：理想 SS（用于估计一次侧电流）----
    Zp0 = ss.Rp + 1j*w*ss.Lp + 1/(1j*w*ss.Cp);
    Zs0 = ss.Rs + 1j*w*ss.Ls + 1/(1j*w*ss.Cs);
    Zm  = 1j*w*ss.M;

    denom0 = Zp0*(Zs0 + Zin) - Zm^2;
    if abs(denom0) < eps
        Is0 = 0;
    else
        Is0 = V1_rms * Zm / denom0;
    end
    Ip0 = Is0 * (Zs0 + Zin) / Zm;

    % ---- SS 一次侧磁芯损耗（功率相关）----
    R_core_p = k_core * abs(Ip0)^2;
    Rcore_p(i) = R_core_p;

    % ---- 含磁芯损耗的 SS 阻抗 ----
    Zp = (ss.Rp + R_core_p) + 1j*w*ss.Lp + 1/(1j*w*ss.Cp);
    Zs = ss.Rs + 1j*w*ss.Ls + 1/(1j*w*ss.Cs);

    denom = Zp*(Zs + Zin) - Zm^2;
    if abs(denom) < eps
        Is = 0;
    else
        Is = V1_rms * Zm / denom;
    end
    Ip = Is*(Zs + Zin)/Zm;

    % ---- 输入功率（基波 AC）----
    Pin_AC(i) = real(V1_rms * conj(Ip));

    % ---- 输出功率（基波 AC）----
    Vin = Is * Zin;
    [~,~,V3t,V4t] = run_single_freq_sim_3D_detail_thevenin( ...
        f, prm, fixed, 1, 0);

    V3 = V3t * Vin;
    V4 = V4t * Vin;

    IL = (V3 - V4)/fixed.R_L;
    Pout_AC(i) = abs(IL)^2 * fixed.R_L;

    Ip_vec(i) = abs(Ip);
    eta_AC(i) = Pout_AC(i)/max(Pin_AC(i),eps);
end
disp('扫频完成。');

%% ========================================================================
%% 6) 全频结果绘图
%% ========================================================================
fMHz = freq_sim/1e6;

figure;
subplot(3,1,1)
plot(fMHz, Pin_AC,'b', fMHz, Pout_AC,'r','LineWidth',1.5);
legend('P_{in,AC,fund}','P_{out,AC,fund}');
ylabel('Power (W)');
grid on;

subplot(3,1,2)
plot(fMHz, eta_AC*100,'k','LineWidth',1.5);
ylabel('Efficiency (%)');
grid on;

subplot(3,1,3)
plot(fMHz, Ip_vec,'m','LineWidth',1.5);
xlabel('Frequency (MHz)');
ylabel('I_{SS,pri,fund} (A RMS)');
grid on;

%% ========================================================================
%% 7) 局部频段：2.1–2.5 MHz
%% ========================================================================
idx = (fMHz >= 2.1) & (fMHz <= 2.5);

figure('Color','w');
plot(fMHz(idx), Pout_AC(idx),'r','LineWidth',2);
grid on; grid minor;
xlabel('Frequency (MHz)');
ylabel('P_{out,AC,fund} (W)');
title('Fundamental AC Output Power (2.1–2.5 MHz)');

%% ========================================================================
%% 8) 单点审计
%% ========================================================================
f0 = 2.42e6;
[~,id] = min(abs(freq_sim - f0));

fprintf('\n====== FUNDAMENTAL + SS PRIMARY CORE LOSS ======\n');
fprintf('Frequency               : %.3f MHz\n', freq_sim(id)/1e6);
fprintf('Pout_AC_fund            : %.2f W\n', Pout_AC(id));
fprintf('Pin_AC_fund             : %.2f W\n', Pin_AC(id));
fprintf('Efficiency              : %.2f %%\n', eta_AC(id)*100);
fprintf('I_SS_primary_fund       : %.3f A\n', Ip_vec(id));
fprintf('R_core_primary          : %.3f Ohm\n', Rcore_p(id));
fprintf('===============================================\n');

%% ========================================================================
%% 9) Coupler 求解函数（保持不变）
%% ========================================================================
function [Z_out, I_in, V3, V4] = ...
    run_single_freq_sim_3D_detail_thevenin(f, params_scaled, fixed, V_th, Z_th)

    delta_R11 = params_scaled(1);
    delta_R21 = params_scaled(2);
    CM = params_scaled(3)*1e-12;
    w = 2*pi*f;

    R_L = fixed.R_L;
    RM  = fixed.RM_fixed;

    Cs1 = fixed.Cs1; Cs2 = fixed.Cs2;
    Rs1 = fixed.Rs1; Rs2 = fixed.Rs2;
    L11 = fixed.L11; L12 = fixed.L12;
    L21 = fixed.L21; L22 = fixed.L22;

    M1 = fixed.M1; M2 = fixed.M2;
    M3 = fixed.M3; M4 = fixed.M4;
    M5 = fixed.M5; M6 = fixed.M6;

    R11 = max(8.1969e-9*f^1.318 + 2.00 + delta_R11, 1e-6);
    R21 = max(5.9395e-8*f^1.181 + 1.66 + delta_R21, 1e-6);
    R12 = R11; R22 = R21;

    n1=2; n2=3; n3=4; n4=5; n5=6; n6=7;

    Y = zeros(8,8,'like',1j);
    YL = 1/(R_L + 1j*w*fixed.L_load);

    Y(n3,n3)=Y(n3,n3)+YL; Y(n4,n4)=Y(n4,n4)+YL;
    Y(n3,n4)=Y(n3,n4)-YL; Y(n4,n3)=Y(n4,n3)-YL;

    YC1 = 1/(Rs1 + 1/(1j*w*Cs1));
    YC2 = 1/(Rs2 + 1/(1j*w*Cs2));

    Y(n1,n1)=Y(n1,n1)+YC1; Y(n6,n6)=Y(n6,n6)+YC1;
    Y(n1,n6)=Y(n1,n6)-YC1; Y(n6,n1)=Y(n6,n1)-YC1;

    Y(n2,n2)=Y(n2,n2)+YC2; Y(n5,n5)=Y(n5,n5)+YC2;
    Y(n2,n5)=Y(n2,n5)-YC2; Y(n5,n2)=Y(n5,n2)-YC2;

    YCM = 1/(RM + 1/(1j*w*CM));
    Y(n1,n1)=Y(n1,n1)+YCM; Y(n2,n2)=Y(n2,n2)+YCM;
    Y(n1,n2)=Y(n1,n2)-YCM; Y(n2,n1)=Y(n2,n1)-YCM;

    Y(n6,n6)=Y(n6,n6)+YCM; Y(n5,n5)=Y(n5,n5)+YCM;
    Y(n6,n5)=Y(n6,n5)-YCM; Y(n5,n6)=Y(n5,n6)-YCM;

    Rdiag = diag([R11 R21 R12 R22]);
    Lmat = [L11 M1 M3 M5;
            M1 L21 M6 M4;
            M3 M6 L12 M2;
            M5 M4 M2 L22];

    Ybranch = (Rdiag + 1j*w*Lmat)\eye(4);
    C = zeros(8,4);
    C(1,1)=1; C(2,1)=-1;
    C(8,2)=1; C(7,2)=-1;
    C(3,3)=1; C(4,3)=-1;
    C(6,4)=1; C(5,4)=-1;

    Y = Y + C*Ybranch*C.';

    k = [2 3 4 5 6 7]; s=[1 8];
    Vk = -Y(k,k)\(Y(k,s)*[V_th;0]);

    Vs = [V_th;0];
    Is = Y(s,k)*Vk + Y(s,s)*Vs;

    I_in = Is(1);
    Z_out = V_th/I_in;

    V = zeros(8,1,'like',1j);
    V(k)=Vk; V(s)=Vs;
    V3 = V(n3); V4 = V(n4);
end

%% ========================================================================
%% 10) 绘制系统输入阻抗 (双纵坐标：模值 + 相位)
%% ========================================================================

% 筛选 1.5MHz - 3.0MHz 的频率范围
idx_range = (freq_sim >= 1.5e6) & (freq_sim <= 3.0e6);
f_plot = freq_sim(idx_range);
Nf_plot = numel(f_plot);

Z_in_sys = zeros(Nf_plot, 1);

% 计算阻抗数据
for i = 1:Nf_plot
    f = f_plot(i);
    w = 2*pi*f;
    prm = interp1(freq_sim, params_tbl, f, 'linear', 'extrap');
    prm(:,3) = max(prm(:,3),1e-3);
    [Zin_coupler,~,~,~] = run_single_freq_sim_3D_detail_thevenin(f, prm, fixed, 1, 0);

    % 估计电流以包含 R_core 损耗
    Zm = 1j*w*ss.M;
    Zp0 = ss.Rp + 1j*w*ss.Lp + 1/(1j*w*ss.Cp);
    Zs0 = ss.Rs + 1j*w*ss.Ls + 1/(1j*w*ss.Cs);
    Is0 = V1_rms * Zm / (Zp0*(Zs0 + Zin_coupler) - Zm^2);
    Ip0 = Is0 * (Zs0 + Zin_coupler) / Zm;
    R_core_p_local = k_core * abs(Ip0)^2;

    % 最终阻抗计算
    Zp = (ss.Rp + R_core_p_local) + 1j*w*ss.Lp + 1/(1j*w*ss.Cp);
    Zs = ss.Rs + 1j*w*ss.Ls + 1/(1j*w*ss.Cs);
    Z_in_sys(i) = Zp - (Zm^2) / (Zs + Zin_coupler);
end

% --- 开始绘图 ---
figure('Color', 'w', 'Position', [200, 200, 800, 500]);
hold on;

% 左侧纵坐标：阻抗模值
yyaxis left
plot(f_plot/1e6, abs(Z_in_sys), 'b-', 'LineWidth', 2); % 蓝色实线
ylabel('Magnitude |Z_{in}| (\Omega)', 'Color', 'b');
set(gca, 'YColor', 'b'); % 设置坐标轴颜色与曲线一致
grid on; grid minor;

% 右侧纵坐标：相位
yyaxis right
plot(f_plot/1e6, angle(Z_in_sys)*180/pi, 'r--', 'LineWidth', 2); % 红色虚线
ylabel('Phase (deg)', 'Color', 'r');
set(gca, 'YColor', 'r'); % 设置坐标轴颜色与曲线一致

% 辅助线：0度相位线
plot(f_plot/1e6, zeros(size(f_plot)), 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');

% 公共设置
xlabel('Frequency (MHz)');
xlim([1.5, 3.0]);
legend('|Z_{in}| (Left)', 'Phase (Right)', 'Location', 'northeast');

fprintf('绘图完成：蓝色实线为阻抗模值，红色虚线为阻抗相位。\n');

%% ========================================================================
%% 11) 绘制输出功率与平滑效率曲线（1.5 - 3.0 MHz）
%% ========================================================================

% --- 数据准备 ---
f_plot = freq_sim / 1e6;       % 频率 (MHz)
P_out  = Pout_AC;              % 输出功率 (W)
Eff    = eta_AC * 100;         % 效率 (%)

% --- 效率平滑处理 ---
% 增加平滑度：将窗口大小从 10 调至 20
Eff_smooth = smoothdata(Eff, 'gaussian', 20); 

figure('Color', 'w', 'Position', [200, 200, 800, 500]);
hold on;

% --- 左轴：输出功率 ---
yyaxis left
plot(f_plot, P_out, 'r-', 'LineWidth', 2); % 红色实线
ylabel('Output Power P_{out} (W)', 'Color', 'r');
set(gca, 'YColor', 'r'); % 轴文字设为红色
ylim([0, 800]);          % 根据需要手动设置功率量程

% --- 右轴：效率 ---
yyaxis right
plot(f_plot, Eff_smooth, 'k--', 'LineWidth', 2.5); % 黑色加粗虚线
ylabel('Efficiency (%)', 'Color', 'k');
set(gca, 'YColor', 'k'); % 轴文字设为黑色
ylim([75, 100]);         % 固定效率范围 75% - 100%

% --- 坐标轴与范围设置 ---
xlabel('Frequency (MHz)');
xlim([1.5, 3.0]);        % 修改此处：横坐标至 3.0 MHz
grid on; 
grid minor;
ax = gca;
ax.GridColor = [0.15, 0.15, 0.15]; % 强制网格为深灰色，避免变成红色

% --- 图例与边框 ---
legend('Output Power (Left)', 'Smoothed Efficiency (Right)', 'Location', 'southwest');
box on; % 补齐上沿边框

fprintf('图表更新完成：横坐标范围已限制为 1.5 - 3.0 MHz。\n');