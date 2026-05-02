clear; clc;

%% ================== 1) 用户指定工况 ==================
f_target  = 2.28e6;      % 频率 (Hz)
Vin_rms   = 20;          % 输入电压 |V_S| (RMS, V)
RL_target = 12;          % 负载电阻 (Ohm)
L_load    = 0.2e-6;        % 负载电感 (H) —— 2 uH

%% ================== 2) 只载入“已提取参数” ==================
% 文件中包含：
%   freq_exp
%   extracted_params_smooth = [ΔR11, ΔR21, CM_pF]
load('extracted_params_3D_smooth.mat', ...
     'freq_exp', 'extracted_params_smooth');

% 插值得到目标频率对应的参数
prm_target = interp1( ...
    freq_exp, ...
    extracted_params_smooth, ...
    f_target, ...
    'linear', 'extrap' );

%% ================== 3) 固定系统参数 ==================
fixed_params = struct();

fixed_params.R_L    = RL_target;
fixed_params.L_load = L_load;      % <<< 关键：感性负载
fixed_params.V_S    = Vin_rms;

fixed_params.Cs1 = 28.4e-12;
fixed_params.Cs2 = 28.4e-12;
fixed_params.Rs1 = 0.0;
fixed_params.Rs2 = 0.0;

fixed_params.L11 = 39.99e-6;
fixed_params.L12 = fixed_params.L11;
fixed_params.L21 = 34.07e-6;
fixed_params.L22 = fixed_params.L21;

% 互感
k1=0.2422; k2=0.2422; k3=0.65; k4=0.62; k5=0.234; k6=0.234;
fixed_params.M1 = -k1*sqrt(fixed_params.L11*fixed_params.L21);
fixed_params.M2 = -k2*sqrt(fixed_params.L12*fixed_params.L22);
fixed_params.M3 =  k3*sqrt(fixed_params.L11*fixed_params.L12);
fixed_params.M4 =  k4*sqrt(fixed_params.L21*fixed_params.L22);
fixed_params.M5 = -k5*sqrt(fixed_params.L11*fixed_params.L22);
fixed_params.M6 = -k6*sqrt(fixed_params.L21*fixed_params.L12);

fixed_params.RM_fixed       = 0.0;
fixed_params.R_p_stabilizer = 1e15;

%% ================== 4) 单频点仿真 ==================
[Zin, I_in, V3, V4] = run_single_freq_sim_3D_detail( ...
    f_target, prm_target, fixed_params );

%% ================== 5) 计算功率（严格物理定义） ==================
omega = 2*pi*f_target;
Z_load = RL_target + 1j*omega*L_load;

IL   = (V3 - V4) / Z_load;           % 负载电流（相量）
Pin  = real(Vin_rms * conj(I_in));   % 输入有功功率
Pout = abs(IL)^2 * RL_target;        % ⚠️ 只在 R_L 上耗散
eta  = Pout / max(Pin, eps);

%% ================== 6) 打印结果 ==================
fprintf('\n===== 单点验证（R + L 负载，不使用 NN）=====\n');
fprintf('Frequency        : %.6f MHz\n', f_target/1e6);
fprintf('|V_S| (RMS)      : %.2f V\n', Vin_rms);
fprintf('R_L              : %.2f Ohm\n', RL_target);
fprintf('L_load           : %.2f uH\n', L_load*1e6);
fprintf('--------------------------------\n');
fprintf('|I_in| (RMS)     : %.6f A\n', abs(I_in));
fprintf('∠I_in            : %.2f deg\n', angle(I_in)*180/pi);
fprintf('|I_L| (RMS)      : %.6f A\n', abs(IL));
fprintf('--------------------------------\n');
fprintf('Z_in             : %.4f %+ .4fi Ohm\n', real(Zin), imag(Zin));
fprintf('|Z_in|           : %.4f Ohm\n', abs(Zin));
fprintf('--------------------------------\n');
fprintf('Pin              : %.4f W\n', Pin);
fprintf('Pout             : %.4f W\n', Pout);
fprintf('Efficiency η     : %.2f %%\n', eta*100);
fprintf('==============================================\n');

%% =======================================================================
%                         辅助函数
% =======================================================================
function [Z_out, I_in, V3, V4] = run_single_freq_sim_3D_detail(f, params_scaled, fixed_params)

    % -------- 参数提取 --------
    delta_R11 = double(params_scaled(1));
    delta_R21 = double(params_scaled(2));
    CM        = double(params_scaled(3))*1e-12;
    omega     = 2*pi*f;

    % -------- 固定参数 --------
    R_L    = fixed_params.R_L;
    L_load = fixed_params.L_load;
    RM     = fixed_params.RM_fixed;
    Cs1    = fixed_params.Cs1;   Cs2 = fixed_params.Cs2;
    Rs1    = fixed_params.Rs1;   Rs2 = fixed_params.Rs2;
    L11    = fixed_params.L11;   L12 = fixed_params.L12;
    L21    = fixed_params.L21;   L22 = fixed_params.L22;
    M1     = fixed_params.M1;    M2  = fixed_params.M2;
    M3     = fixed_params.M3;    M4  = fixed_params.M4;
    M5     = fixed_params.M5;    M6  = fixed_params.M6;
    Vs     = fixed_params.V_S;

    % -------- 频变电阻 --------
    R11_base = 8.1969e-9*f^1.318 + 2.00;
    R21_base = 5.9395e-8*f^1.181 + 1.66;
    R11_f = max(R11_base + delta_R11, 1e-6);
    R21_f = max(R21_base + delta_R21, 1e-6);
    R12_f = R11_f;  R22_f = R21_f;

    % -------- 节点编号 --------
    n1=2; n2=3; n3=4; n4=5; n5=6; n6=7;

    % -------- 导纳矩阵 --------
    Y_full = zeros(8,8,'like',1j);

    % === 负载：R + jωL ===
    Z_load = R_L + 1j*omega*L_load;
    Y_load = 1 / Z_load;
    Y_full(n3,n3) = Y_full(n3,n3) + Y_load;
    Y_full(n4,n4) = Y_full(n4,n4) + Y_load;
    Y_full(n3,n4) = Y_full(n3,n4) - Y_load;
    Y_full(n4,n3) = Y_full(n4,n3) - Y_load;

    % === 串联补偿电容 ===
    Y_Cs1 = 1 / (Rs1 + 1/(1j*omega*Cs1));
    Y_Cs2 = 1 / (Rs2 + 1/(1j*omega*Cs2));
    Y_full(n1,n1)=Y_full(n1,n1)+Y_Cs1; Y_full(n6,n6)=Y_full(n6,n6)+Y_Cs1;
    Y_full(n1,n6)=Y_full(n1,n6)-Y_Cs1; Y_full(n6,n1)=Y_full(n6,n1)-Y_Cs1;
    Y_full(n2,n2)=Y_full(n2,n2)+Y_Cs2; Y_full(n5,n5)=Y_full(n5,n5)+Y_Cs2;
    Y_full(n2,n5)=Y_full(n2,n5)-Y_Cs2; Y_full(n5,n2)=Y_full(n5,n2)-Y_Cs2;

    % === 互容 ===
    Y_CM = 1/(RM + 1/(1j*omega*CM));
    Y_full(n1,n1)=Y_full(n1,n1)+Y_CM; Y_full(n2,n2)=Y_full(n2,n2)+Y_CM;
    Y_full(n1,n2)=Y_full(n1,n2)-Y_CM; Y_full(n2,n1)=Y_full(n2,n1)-Y_CM;
    Y_full(n6,n6)=Y_full(n6,n6)+Y_CM; Y_full(n5,n5)=Y_full(n5,n5)+Y_CM;
    Y_full(n6,n5)=Y_full(n6,n5)-Y_CM; Y_full(n5,n6)=Y_full(n5,n6)-Y_CM;

    % === 耦合线圈 ===
    R_diag = diag([R11_f,R21_f,R12_f,R22_f]);
    L_mat  = [L11,M1,M3,M5; M1,L21,M6,M4; M3,M6,L12,M2; M5,M4,M2,L22];
    Y_branch = (R_diag + 1j*omega*L_mat) \ eye(4);

    C = zeros(8,4);
    C(1,1)=1; C(2,1)=-1; C(8,2)=1; C(7,2)=-1;
    C(3,3)=1; C(4,3)=-1; C(6,4)=1; C(5,4)=-1;
    Y_full = Y_full + C*Y_branch*C.';

    % === 求解 ===
    k_idx=[2,3,4,5,6,7]; s_idx=[1,8];
    V_s=[Vs;0];
    V_k = Y_full(k_idx,k_idx)\(-Y_full(k_idx,s_idx)*V_s);
    I_s = Y_full(s_idx,k_idx)*V_k + Y_full(s_idx,s_idx)*V_s;

    I_in = I_s(1);
    Z_out = Vs / I_in;

    V_all = zeros(8,1,'like',1j);
    V_all(k_idx)=V_k; V_all(s_idx)=V_s;
    V3 = V_all(n3); V4 = V_all(n4);
end
