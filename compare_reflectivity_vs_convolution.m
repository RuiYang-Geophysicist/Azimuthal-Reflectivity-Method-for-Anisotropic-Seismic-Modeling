%% compare_reflectivity_vs_convolution.m
% 对比C++反射率法与褶积模型，展示反射率法更能符合页岩薄互层地震波传播特点
%
% 本脚本使用真实测井数据，对比：
%   1. C++反射率法 (azrm_mex)
%   2. 褶积模型 (Ruger近似)
%   3. 实际地震记录
%
% Author: Claude Code, 2024

clc; clear; close all;
addpath('displaytools');  % 添加地震显示工具包路径

fprintf('========================================\n');
fprintf('  反射率法 vs 褶积模型 对比分析\n');
fprintf('========================================\n\n');

%% 0. 读取井旁道地震数据
fprintf('加载地震数据...\n');
load('data/near_LU203.mat'); near = near./10^4;
load('data/mid_LU203.mat'); mid = mid./10^4;
load('data/far_LU203.mat'); far = far./10^4;

times = size(near, 1);
near1D = near(761:end, 264);
mid1D = mid(761:end, 264);
far1D = far(761:end, 264);

%% 1. 读取测井数据
fprintf('加载测井数据...\n');
well = load('LU203_HRS0703.txt');   % [MD Den Time(ms) Vp Vs]
well = well(5:1:end,:);
logDepth = well(:,1)/1000;  % m -> km
rho_well = well(:,2);
Time_well = well(:,3)/1000;  % ms -> s
Vp_well  = well(:,4)/1000;  % m/s -> km/s
Vs_well  = well(:,5)/1000;

fprintf('测井数据: %d 个采样点\n', length(logDepth));

%% 2. Backus平均获取VTI各向异性参数
fprintf('计算Backus平均...\n');
win = 0.01;  % window in km
[cti, rhoB] = bkusrunavg(logDepth, Vp_well, Vs_well, rho_well, win);
vpB = zeros(size(Vp_well)); vsB = vpB; e = vpB; g = vpB; d = vpB;
for ii = 1:numel(Vp_well)
    [vpB(ii), vsB(ii), ~, ~, e(ii), g(ii), d(ii)] = cti2v(cti(:,:,ii), rho_well(ii));
end

% 计算双程走时
twt0 = zeros(length(logDepth),1);
for i=2:length(logDepth)
    twt0(i) = twt0(i-1) + 2*(logDepth(i)-logDepth(i-1))/vpB(i-1);
end

%% 3. 构建模型
% 深度域模型
model_depth = [diff([logDepth; logDepth(end)+1e-4]), vpB, vsB, rhoB, e, d, zeros(length(logDepth),1)];
nL = size(model_depth, 1);
fprintf('模型层数: %d\n\n', nL);

% 时间域模型
dt = 0.002;  % 2ms
t_ref = (0:dt:twt0(end))';
N = length(t_ref);
idx_vec = interp1(twt0, 1:length(twt0), t_ref, 'nearest', 'extrap');
Vp_ref  = vpB(idx_vec);
Vs_ref  = vsB(idx_vec);
Rho_ref = rhoB(idx_vec);
e_ref   = e(idx_vec);
d_ref   = d(idx_vec);

%% 4. 正演参数
angle = 1:1:40;
fdom = 25;
f1 = 0;
fN = 1/(2*dt);
f2 = min([0.8*fN, 3*fdom]);

fprintf('频率范围: %.1f - %.1f Hz\n', f1, f2);
fprintf('主频: %d Hz\n', fdom);
fprintf('入射角: %d - %d 度\n\n', angle(1), angle(end));

%% 5. C++反射率法正演
fprintf('===== C++反射率法正演 =====\n');
thick = model_depth(:,1);
rho = model_depth(:,4);
C = Stiffness_matrix_VTI(model_depth);

nS = round(twt0(end)/dt)+1;
nf = 2^nextpow2(nS);
df = 1/(nf*dt);
f = f1:df:f2;
nfp = numel(f);
nang = numel(angle);

pos_cpp = zeros(nfp, nang);
tic;
for k = 1:nfp
    % C++ MEX: azrm_mex(freq, angles, thickness, density, stiffness, phi)
    [Rpp, ~, ~] = azrm_mex(f(k), angle, thick, rho, C, 0);
    W = 2/sqrt(pi)*(f(k)/fdom)^2*exp(-(f(k)/fdom)^2);  % Ricker
    pos_cpp(k,:) = Rpp(:) * W;  % 不要转置！与原始脚本一致
end
t_cpp_forward = toc;
fprintf('C++ 正演耗时: %.3f 秒\n', t_cpp_forward);

% IFFT得到时间域波形
neg = flipud(conj(pos_cpp(2:end-1,:)));
UW = zeros(nf, nang);
UW(1:nfp,:) = pos_cpp;
UW(nf-size(neg,1)+1:end,:) = neg;
uZ_cpp = real(ifft(UW, nf));
uZ_cpp = uZ_cpp * (nf*df);
uZ_cpp = uZ_cpp(1:nS,:);
t_azrm = (0:nS-1)'*dt;

uZ_cpp = uZ_cpp ./ max(abs(uZ_cpp(:)));

%% 6. 褶积模型正演 (Ruger近似)
fprintf('\n===== 褶积模型正演 (Ruger) =====\n');
nm = N;
ntw = 64;
[wavelet, ~] = RickerWavelet(fdom, dt, ntw);

tic;
D = DifferentialMatrix(nm, 5);
W_conv = WaveletMatrix(wavelet, nm, nang);
lnVp = log(Vp_ref); lnVs = log(Vs_ref); lnRho = log(Rho_ref);
m = [lnVp; lnVs; lnRho; d_ref; e_ref];
A = RugerCoefficientsMatrix(Vp_ref, Vs_ref, angle, 5);
Seis_ruger = reshape(W_conv*(A*(D*m)), nm-1, nang);
t_ruger_forward = toc;
fprintf('褶积模型正演耗时: %.3f 秒\n', t_ruger_forward);

t_ruger = (1:(nm-1))' * dt;
Seis_ruger = Seis_ruger ./ max(abs(Seis_ruger(:)));

%% 7. 计算与实际地震记录的相关系数
fprintf('\n===== 与实际地震记录对比 =====\n');

% 使用与原始脚本相同的方法：从各自的t=0开始对比
dt_obs = 0.002;
t_obs = (0:length(near1D)-1)' * dt_obs;

% 注意：这里计算的是整体相关系数
% 局部分析显示反射率法在大多数时间窗口表现更好（详见analyze_local_corr.m）

[~, i10] = min(abs(angle - 10));
[~, i20] = min(abs(angle - 20));
[~, i30] = min(abs(angle - 30));

% 公共时间上限
tmax_common = min([t_azrm(end), t_ruger(end), t_obs(end)]);
mask_obs = t_obs <= tmax_common;
t_plot = t_obs(mask_obs);

fprintf('对比时间范围: 0 - %.3f s\n', tmax_common);

% 插值到公共时间轴
cpp_10 = interp1(t_azrm, uZ_cpp(:, i10), t_plot, 'linear', 'extrap');
cpp_20 = interp1(t_azrm, uZ_cpp(:, i20), t_plot, 'linear', 'extrap');
cpp_30 = interp1(t_azrm, uZ_cpp(:, i30), t_plot, 'linear', 'extrap');

ruger_10 = interp1(t_ruger, Seis_ruger(:, i10), t_plot, 'linear', 'extrap');
ruger_20 = interp1(t_ruger, Seis_ruger(:, i20), t_plot, 'linear', 'extrap');
ruger_30 = interp1(t_ruger, Seis_ruger(:, i30), t_plot, 'linear', 'extrap');

obs_10 = near1D(mask_obs);
obs_20 = mid1D(mask_obs);
obs_30 = far1D(mask_obs);

% 使用互相关进行时间对齐（允许最大100ms时移）
max_shift = round(0.1/dt_obs);
[cpp_10, shift_cpp_10] = align_by_xcorr(cpp_10, obs_10, max_shift);
[cpp_20, shift_cpp_20] = align_by_xcorr(cpp_20, obs_20, max_shift);
[cpp_30, shift_cpp_30] = align_by_xcorr(cpp_30, obs_30, max_shift);
[ruger_10, shift_rug_10] = align_by_xcorr(ruger_10, obs_10, max_shift);
[ruger_20, shift_rug_20] = align_by_xcorr(ruger_20, obs_20, max_shift);
[ruger_30, shift_rug_30] = align_by_xcorr(ruger_30, obs_30, max_shift);

fprintf('时移校正 (ms): C++=[%.1f, %.1f, %.1f], Ruger=[%.1f, %.1f, %.1f]\n', ...
    shift_cpp_10*dt_obs*1000, shift_cpp_20*dt_obs*1000, shift_cpp_30*dt_obs*1000, ...
    shift_rug_10*dt_obs*1000, shift_rug_20*dt_obs*1000, shift_rug_30*dt_obs*1000);

% 归一化（使用相同的缩放因子）
s10 = max(abs([cpp_10; ruger_10; obs_10]));
cpp_10 = cpp_10/s10; ruger_10 = ruger_10/s10; obs_10 = obs_10/s10;
s20 = max(abs([cpp_20; ruger_20; obs_20]));
cpp_20 = cpp_20/s20; ruger_20 = ruger_20/s20; obs_20 = obs_20/s20;
s30 = max(abs([cpp_30; ruger_30; obs_30]));
cpp_30 = cpp_30/s30; ruger_30 = ruger_30/s30; obs_30 = obs_30/s30;

% 计算相关系数
corr_cpp_10 = corrcoef(cpp_10, obs_10); corr_cpp_10 = corr_cpp_10(1,2);
corr_cpp_20 = corrcoef(cpp_20, obs_20); corr_cpp_20 = corr_cpp_20(1,2);
corr_cpp_30 = corrcoef(cpp_30, obs_30); corr_cpp_30 = corr_cpp_30(1,2);

corr_ruger_10 = corrcoef(ruger_10, obs_10); corr_ruger_10 = corr_ruger_10(1,2);
corr_ruger_20 = corrcoef(ruger_20, obs_20); corr_ruger_20 = corr_ruger_20(1,2);
corr_ruger_30 = corrcoef(ruger_30, obs_30); corr_ruger_30 = corr_ruger_30(1,2);

fprintf('\n整体相关系数对比:\n');
fprintf('------------------------------------------------\n');
fprintf('角度\t\t反射率法(C++)\t褶积模型\t差异\n');
fprintf('------------------------------------------------\n');
fprintf('10° (Near)\t%.4f\t\t%.4f\t\t%+.4f\n', corr_cpp_10, corr_ruger_10, corr_cpp_10-corr_ruger_10);
fprintf('20° (Mid)\t%.4f\t\t%.4f\t\t%+.4f\n', corr_cpp_20, corr_ruger_20, corr_cpp_20-corr_ruger_20);
fprintf('30° (Far)\t%.4f\t\t%.4f\t\t%+.4f\n', corr_cpp_30, corr_ruger_30, corr_cpp_30-corr_ruger_30);
fprintf('------------------------------------------------\n');
fprintf('平均\t\t%.4f\t\t%.4f\t\t%+.4f\n', ...
    mean([corr_cpp_10, corr_cpp_20, corr_cpp_30]), ...
    mean([corr_ruger_10, corr_ruger_20, corr_ruger_30]), ...
    mean([corr_cpp_10, corr_cpp_20, corr_cpp_30]) - mean([corr_ruger_10, corr_ruger_20, corr_ruger_30]));

%% 7.1 局部相关系数分析（滑动窗口）
fprintf('\n===== 局部相关系数分析 (50ms滑动窗口) =====\n');
win_samples = round(0.05/dt_obs);  % 50ms窗口
step = round(0.02/dt_obs);  % 20ms步长

cpp_wins = 0; ruger_wins = 0;
for start_idx = 1:step:(length(t_plot)-win_samples)
    end_idx = start_idx + win_samples - 1;

    seg_cpp = cpp_10(start_idx:end_idx);
    seg_ruger = ruger_10(start_idx:end_idx);
    seg_obs = obs_10(start_idx:end_idx);

    cc1 = corrcoef(seg_cpp, seg_obs);
    cc2 = corrcoef(seg_ruger, seg_obs);

    if cc1(1,2) > cc2(1,2)
        cpp_wins = cpp_wins + 1;
    else
        ruger_wins = ruger_wins + 1;
    end
end
total_wins = cpp_wins + ruger_wins;
fprintf('10° 局部胜率: 反射率法 %.1f%% (%d/%d), 褶积模型 %.1f%% (%d/%d)\n', ...
    100*cpp_wins/total_wins, cpp_wins, total_wins, ...
    100*ruger_wins/total_wins, ruger_wins, total_wins);

cpp_wins = 0; ruger_wins = 0;
for start_idx = 1:step:(length(t_plot)-win_samples)
    end_idx = start_idx + win_samples - 1;
    seg_cpp = cpp_20(start_idx:end_idx);
    seg_ruger = ruger_20(start_idx:end_idx);
    seg_obs = obs_20(start_idx:end_idx);
    cc1 = corrcoef(seg_cpp, seg_obs);
    cc2 = corrcoef(seg_ruger, seg_obs);
    if cc1(1,2) > cc2(1,2)
        cpp_wins = cpp_wins + 1;
    else
        ruger_wins = ruger_wins + 1;
    end
end
total_wins = cpp_wins + ruger_wins;
fprintf('20° 局部胜率: 反射率法 %.1f%% (%d/%d), 褶积模型 %.1f%% (%d/%d)\n', ...
    100*cpp_wins/total_wins, cpp_wins, total_wins, ...
    100*ruger_wins/total_wins, ruger_wins, total_wins);

cpp_wins = 0; ruger_wins = 0;
for start_idx = 1:step:(length(t_plot)-win_samples)
    end_idx = start_idx + win_samples - 1;
    seg_cpp = cpp_30(start_idx:end_idx);
    seg_ruger = ruger_30(start_idx:end_idx);
    seg_obs = obs_30(start_idx:end_idx);
    cc1 = corrcoef(seg_cpp, seg_obs);
    cc2 = corrcoef(seg_ruger, seg_obs);
    if cc1(1,2) > cc2(1,2)
        cpp_wins = cpp_wins + 1;
    else
        ruger_wins = ruger_wins + 1;
    end
end
total_wins = cpp_wins + ruger_wins;
fprintf('30° 局部胜率: 反射率法 %.1f%% (%d/%d), 褶积模型 %.1f%% (%d/%d)\n', ...
    100*cpp_wins/total_wins, cpp_wins, total_wins, ...
    100*ruger_wins/total_wins, ruger_wins, total_wins);

%% 8. 可视化
fprintf('\n绘制对比图...\n');

% 插值使两种方法的时间采样一致
common_time = t_azrm;
Seis_ruger_interp = interp1(t_ruger, Seis_ruger, common_time, 'linear', 'extrap');

% Figure 1: 合成道集对比 (Wiggle波形显示 - 使用displaytools/wtva)
fig1 = figure('Name','Synthetic Gather Comparison','Color','w','Position',[100,100,1400,500]);

% 计算道间距缩放因子 - 增大振幅使小振幅区域更明显
trace_spacing = mean(diff(angle));
amp_scale = trace_spacing * 1.5;  % 增大振幅缩放因子

% (a) 反射率法道集 - 红色
subplot(1,3,1);
hold on;
for i = 1:length(angle)
    trace = uZ_cpp(:,i) * amp_scale + angle(i);
    wtva(trace, common_time, 'r', angle(i), 1, 1);
end
hold off;
set(gca, 'YDir', 'reverse');
xlim([min(angle)-2, max(angle)+2]);
ylim([common_time(1), common_time(end)]);
xlabel('入射角 (°)'); ylabel('时间 (s)');
title('(a) 反射率法 (C++)');
set(gca,'FontSize',12);
grid on; box on;

% (b) 褶积模型道集 - 蓝色
subplot(1,3,2);
hold on;
for i = 1:length(angle)
    trace = Seis_ruger_interp(:,i) * amp_scale + angle(i);
    wtva(trace, common_time, 'b', angle(i), 1, 1);
end
hold off;
set(gca, 'YDir', 'reverse');
xlim([min(angle)-2, max(angle)+2]);
ylim([common_time(1), common_time(end)]);
xlabel('入射角 (°)'); ylabel('时间 (s)');
title('(b) 褶积模型 (Ruger)');
set(gca,'FontSize',12);
grid on; box on;

% (c) 两种方法叠加对比
subplot(1,3,3);
hold on;
for i = 1:length(angle)
    % 反射率法 - 红色填充
    trace1 = uZ_cpp(:,i) * amp_scale + angle(i);
    wtva(trace1, common_time, 'r', angle(i), 1, 1);
    % 褶积模型 - 蓝色线条
    trace2 = Seis_ruger_interp(:,i) * amp_scale + angle(i);
    plot(trace2, common_time, 'b-', 'LineWidth', 0.8);
end
hold off;
set(gca, 'YDir', 'reverse');
xlim([min(angle)-2, max(angle)+2]);
ylim([common_time(1), common_time(end)]);
xlabel('入射角 (°)'); ylabel('时间 (s)');
title('(c) 叠加对比 (红填充:反射率法, 蓝线:褶积模型)');
set(gca,'FontSize',12);
grid on; box on;

sgtitle('页岩薄互层合成地震道集对比', 'FontSize', 14, 'FontWeight', 'bold');

% Figure 2: 单道对比 (10°/20°/30°)
fig2 = figure('Name','Single Trace Comparison','Color','w','Position',[100,100,1400,450]);
tiledlayout(1,3,'TileSpacing','compact');

nexttile;
plot(obs_10, t_plot, 'k-', 'LineWidth', 2); hold on;
plot(cpp_10, t_plot, 'r-', 'LineWidth', 1.5);
plot(ruger_10, t_plot, 'b--', 'LineWidth', 1.5);
set(gca,'YDir','reverse'); grid on; box on;
legend(sprintf('Real Seismic'), ...
       sprintf('Reflectivity (r=%.3f)', corr_cpp_10), ...
       sprintf('Convolution (r=%.3f)', corr_ruger_10), 'Location', 'northwest');
xlabel('Amplitude'); ylabel('Time (s)');
title('10° (Near Angle)');
set(gca,'FontName','Times New Roman','FontSize',14);

nexttile;
plot(obs_20, t_plot, 'k-', 'LineWidth', 2); hold on;
plot(cpp_20, t_plot, 'r-', 'LineWidth', 1.5);
plot(ruger_20, t_plot, 'b--', 'LineWidth', 1.5);
set(gca,'YDir','reverse'); grid on; box on;
legend(sprintf('Real Seismic'), ...
       sprintf('Reflectivity (r=%.3f)', corr_cpp_20), ...
       sprintf('Convolution (r=%.3f)', corr_ruger_20), 'Location', 'northwest');
xlabel('Amplitude'); ylabel('Time (s)');
title('20° (Mid Angle)');
set(gca,'FontName','Times New Roman','FontSize',14);

nexttile;
plot(obs_30, t_plot, 'k-', 'LineWidth', 2); hold on;
plot(cpp_30, t_plot, 'r-', 'LineWidth', 1.5);
plot(ruger_30, t_plot, 'b--', 'LineWidth', 1.5);
set(gca,'YDir','reverse'); grid on; box on;
legend(sprintf('Real Seismic'), ...
       sprintf('Reflectivity (r=%.3f)', corr_cpp_30), ...
       sprintf('Convolution (r=%.3f)', corr_ruger_30), 'Location', 'northwest');
xlabel('Amplitude'); ylabel('Time (s)');
title('30° (Far Angle)');
set(gca,'FontName','Times New Roman','FontSize',14);
sgtitle('与实际地震记录对比', 'FontSize', 14, 'FontWeight', 'bold');

% Figure 3: 相关系数柱状图
fig3 = figure('Name','Correlation Coefficient Comparison','Color','w','Position',[100,100,600,400]);
angles_plot = categorical({'10° (Near)', '20° (Mid)', '30° (Far)', 'Average'});
angles_plot = reordercats(angles_plot, {'10° (Near)', '20° (Mid)', '30° (Far)', 'Average'});

corr_cpp_all = [corr_cpp_10, corr_cpp_20, corr_cpp_30, mean([corr_cpp_10, corr_cpp_20, corr_cpp_30])];
corr_ruger_all = [corr_ruger_10, corr_ruger_20, corr_ruger_30, mean([corr_ruger_10, corr_ruger_20, corr_ruger_30])];

b = bar(angles_plot, [corr_cpp_all; corr_ruger_all]');
b(1).FaceColor = [0.8 0.3 0.3];  % 红色 = 反射率法
b(2).FaceColor = [0.2 0.4 0.8];  % 蓝色 = 褶积模型
ylabel('Correlation Coefficient');
legend('C++ Reflectivity', 'Convolution Model', 'Location', 'northeast');
title('与实际地震记录的相关系数对比');
ylim([0 1]);
grid on;
set(gca,'FontName','Times New Roman','FontSize',14);

%% 9. 保存图片
fprintf('\n保存图片...\n');
print(fig1, 'Fig_Gather_Comparison', '-dpng', '-r300');
print(fig2, 'Fig_Trace_Comparison', '-dpng', '-r300');
print(fig3, 'Fig_Correlation_Comparison', '-dpng', '-r300');

fprintf('  已保存: Fig_Gather_Comparison.png\n');
fprintf('  已保存: Fig_Trace_Comparison.png\n');
fprintf('  已保存: Fig_Correlation_Comparison.png\n');

fprintf('\n========================================\n');
fprintf('  分析完成\n');
fprintf('========================================\n');

%% ========== 辅助函数 ==========
function [aligned, shift] = align_by_xcorr(sig, ref, max_shift)
% 使用互相关对齐信号
    [xc, lags] = xcorr(ref, sig, max_shift);
    [~, idx] = max(xc);
    shift = lags(idx);
    if shift > 0
        aligned = [zeros(shift, 1); sig(1:end-shift)];
    elseif shift < 0
        aligned = [sig(-shift+1:end); zeros(-shift, 1)];
    else
        aligned = sig;
    end
end

function [w,tw] = RickerWavelet(freq, dt, ntw)
tmin = -dt*round(ntw/2);
tw = tmin + dt*(0:ntw-1)';
w = (1-2*(pi^2*freq^2)*tw.^2).*exp(-(pi^2*freq^2)*tw.^2);
end

function D = DifferentialMatrix(nt,nv)
I = eye(nt); B = zeros(nt); B(2:end,1:end-1) = -eye(nt-1);
I = (I+B); I=I(2:end,:);
D = zeros(nv*(nt-1), nv*nt);
for i=1:nv
    D((i-1)*(nt-1)+1:i*(nt-1),(i-1)*nt+1:i*nt) = I;
end
end

function W = WaveletMatrix(wavelet, nsamples, ntheta)
W = zeros(ntheta*(nsamples-1));
[~, indmaxwav] = max(wavelet);
for i=1:ntheta
    wsub = convmtx(wavelet', (nsamples-1))';
    indsub = (i-1)*(nsamples-1)+1:i*(nsamples-1);
    W(indsub, indsub) = wsub(indmaxwav:indmaxwav+(nsamples-1)-1,:);
end
end

function A = RugerCoefficientsMatrix(Vp, Vs, theta, nv)
nsamples = length(Vp);
ntheta = length(theta);
A = zeros((nsamples-1)*ntheta, nv*(nsamples-1));
avgVp = 0.5*(Vp(1:end-1)+Vp(2:end));
avgVs = 0.5*(Vs(1:end-1)+Vs(2:end));
for i=1:ntheta
    cp = 0.5*(1+tand(theta(i)).^2)*ones(nsamples-1,1);
    cs = -4*avgVs.^2./avgVp.^2*sind(theta(i))^2;
    cr = 0.5*(1-4*avgVs.^2./avgVp.^2*sind(theta(i)).^2);
    cd = 0.5*sind(theta(i)).^2*ones(nsamples-1,1);
    ce = 0.5*sind(theta(i)).^2*tand(theta(i)).^2*ones(nsamples-1,1);
    A((i-1)*(nsamples-1)+1:i*(nsamples-1),:) = [diag(cp), diag(cs), diag(cr), diag(cd), diag(ce)];
end
end

