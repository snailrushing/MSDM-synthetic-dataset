function SIMtankerDeckMotion()
% SIMTANKERDECKMOTION  Tanker 4-DOF deck motion simulation and CSV export.
%
% Scenario:
%   Vessel  : Tanker (WAMIT, tanker.mat)
%   Sea     : Sea state 6 (Hs = 5 m, JONSWAP, w0 = 0.5 rad/s)
%   Waves   : Oblique seas, beta_wave = 135 deg (0 = following, 180 = head)
%   Motion  : Fixed heading 0 deg, 30 kn (15.43 m/s) forward speed
%   Deck pt : r_p = [-123, 0, -5] (FSD, from CO: stern, centerline, 5 m up)
%
% Outputs:
%   One CSV per independent wave realization, split by realization into
%   train/val/test subfolders under cfg.data.outputDir (70/15/15 by default):
%       generated_data/train/tanker_deck_motion_001.csv ...
%       generated_data/val/tanker_deck_motion_071.csv ...
%       generated_data/test/tanker_deck_motion_086.csv ...
%   Each CSV has columns
%       time_s, heave_m, heave_vel_mps, roll_deg, pitch_deg
%   Heave and heave velocity are positive UP; roll/pitch in degrees.
%
%   Each realization differs only by the random wave phases (and randomized
%   component directions/frequencies), driven by rng(seed) where
%       seed = cfg.data.seedStart + (run - 1).
%
% Notes:
%   - The heading autopilot holds heading but does not prevent lateral drift
%     (beam seas push the ship sideways), which is expected for "fixed heading".
%   - tauX = 2.4e6 N balances the ITTC surge resistance at ~30 kn. If the
%     printed mean speed is off target, adjust tauX and rerun.
%
% Reference: Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and
% Motion Control, 3rd ed., John Wiley & Sons Ltd.

clear headingAutopilot;   % clear persistent autopilot state
close all;

% -------------------------------------------------------------------------
% 1. Load vessel
% -------------------------------------------------------------------------
matFile = 'tanker.mat';
matPath = resolveHydroVesselFile(matFile);
data = load(matPath, 'vessel');
vessel = data.vessel;
vessel.matFile = matFile;

% -------------------------------------------------------------------------
% 2. Configuration 参数配置
% -------------------------------------------------------------------------
% 波浪谱类型（可切换）：
%   'Modified PM'  - ITTC 修正 Pierson-Moskowitz 谱（充分成长海浪）
%   'JONSWAP'      - JONSWAP 谱（风区受限、有尖峰，gamma 可调）
%   'Torsethaugen' - 双峰谱（涌浪 + 风浪，w0<=0.6 时退化为单峰）
cfg.environment.spectrumType      = 'JONSWAP';    % 波浪谱类型
cfg.environment.spreadingFlag     = 0;            % 方向扩散标志：0=长峰波（单向斜浪），1=短峰波（方向谱）
cfg.environment.Hs                = 3.0;          % 有义波高（m），5 m 对应 6 级海况,3m对应5级
cfg.environment.w0                = 0.63;         % 谱峰频率（rad/s），0.63对应5级海况
cfg.environment.gamma             = 3.3;          % JONSWAP 谱峰升高因子（其他谱忽略此参数）
cfg.environment.beta_wave         = deg2rad(135); % 浪向角（rad），135° 斜浪（0=随浪，180=顶浪）
cfg.environment.Vc                = 0;            % 海流速度（m/s），0 表示无海流
cfg.environment.betaVc            = deg2rad(0);   % 海流方向（rad）
cfg.environment.numFreqIntervals  = 100;          % 频率离散区间数（决定波浪分量的数量）
cfg.environment.numDirections     = 24;           % 方向离散数（长峰波 spreadingFlag=0 时不起作用）

cfg.simulation.T_final = 300;      % 仿真时长（s）
cfg.simulation.h       = 0.05;     % RK4 积分步长（s），也是数据采样间隔
cfg.RAO_update_period  = 0.1;      % 波浪力更新周期（s），即 10 Hz 更新

cfg.initial.eta = zeros(6,1);                % 初始位置/姿态 [x,y,z,phi,theta,psi]（NED 系，全零）
cfg.initial.nu  = [15.43; 0; 0; 0; 0; 0];    % 初始速度 [u,v,w,p,q,r]，15.43 m/s = 30 kn 前向航速

cfg.control.mode               = 'headingAutopilot';  % 控制模式：航向自动驾驶
cfg.control.setpointChangeTime = 50;        % 设定点切换时刻（s），此处参考航向恒定，实际无作用
cfg.control.heading.psi_ref       = 0;                    % 参考航向角（rad）
cfg.control.heading.psi_ref_after = 0;                    % 切换后参考航向角（rad），0=固定航向
cfg.control.heading.r_max         = deg2rad(2.0);         % 最大期望转艏角速率（rad/s）
cfg.control.heading.wn            = 1.0;                  % 航向闭环固有频率（rad/s）
cfg.control.heading.zeta          = 1.0;                  % 航向闭环相对阻尼比（-）
cfg.control.heading.tauX          = 2.4e6;  % 纵荡推力（N），用于平衡 30 kn 航速下的阻力

cfg.damping.kappa_126      = [0.05 0.05 0.05];  % 纵荡/横荡/艏摇的附加线性阻尼比
cfg.damping.delta_zeta_345 = [0 0.1 0];         % 垂荡/横摇/纵摇的阻尼修正量
cfg.damping.nonlinear_456  = [5 0 0];           % 横摇/纵摇/艏摇的非线性阻尼系数

% 数据生成设置：每次运行对应一个独立的波浪实现（随机相位不同）
% 随机种子取 seedStart .. seedStart+numRuns-1
cfg.data.numRuns   = 100;             % 生成轨迹数量（独立实现的条数）
cfg.data.seedStart = 1;               % 第一条实现的随机种子
cfg.data.outputDir = 'generated_data';% 输出 CSV 的文件夹
cfg.data.splitNames = {'train', 'val', 'test'};  % 数据集划分名称
cfg.data.splitRatio = [0.70, 0.15, 0.15];        % 训练/验证/测试的比例

vessel.kappa_4 = cfg.damping.nonlinear_456(1);   % 横摇非线性阻尼系数
vessel.kappa_5 = cfg.damping.nonlinear_456(2);   % 纵摇非线性阻尼系数
vessel.kappa_6 = cfg.damping.nonlinear_456(3);   % 艏摇非线性阻尼系数

% -------------------------------------------------------------------------
% 3. Wave model and power-based maneuvering model
% -------------------------------------------------------------------------
% Validate spectrum type and, for JONSWAP, the DNV recommended peak range.
validSpectra = {'Modified PM', 'JONSWAP', 'Torsethaugen'};
if ~any(strcmp(cfg.environment.spectrumType, validSpectra))
    error('spectrumType must be ''Modified PM'', ''JONSWAP'' or ''Torsethaugen''.');
end
if strcmp(cfg.environment.spectrumType, 'JONSWAP')
    w0sqrtHs = cfg.environment.w0 * sqrt(cfg.environment.Hs);
    if w0sqrtHs < 1.25 || w0sqrtHs > 1.75
        warning('JONSWAP: DNV recommends 1.25 <= w0*sqrt(Hs) <= 1.75 (current = %.3f).', ...
            w0sqrtHs);
    end
end

% The maneuvering model is deterministic, so it is computed once outside the
% realization loop.
vessel = computeManeuveringModel(vessel, cfg.environment.w0, ...
    cfg.damping.kappa_126, cfg.damping.delta_zeta_345, 0);

% -------------------------------------------------------------------------
% 4. Output folders (train/val/test split by realization)
% -------------------------------------------------------------------------
% Every realization is independent (i.i.d. random phases), so a fixed
% contiguous assignment run 1..nTrain -> train, etc. is a valid random split.
% Splitting by realization (not by time window) prevents train/test leakage.
nSplit = numel(cfg.data.splitNames);
splitCounts = round(cfg.data.numRuns * cfg.data.splitRatio);
splitCounts(nSplit) = cfg.data.numRuns - sum(splitCounts(1:nSplit-1));
splitEdges  = cumsum(splitCounts);   % last run index of each split

for s = 1:nSplit
    splitDir = fullfile(cfg.data.outputDir, cfg.data.splitNames{s});
    if ~exist(splitDir, 'dir')
        mkdir(splitDir);
    end
end

r_p = [-123; 0; -5];   % FSD lever arm from CO to deck point

% -------------------------------------------------------------------------
% 5. Realization loop: one simulation + one CSV per random seed
% -------------------------------------------------------------------------
for run = 1:cfg.data.numRuns

    clear headingAutopilot;          % reset persistent autopilot state
    seed = cfg.data.seedStart + run - 1;
    rng(seed);                       % distinct random wave realization

    environment = cfg.environment;   % fresh copy each run
    [environment, waveModel] = waveInitialization(vessel, environment);

    % Main time loop (replicates SIMhydroVessel with heading autopilot)
    t = 0:cfg.simulation.h:cfg.simulation.T_final;
    nTimeSteps = length(t);
    nextRAOtime = 0;

    simdata = zeros(nTimeSteps, 25);
    x = [cfg.initial.nu; cfg.initial.eta];

    for i = 1:nTimeSteps

        nu  = x(1:6);
        eta = x(7:12);

        if t(i) > cfg.control.setpointChangeTime
            cfg.control.heading.psi_ref = cfg.control.heading.psi_ref_after;
        end

        % Heading autopilot PID (yaw moment)
        tauN = headingAutopilot(eta(6), nu(6), cfg.control.heading.psi_ref, ...
            cfg.control.heading.r_max, vessel.M(6,6), ...
            cfg.control.heading.wn, cfg.control.heading.zeta, cfg.simulation.h);
        tau = [cfg.control.heading.tauX; 0; 0; 0; 0; tauN];

        % 6-DOF first-order wave forces (computed every 0.1 s)
        if t(i) >= nextRAOtime
            U = sqrt(nu(1)^2 + nu(2)^2);
            [tau_wave1, waveElevation] = waveForceRAO(t(i), waveModel, ...
                U, eta(6), environment.beta_wave);
            nextRAOtime = nextRAOtime + cfg.RAO_update_period;
        end

        simdata(i,:) = [eta; nu; tau; tau_wave1; waveElevation];

        % RK4 integration
        tau_ext = tau + tau_wave1;
        x = rk4(@hydroVessel, cfg.simulation.h, x, tau_ext, vessel, ...
            environment.Vc, environment.betaVc);

    end

    % Deck-point motion at r_p (4 DOF)
    heave     = zeros(nTimeSteps,1);  % heave (m, up positive)
    heave_vel = zeros(nTimeSteps,1);  % heave velocity (m/s, up positive)
    roll      = zeros(nTimeSteps,1);  % roll (rad)
    pitch     = zeros(nTimeSteps,1);  % pitch (rad)

    for i = 1:nTimeSteps
        eta = simdata(i,1:6)';
        nu  = simdata(i,7:12)';
        phi = eta(4); theta = eta(5); psi = eta(6);

        R = Rzyx(phi, theta, psi);             % body -> NED rotation

        % Point position in NED: p_NED = eta_pos + R*r_p  (z down positive)
        p_NED = eta(1:3) + R*r_p;
        heave(i) = -(p_NED(3) - r_p(3));       % heave rel. to t=0, up positive

        % Point velocity: v_p = nu1 + omega x r_p (body), then rotate to NED
        v_p_body = nu(1:3) + Smtrx(nu(4:6))*r_p;
        v_p_NED  = R*v_p_body;
        heave_vel(i) = -v_p_NED(3);            % up positive

        roll(i)  = phi;
        pitch(i) = theta;
    end

    % Save CSV for this realization into its train/val/test split folder
    splitIdx = find(run <= splitEdges, 1, 'first');
    splitDir = fullfile(cfg.data.outputDir, cfg.data.splitNames{splitIdx});
    out = [t(:), heave, heave_vel, rad2deg(roll), rad2deg(pitch)];
    fname = fullfile(splitDir, sprintf('tanker_deck_motion_%03d.csv', run+100));
    fid = fopen(fname, 'w');
    fprintf(fid, 'time_s,heave_m,heave_vel_mps,roll_deg,pitch_deg\n');
    for i = 1:nTimeSteps
        fprintf(fid, '%.4f,%.6f,%.6f,%.6f,%.6f\n', out(i,:));
    end
    fclose(fid);

    U = sqrt(simdata(:,7).^2 + simdata(:,8).^2);   % speed (m/s)
    fprintf('Saved %s (seed %d, %d samples, mean speed %.2f kn)\n', ...
        fname, seed, nTimeSteps, mean(U(round(nTimeSteps/2):end))/0.51444);

end

fprintf('Done: %d realizations written to %s/\n', cfg.data.numRuns, cfg.data.outputDir);
for s = 1:nSplit
    fprintf('  %-5s: %d files\n', cfg.data.splitNames{s}, splitCounts(s));
end

% -------------------------------------------------------------------------
% 6. Plot (last realization)
% -------------------------------------------------------------------------
figure;
subplot(4,1,1); plot(t, heave);     ylabel('Heave (m)');        grid on;
subplot(4,1,2); plot(t, heave_vel); ylabel('Heave vel (m/s)');  grid on;
subplot(4,1,3); plot(t, rad2deg(roll));  ylabel('Roll (deg)');  grid on;
subplot(4,1,4); plot(t, rad2deg(pitch)); ylabel('Pitch (deg)'); grid on;
xlabel('Time (s)');
if ~isoctave
    sgtitle('Tanker deck point r_p = [-123, 0, -5]: SS6, oblique seas (135^\circ), 30 kn');
end

end
