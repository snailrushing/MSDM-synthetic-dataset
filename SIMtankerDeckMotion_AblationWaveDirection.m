function SIMtankerDeckMotion_AblationWaveDirection(varargin)
% SIMTANKERDECKMOTION_ABLATIONWAVEDIRECTION  Tanker deck-motion ablation set.
%
% Generates a tanker deck-motion ablation dataset for 3 sea states (SS 4/5/6)
% and 8 wave directions offset by 22.5 deg from the directions used in
% SIMtankerDeckMotion_AllStates (0:45:315). The 8 held-out directions are
%   22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5 deg,
% giving 3 x 8 = 24 operating conditions. Each condition is simulated once
% (single random wave-phase realization), producing 24 CSV files written to
% generated_data_ablation_wave_direction/.
%
% Scenario:
%   Vessel  : Tanker (WAMIT, tanker.mat)
%   Seas    : Sea states 4, 5, 6 with JONSWAP (Hs, w0) = (2, 0.79), (3, 0.63),
%             (5, 0.5) m / rad/s (gamma = 3.3)
%   Waves   : 8 directions beta_wave = 22.5:45:337.5 deg (0 = following,
%             180 = head)
%   Motion  : Fixed heading 0 deg, 30 kn (15.43 m/s) forward speed
%   Deck pt : r_p = [-123, 0, -5] (FSD, from CO: stern, centerline, 5 m up)
%
% Outputs (one CSV per condition):
%   tanker_seastate_A_angle_B.csv
%     A : sea-state level 4, 5 or 6
%     B : wave direction 22.5, 67.5, ..., 337.5 deg
%
%   Each CSV has the columns
%       time_s, heave_m, heave_vel_mps, roll_deg, pitch_deg,
%       sea_state_norm, sin_beta_wave, cos_beta_wave
%   where heave and heave velocity are positive UP, roll/pitch are in degrees,
%   sea_state_norm = sea-state level / 15, and sin/cos are taken of the wave
%   direction angle beta_wave (rad).
%
% Optional name/value overrides (mainly for quick testing):
%   'seaStates'      1x3 integer vector (default [4 5 6])
%   'Hs'             1x3 vector, significant wave heights (default [2 3 5])
%   'w0'             1x3 vector, peak frequencies (default [0.79 0.63 0.5])
%   'waveDirections' vector of directions in deg (default 22.5:45:337.5)
%   'outputDir'      output root folder (default 'generated_data_ablation_wave_direction')
%   'seedStart'      first random seed (default 1)
%
% Reference: Fossen, T. I. (2027). Handbook of Marine Craft Hydrodynamics and
% Motion Control, 3rd ed., John Wiley & Sons Ltd.

clear headingAutopilot;   % clear persistent autopilot state
close all;

tStart = tic;

% -------------------------------------------------------------------------
% 0. Configuration defaults, overridable by name/value arguments
% -------------------------------------------------------------------------
seaStates      = [4 5 6];              % sea-state levels
Hs             = [2 3 5];              % significant wave height (m)
w0             = [0.79 0.63 0.5];      % spectral peak frequency (rad/s)
waveDirections = 22.5:45:337.5;        % held-out wave directions (deg, 0 = following)
outputDir      = 'generated_data_ablation_wave_direction';
seedStart      = 1;

for k = 1:2:numel(varargin)
    if k+1 > numel(varargin)
        error('SIMtankerDeckMotion_AblationWaveDirection:OptionValue', ...
            'Option "%s" is missing its value.', varargin{k});
    end
    switch lower(varargin{k})
        case 'seastates',       seaStates      = varargin{k+1};
        case 'hs',              Hs             = varargin{k+1};
        case 'w0',              w0             = varargin{k+1};
        case 'wavedirections',  waveDirections = varargin{k+1};
        case 'outputdir',       outputDir      = varargin{k+1};
        case 'seedstart',       seedStart      = varargin{k+1};
        otherwise
            error('SIMtankerDeckMotion_AblationWaveDirection:UnknownOption', ...
                'Unknown option "%s".', varargin{k});
    end
end

if numel(seaStates) ~= numel(Hs) || numel(seaStates) ~= numel(w0)
    error('seaStates, Hs and w0 must have the same number of elements.');
end

% -------------------------------------------------------------------------
% 1. Load vessel
% -------------------------------------------------------------------------
matFile = 'tanker.mat';
matPath = resolveHydroVesselFile(matFile);
data = load(matPath, 'vessel');
vesselBase = data.vessel;
vesselBase.matFile = matFile;

% -------------------------------------------------------------------------
% 2. Simulation configuration
% -------------------------------------------------------------------------
cfg.environment.spectrumType      = 'JONSWAP';    % wave spectrum type
cfg.environment.spreadingFlag     = 0;            % 0 = long-crested (unidirectional)
cfg.environment.gamma             = 3.3;          % JONSWAP peakedness factor
cfg.environment.Vc                = 0;            % current speed (m/s)
cfg.environment.betaVc            = deg2rad(0);   % current direction (rad)
cfg.environment.numFreqIntervals  = 100;          % number of wave components
cfg.environment.numDirections     = 24;           % ignored for long-crested seas

cfg.simulation.T_final = 300;      % simulation length (s)
cfg.simulation.h       = 0.04;     % RK4 step = 1/25 s -> 25 Hz sampling
cfg.RAO_update_period  = 0.1;      % wave-force update period (s)

U_kn = 30;                          % target forward speed (kn)
cfg.initial.eta = zeros(6,1);
cfg.initial.nu  = [U_kn*0.514444; 0; 0; 0; 0; 0];   % 30 kn in m/s

cfg.control.mode               = 'headingAutopilot';
cfg.control.setpointChangeTime = 50;        % reference-heading switch time (unused)
cfg.control.heading.psi_ref       = 0;      % reference heading (rad)
cfg.control.heading.psi_ref_after = 0;      % 0 = fixed heading
cfg.control.heading.r_max         = deg2rad(2.0);
cfg.control.heading.wn            = 1.0;
cfg.control.heading.zeta          = 1.0;
cfg.control.heading.tauX          = 2.4e6;   % surge thrust (N), balances the ITTC
                                             % resistance at ~30 kn

cfg.damping.kappa_126      = [0.05 0.05 0.05];  % surge/sway/yaw viscous damping ratio
cfg.damping.delta_zeta_345 = [0 0.1 0];         % heave/roll/pitch damping-ratio increment
cfg.damping.nonlinear_456  = [5 0 0];           % roll/pitch/yaw nonlinear damping coeffs

vesselBase.kappa_4 = cfg.damping.nonlinear_456(1);
vesselBase.kappa_5 = cfg.damping.nonlinear_456(2);
vesselBase.kappa_6 = cfg.damping.nonlinear_456(3);

% -------------------------------------------------------------------------
% 3. Output folder
% -------------------------------------------------------------------------
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

r_p = [-123; 0; -5];   % FSD lever arm from CO to deck point

numSeaStates  = numel(seaStates);
numDirections = numel(waveDirections);
seedCounter   = seedStart;

% -------------------------------------------------------------------------
% 4. Double loop: sea state -> wave direction (one realization each)
% -------------------------------------------------------------------------
for sIdx = 1:numSeaStates
    seaState = seaStates(sIdx);
    Hs_s = Hs(sIdx);
    w0_s  = w0(sIdx);

    % Informational JONSWAP peak-range check (DNV).
    w0sqrtHs = w0_s * sqrt(Hs_s);
    if w0sqrtHs < 1.25 || w0sqrtHs > 1.75
        fprintf('Note: SS%d JONSWAP w0*sqrt(Hs) = %.3f (DNV range is [1.25, 1.75]).\n', ...
            seaState, w0sqrtHs);
    end

    % The power-based maneuvering model depends only on the spectral peak
    % frequency, so it is computed once per sea state (3 times in total).
    vessel = computeManeuveringModel(vesselBase, w0_s, ...
        cfg.damping.kappa_126, cfg.damping.delta_zeta_345, 0);

    cfg.environment.Hs = Hs_s;
    cfg.environment.w0 = w0_s;

    for dIdx = 1:numDirections
        betaDeg   = waveDirections(dIdx);
        beta_wave = deg2rad(betaDeg);
        cfg.environment.beta_wave = beta_wave;

        % Constant per-condition output columns.
        seaStateNorm = seaState / 15;
        sinBeta = sin(beta_wave);
        cosBeta = cos(beta_wave);

        clear headingAutopilot;      % reset persistent autopilot state
        rng(seedCounter);            % distinct random wave realization
        seedCounter = seedCounter + 1;

        environment = cfg.environment;   % fresh copy
        [environment, waveModel] = waveInitialization(vessel, environment);

        % Main time loop (replicates SIMhydroVessel with heading autopilot).
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

            % Heading autopilot PID (yaw moment).
            tauN = headingAutopilot(eta(6), nu(6), cfg.control.heading.psi_ref, ...
                cfg.control.heading.r_max, vessel.M(6,6), ...
                cfg.control.heading.wn, cfg.control.heading.zeta, cfg.simulation.h);
            tau = [cfg.control.heading.tauX; 0; 0; 0; 0; tauN];

            % 6-DOF first-order wave forces (computed every 0.1 s).
            if t(i) >= nextRAOtime
                U = sqrt(nu(1)^2 + nu(2)^2);
                [tau_wave1, waveElevation] = waveForceRAO(t(i), waveModel, ...
                    U, eta(6), environment.beta_wave);
                nextRAOtime = nextRAOtime + cfg.RAO_update_period;
            end

            simdata(i,:) = [eta; nu; tau; tau_wave1; waveElevation];

            % RK4 integration.
            tau_ext = tau + tau_wave1;
            x = rk4(@hydroVessel, cfg.simulation.h, x, tau_ext, vessel, ...
                environment.Vc, environment.betaVc);
        end

        % Deck-point motion at r_p (4 DOF).
        heave     = zeros(nTimeSteps,1);  % heave (m, up positive)
        heave_vel = zeros(nTimeSteps,1);  % heave velocity (m/s, up positive)
        roll      = zeros(nTimeSteps,1);  % roll (rad)
        pitch     = zeros(nTimeSteps,1);  % pitch (rad)

        for i = 1:nTimeSteps
            eta = simdata(i,1:6)';
            nu  = simdata(i,7:12)';
            phi = eta(4); theta = eta(5); psi = eta(6);

            R = Rzyx(phi, theta, psi);           % body -> NED rotation

            % Point position in NED: p_NED = eta_pos + R*r_p (z down positive).
            p_NED = eta(1:3) + R*r_p;
            heave(i) = -(p_NED(3) - r_p(3));     % heave rel. to t=0, up positive

            % Point velocity: v_p = nu1 + omega x r_p (body), then to NED.
            v_p_body = nu(1:3) + Smtrx(nu(4:6))*r_p;
            v_p_NED  = R*v_p_body;
            heave_vel(i) = -v_p_NED(3);          % up positive

            roll(i)  = phi;
            pitch(i) = theta;
        end

        out = [t(:), heave, heave_vel, rad2deg(roll), rad2deg(pitch), ...
               seaStateNorm*ones(nTimeSteps,1), ...
               sinBeta*ones(nTimeSteps,1), cosBeta*ones(nTimeSteps,1)];
        angleStr = sprintf('%g', betaDeg);
        fname = fullfile(outputDir, sprintf('tanker_seastate_%d_angle_%s.csv', ...
            seaState, angleStr));
        fid = fopen(fname, 'w');
        fprintf(fid, ['time_s,heave_m,heave_vel_mps,roll_deg,pitch_deg,' ...
            'sea_state_norm,sin_beta_wave,cos_beta_wave\n']);
        fprintf(fid, '%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', out');
        fclose(fid);

        U = sqrt(simdata(:,7).^2 + simdata(:,8).^2);   % speed (m/s)
        meanSpeedKn = mean(U(round(nTimeSteps/2):end)) / 0.514444;
        fprintf('Done SS%d angle %s deg (mean speed %.2f kn)\n', ...
            seaState, angleStr, meanSpeedKn);
    end
end

fprintf('Done: %d conditions written to %s/ (%.1f s)\n', ...
    numSeaStates*numDirections, outputDir, toc(tStart));

end
