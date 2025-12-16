clear; clc; close all;
fc = 28e9;              % 28 GHz
c = 3e8; 
lambda = c/fc; 
d = lambda/2;           % Antenna spacing
N = 64;                 % 64-element Phased Array
Pt_dbm = 20;            % Transmit power 20 dBm
Pt = 10^((Pt_dbm-30)/10);
noise_pwr = 10^((-90-30)/10);

% simulated target motion
time_steps = 50; 
velocity_angular = 0.5; % Degrees per step (User moving)


% Randomize Direct Path (LOS) between -30 and 30 degrees
ang_los_true = -30 + 60*rand(); 

% Randomize Reflected Path (NLOS) to be distinct from LOS

offset = 20 + 40*rand(); 
if rand() > 0.5, ang_ref_true = ang_los_true + offset;
else, ang_ref_true = ang_los_true - offset; end

% Randomize Complex Gains
% LOS is strong (normalized 1.0), Reflected is random 0.3 to 0.7 magnitude
gain_los = 1.0 * exp(1j * 2*pi*rand());
gain_ref = (0.3 + 0.4*rand()) * exp(1j * 2*pi*rand());

fprintf('--- Initial Random Channel here :::::::::: ---\n');
fprintf('LOS Angle: %.1f deg | Ref Angle: %.1f deg\n', ang_los_true, ang_ref_true);
fprintf('------------------------------\n');


est_ang_los = ang_los_true;
est_ang_ref = ang_ref_true;


w_current = get_constructive_weights(N, d, lambda, est_ang_los, est_ang_ref, gain_los, gain_ref);
p_ref_los = abs(gain_los)^2; 
p_ref_ref = abs(gain_ref)^2; 

% tracking loop starts below
history_true_los = zeros(1, time_steps);
history_est_los  = zeros(1, time_steps);
snr_tracking     = zeros(1, time_steps);
snr_no_tracking  = zeros(1, time_steps);

% Static weight for comparison (if we didn't track)
w_static = w_current; 

fprintf('\nStarting Tracking Loop...\n');

for t = 1:time_steps
    
    ang_los_true = ang_los_true + velocity_angular;
    ang_ref_true = ang_ref_true + velocity_angular; % rigid motion target
    
    h_true = get_channel(N, d, lambda, [ang_los_true, ang_ref_true], [gain_los, gain_ref]);
    
    
    %Measure Per-Beam Power 
    
    % Current beam pointing at where we *think* the user is
    w_probe_los = get_steering_vector(N, d, lambda, est_ang_los);
    
    % Measure power received using this slightly misaligned beam
    p_measured_los = measure_power(h_true, w_probe_los, Pt);
    
    % Calculate Power Loss Ratio
    power_ratio = p_measured_los / p_ref_los; 
    
    % Estimate Angle Shift (Inverting the Beam Pattern)

    estimated_shift_deg = invert_beam_pattern(power_ratio, N);
    

    
    test_angle_pos = est_ang_los + estimated_shift_deg;
    w_test_pos = get_steering_vector(N, d, lambda, test_angle_pos);
    p_test_pos = measure_power(h_true, w_test_pos, Pt);
    
    if p_test_pos > p_measured_los
        est_ang_los = est_ang_los + estimated_shift_deg;
        est_ang_ref = est_ang_ref + estimated_shift_deg; % Apply same shift to ref
    else
        est_ang_los = est_ang_los - estimated_shift_deg;
        est_ang_ref = est_ang_ref - estimated_shift_deg;
    end
    
    
    % Update Constructive Weights with new angles
    w_current = get_constructive_weights(N, d, lambda, est_ang_los, est_ang_ref, gain_los, gain_ref);
    
    % Calculate SNR with tracking
    % 
    sig_track = h_true' * w_current * sqrt(Pt);
    snr_tracking(t) = 10*log10(abs(sig_track)^2 / noise_pwr);
    
    % calculate SNR without tracking
    sig_static = h_true' * w_static * sqrt(Pt);
    snr_no_tracking(t) = 10*log10(abs(sig_static)^2 / noise_pwr);
    

    history_true_los(t) = ang_los_true;
    history_est_los(t)  = est_ang_los;
end

% plots 
figure;
subplot(2,1,1);
plot(1:time_steps, history_true_los, 'k-', 'LineWidth', 2); hold on;
plot(1:time_steps, history_est_los, 'r--', 'LineWidth', 2);
legend('True User Angle', 'Estimated (Tracked) Angle');
ylabel('Angle (Deg)'); title('Tracking Accuracy'); grid on;

subplot(2,1,2);
plot(1:time_steps, snr_tracking, 'r-', 'LineWidth', 2); hold on;
plot(1:time_steps, snr_no_tracking, 'b--', 'LineWidth', 1.5);
legend('Per-Beam Tracking', 'Static Beam (No Tracking)');
ylabel('SNR (dB)'); xlabel('Time Step'); 
title('SNR Performance During Motion'); grid on;
ylim([min(snr_no_tracking)-5, max(snr_tracking)+5]);



function shift_deg = invert_beam_pattern(power_ratio, N)
    % Inverts the array factor G(theta) to find misalignment theta.

    
    if power_ratio >= 1.0
        shift_deg = 0;
        return;
    end
    
    % lut for the main lobe of a 64-element array
    test_angles = 0:0.05:5; % Check 0 to 5 degrees misalignment
    
    pattern = abs(sin(N*pi*sin(deg2rad(test_angles))./2) ./ (N*sin(pi*sin(deg2rad(test_angles))./2))).^2;
    
    % Find closest match to the observed power ratio
    [~, idx] = min(abs(pattern - power_ratio));
    shift_deg = test_angles(idx);
end

function w = get_constructive_weights(N, d, lambda, ang1, ang2, g1, g2)
    % Re-calculates weights combining both paths constructively
    w1 = get_steering_vector(N, d, lambda, ang1);
    w2 = get_steering_vector(N, d, lambda, ang2);
    
    ratio = (g2/g1); 
    
    w_raw = w1 + ratio * w2;
    w = w_raw / norm(w_raw);
end

function h = get_channel(N, d, lambda, angles, gains)
    h = zeros(N, 1);
    for k = 1:length(angles)
        a_vec = get_steering_vector(N, d, lambda, angles(k));
        h = h + gains(k) * a_vec;
    end
end

function a = get_steering_vector(N, d, lambda, angle_deg)
    n = (0:N-1).';
    psi = 2*pi * d/lambda * sind(angle_deg);
    a = exp(1j * n * psi) / sqrt(N); 
end

function p = measure_power(h, w, Pt)
    signal = h' * w * sqrt(Pt);
    p = abs(signal)^2;
end
