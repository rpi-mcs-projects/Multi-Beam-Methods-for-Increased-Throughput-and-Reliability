clear all; 
close all;
clc;

%% Phased array setup
isA = phased.CosineAntennaElement( ...
    'CosinePower',[1 0], ...
    'FrequencyRange',[1e9 50e9] ...
);

Nx = 8; Ny = 8;
max_reflect = 6;
fc = 28e9;               
c = physconst('LightSpeed');  
lam = c/fc; 
fs = 1e6;
sceneFile = "env.stl";
% sceneFile = "env_outdoor.stl";
% uncoment for outdoor
 
txPos = [ 1.4008; 1.5644; 0.10225 ];
rxPos = [ -1.5951; -1.6155; 3.8138];

% rxPos = [ -1.5951; -1.6155; 0.63309];
% uncoment for outdoor


sURA = phased.URA('Element',isA,...
    'Size',[Nx,Ny],...
    'ElementSpacing',0.45*lam, 'ArrayNormal','z');

%% 3D Env Simulation Setup

tx = txsite('cartesian', ...
    'Antenna', sURA, ...
    'AntennaPosition', txPos, ...
    'TransmitterFrequency', fc);

rx = rxsite('cartesian', ...
    'AntennaPosition', rxPos);

sv = siteviewer(SceneModel=sceneFile); 
show(tx);
show(rx);

pm = propagationModel("raytracing", ...
    "CoordinateSystem", "cartesian", ... 
    "Method", "sbr", ...
    "MaxNumReflections", max_reflect);

rays = raytrace(tx, rx, pm, "Map", sceneFile);


%% Ray trace Channel Anaysis

if ~isempty(rays{1})
    plot(rays{1});
else
    warning('No rays found...');
end


numRays = length(rays{1});
H_Rays = zeros(1, numRays); %channels

angles = zeros(2,numRays);

svObj = phased.SteeringVector('SensorArray', sURA);

for k = 1:numRays
    currentRay = rays{1}(k);
    
    az = currentRay.AngleOfDeparture(1);
    el = currentRay.AngleOfDeparture(2);

    V_k = svObj(fc, [az; el]); 

    angles(:, k) = [az; el];

    pathGain_linear = 10^(-currentRay.PathLoss/20);
    total_phase = deg2rad(currentRay.PhaseShift); 
    P_k = pathGain_linear * exp(-1i * total_phase);

    H_Rays(k) = (V_k' * V_k) * P_k;
end


total_sims = 5;

beam_nums = 1:total_sims;

gains = zeros(total_sims,1);

gains_realtive = zeros(total_sims,1);

sv = phased.SteeringVector('SensorArray', sURA, 'PropagationSpeed', c);


for num = beam_nums
    [best_beams,max_gain] = find_multi_beam(H_Rays,rays,num);
    best_angs = zeros(2,num);

    ws = zeros(64,num);

    for i = 1:num
        best_angs(:,i) = angles(:,best_beams(i));

        w_steer = sv(fc, best_angs(:,i));  

        tx = hamming(Nx); ty = hamming(Ny);
        window2D = tx * ty.'; 
        window2D = window2D(:);
        
        %hamming taper
        windowStrength = 0.6;  
        ws(:,i) = w_steer .* (1 - windowStrength + windowStrength*window2D);
    end

    w = sum(ws, 2);

    % normlaize total transmit power
    w = w / norm(w);

    ar = phased.ArrayResponse('SensorArray', sURA, 'WeightsInputPort', true);

    A = ar(fc, angles, w);
    H_eff = A' .* H_Rays;
    
    H_practical = sum(H_eff);
    
    gains(num) = abs(H_practical);
    
    gains_relative(num) = 10 * log10((gains(num) / gains(1))^2);

end

figure;

plot(gains_relative); hold on

plot(10*log10(1:5));

title("Relative SNR gain for varing beam numbers: Outside");
ylabel("SNR gain relative to 1 beam [dB]");
xlabel("number of beams");

legend('Practical Simulated', 'Theoretical Maximum')


function [best_rays, max_gain] = find_multi_beam(H_Rays, rays, num_beams)

    numRays = length(H_Rays);
    combos = nchoosek(1:numRays, num_beams);
    max_gain = -inf;
    best_combo = [];

    for c = 1:size(combos,1)
        idx = combos(c,:);
        combined = sum(H_Rays(idx));
        gain_mag = abs(combined);

        if gain_mag > max_gain
            max_gain = gain_mag;
            best_combo = idx;
        end
    end

    best_rays = best_combo;
end