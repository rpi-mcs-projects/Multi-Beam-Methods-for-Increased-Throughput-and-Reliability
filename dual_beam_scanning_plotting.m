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

txPos = [ 1.4008; 1.5644; 0.10225 ];
rxPos = [ -1.5951; -1.6155; 3.8138];

% rxPos = [ -1.5951; -1.6155; 0.63309];

sURA = phased.URA('Element',isA,...
    'Size',[8,8],...
    'ElementSpacing',0.45*lam, 'ArrayNormal','z');

%% 3D Env Simulation Setup

tx = txsite('cartesian', ...
    'Antenna', sURA, ...
    'AntennaPosition', txPos, ...
    'TransmitterFrequency', fc);

rx = rxsite('cartesian', ...
    'AntennaPosition', rxPos);
% 
% sv = siteviewer(SceneModel=sceneFile); 
% show(tx);
% show(rx);

pm = propagationModel("raytracing", ...
    "CoordinateSystem", "cartesian", ... 
    "Method", "sbr", ...
    "MaxNumReflections", max_reflect);

rays = raytrace(tx, rx, pm, "Map", sceneFile);


%% Ray trace Channel Anaysis

if ~isempty(rays{1})
    % plot(rays{1});
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


%find channels that sum constructively
numBeams = length(H_Rays);
maxGain = -inf;
bestBeamPair = [0, 0];

for i = 1:numBeams
    for j = i+1:numBeams
        H_i = H_Rays(i);
        H_j = H_Rays(j);
        
      
        phase_i = angle(H_i);
        phase_j = angle(H_j);
        
        combined = H_i + H_j;
        currentGainMag = abs(combined);
        
        if currentGainMag > maxGain
            maxGain = currentGainMag;
            bestBeamPair = [i, j];
        end
    end
end

bestRayA = rays{1}(bestBeamPair(1));
bestRayB = rays{1}(bestBeamPair(2));

%best path
G_single = max(abs(H_Rays));

%gain for constructive sum
G_combined = maxGain;

%gain increse
gainIncrease_dB = 20 * log10(G_combined / G_single);

fprintf('\n=== Dual-Beam Gain Comparison ===\n');
fprintf('Total Rays Found: %d\n', numRays);
fprintf('Strongest Single Beam Gain (Linear Magnitude): %.4e\n', G_single);
fprintf('Combined Dual-Beam Gain (Linear Magnitude): %.4e\n', G_combined);
fprintf('Gain Improvement from Combining 2 Beams: +%.2f dB\n', gainIncrease_dB);
fprintf('Best Constructive Beams: Ray %d (A) and Ray %d (B)\n', bestBeamPair(1), bestBeamPair(2));
fprintf('------------------------------------------------\n');

% sv2 = siteviewer(SceneModel=sceneFile); 
% show(tx);
% show(rx);

% plot([bestRayA, bestRayB]);

%% Practical anglaysis

angle1 = angles(:,bestBeamPair(1));
angle2 = angles(:,bestBeamPair(2));


sv = phased.SteeringVector('SensorArray', sURA, 'PropagationSpeed', c);
w_steer = sv(fc, angle1);  

Nx = 8; Ny = 8;
tx = hamming(Nx); ty = hamming(Ny);
window2D = tx * ty.'; 
window2D = window2D(:);

%hamming taper
windowStrength = 0.65;  
w1 = w_steer .* (1 - windowStrength + windowStrength*window2D);

% steering
w_steer = sv(fc, angle2);

% hamming taper

windowStrength = 0.65; 
w2 = w_steer .* (1 - windowStrength + windowStrength*window2D);


%add both beams
w = w1 + w2;
w = w / norm(w);


azimuth_angles = -180:1:180;

elevation_angles = 0:1:90; 

%plotting lobe diagrams

figure; 
pattern(sURA, fc, azimuth_angles, elevation_angles, ...
        'Weights', w, ...
        'PropagationSpeed', c, ...
        'Type', 'powerdb', 'Normalize', true);
title('Single-Beam Pattern');

elevation_angles = 36; 

figure; 
pattern(sURA, fc, azimuth_angles, elevation_angles, ...
        'Weights', w, ...
        'PropagationSpeed', c, ...
        'Type', 'powerdb', 'Normalize', true);
title('Single-Beam Pattern (Steered to 30 deg)');


ar = phased.ArrayResponse(SensorArray=sURA, WeightsInputPort=true);

resp = ar(fc, angles, w);

power = abs(resp).^2;
powerdB = 10*log10(power / max(power));

mask = abs(angles(2,:) - 50) < 3;

pat_masked = powerdB;
pat_masked(~mask) = NaN; 


theta = deg2rad(angles(1,:));

figure;
polarplot(theta, pat_masked, 'o'); hold on

az = -180:1:180; 
el = 50 * ones(size(az)); 
ang = [az; el];

%manually poltting (equivalent to "patter()")
ar = phased.ArrayResponse(SensorArray=sURA, WeightsInputPort=true);

resp = ar(fc, ang, w);

power = abs(resp).^2;
powerdB = 10*log10(power / max(power));

theta = deg2rad(az);
polarplot(theta, powerdB);

title('Azimuth Cut');

rlim([-60 0]);
thetalim([-180 180]);


%% Theoretical multipath gain for two beams

%phased array in all directions
ar = phased.ArrayResponse('SensorArray', sURA, 'WeightsInputPort', true);

% array response applied to channel
A1 = ar(fc, angles, w);
H_eff1 = A1' .* H_Rays;  
H_practical = sum(H_eff1);
G_practical = abs(H_practical);

A2 = ar(fc, angles, w1/norm(w1)); 
H_eff2 = A2' .* H_Rays; 
H_single = sum(H_eff2);
G_single = abs(H_single);

practicalGain_dB = 20*log10(G_practical / G_single);
