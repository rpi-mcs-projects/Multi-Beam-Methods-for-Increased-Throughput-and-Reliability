import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# System Config
fc = 28e9               # Frequency: 28 GHz
c = 3e8                 
lambd = c / fc          
d = lambd / 2           # antenna spacing
N = 64                  # 64-element Phased Array (ULA)
Pt_dbm = 20             
Pt = 10**((Pt_dbm - 30) / 10)  # power conversion to watts
noise_pwr = 10**((-90 - 30) / 10) # noise floor 9-90 dbm

# Distance for Path Loss
dist = 7 # meters
# incorporate path loss
fspl_lin = (4 * np.pi * dist / lambd)**2
fspl_db = 10 * np.log10(fspl_lin)

print(f"Path Loss applied: {fspl_db:.2f} dB")

# Simulation Settings
time_steps = 120
t_vec = np.arange(time_steps)


path_losses_db = np.array([0, 6, 8, 10]) 
path_amps = 10**(-path_losses_db / 20)


init_angles = np.array([-10.0, 25.0, -45.0, 60.0])

# randomize path phases
np.random.seed(42) 
path_phases = np.random.uniform(0, 2*np.pi, 4)
path_phases[0] = 0 
init_gains = path_amps * np.exp(1j * path_phases)

# helper functions - reference in paper
def get_steering_vector(N, d, lambd, angle_deg):
    n = np.arange(N).reshape(-1, 1)
    psi = 2 * np.pi * d / lambd * np.sin(np.deg2rad(angle_deg))
    return np.exp(1j * n * psi) / np.sqrt(N)

def get_channel(N, d, lambd, angles, gains):
    h = np.zeros((N, 1), dtype=complex)
    for k in range(len(angles)):
        a = get_steering_vector(N, d, lambd, angles[k])
        h += gains[k] * a
    return h

def measure_power(h, w, Pt):
    # Received Signal = h^H * w * sqrt(Pt)
    signal_power_lossless = np.abs(np.vdot(h, w) * np.sqrt(Pt))**2
    return signal_power_lossless / fspl_lin

def invert_beam_pattern(power_ratio, N):
    if power_ratio >= 0.99: return 0.0
    if power_ratio < 0.01: return 5.0 
    
    test_angles = np.arange(0, 6.05, 0.05)
    psi = np.pi * np.sin(np.deg2rad(test_angles))
    with np.errstate(divide='ignore', invalid='ignore'):
        pattern = np.abs(np.sin(N * psi / 2) / (N * np.sin(psi / 2)))
    pattern[0] = 1.0
    
    idx = np.argmin(np.abs(pattern**2 - power_ratio))
    return test_angles[idx]

def get_array_response(w, N, d, lambd, theta_vals):
    # Calculate Array Factor for a set of angles theta_vals
    resp = np.zeros(len(theta_vals), dtype=complex)
    for i, ang in enumerate(theta_vals):
        a = get_steering_vector(N, d, lambd, ang)
        resp[i] = np.vdot(a, w) # w^H * a
    return resp

est_angles = init_angles.copy()
est_ang_los = est_angles[0] 

# tracking reference
w_los_init = get_steering_vector(N, d, lambd, est_ang_los)
h_init = get_channel(N, d, lambd, init_angles, init_gains)
p_ref_tracking = measure_power(h_init, w_los_init, Pt) # Now includes path loss

# Store data for plotting
snapshots = {}
snapshot_times = [10, 25, 55, 90] # times for different blocking stages. Seen in SNR plot

snr_mmReliable = []
snr_single = []

# main simulation loop
for t in t_vec:
    
    velocity = 0.15 
    curr_angles = init_angles + velocity * t
    curr_gains = init_gains.copy()
    
    # blockage
    blockage_state = "Clear"
    if 20 <= t <= 35: 
        curr_gains[1] *= 0.01  # Block NLOS 1
        blockage_state = "NLOS 1 Blocked"
    if 50 <= t <= 65:          # Block NLOS 2 & 3
        curr_gains[2] *= 0.01
        curr_gains[3] *= 0.01
        blockage_state = "NLOS 2&3 Blocked"
    if 80 <= t <= 100: 
        curr_gains[0] *= 0.001 # Block LOS (Deep Fade)
        blockage_state = "LOS Blocked"
        
    h_curr = get_channel(N, d, lambd, curr_angles, curr_gains)
    
    # tracking logic (steering vector calculations)
    w_track_probe = get_steering_vector(N, d, lambd, est_ang_los)
    p_track = measure_power(h_curr, w_track_probe, Pt)
    
    ratio_track = p_track / p_ref_tracking
    shift = 0
    if ratio_track > 0.1: 
        shift = invert_beam_pattern(ratio_track, N)
        w_test = get_steering_vector(N, d, lambd, est_ang_los + shift)
        if measure_power(h_curr, w_test, Pt) > p_track: direction = 1
        else: direction = -1
        shift *= direction
    
    est_angles += shift
    est_ang_los = est_angles[0] 
    
    # Optimal beam scaling
    w_mm = np.zeros((N, 1), dtype=complex)
    for k in range(4):
        vk = get_steering_vector(N, d, lambd, est_angles[k])
        w_mm += vk * np.conj(curr_gains[k]) 
        
    w_mm = w_mm / np.linalg.norm(w_mm) 
    
    # Single Beam reference
    w_s = get_steering_vector(N, d, lambd, est_angles[0])
    
    # SNR calcs
    sig_pwr_mm = measure_power(h_curr, w_mm, Pt)
    sig_pwr_s = measure_power(h_curr, w_s, Pt)
    
    snr_mmReliable.append(10 * np.log10(sig_pwr_mm / noise_pwr))
    snr_single.append(10 * np.log10(sig_pwr_s / noise_pwr))

    # data capture
    if t in snapshot_times:
        snapshots[t] = {
            'w_mm': w_mm.copy(),
            'angles': est_angles.copy(),
            'gains': curr_gains.copy(),
            'state': blockage_state
        }

# the following is plotting logic
plt.figure(figsize=(10, 6))
plt.plot(t_vec, snr_mmReliable, 'r-', linewidth=2, label='Optimal Scaled Multi-Beam')
plt.plot(t_vec, snr_single, 'b--', linewidth=2, label='Single Beam')

plt.axvspan(20, 35, color='orange', alpha=0.1, label='NLOS 1 Blocked')
plt.axvspan(50, 65, color='green', alpha=0.1, label='NLOS 2&3 Blocked')
plt.axvspan(80, 100, color='gray', alpha=0.3, label='LOS Blocked')

# Outage Threshold: paper has it as 5-6 dB
plt.axhline(6, color='k', linestyle=':', label='Outage Threshold (6 dB)')

plt.title(f"SNR with Tracking & Optimal Beam Scaling")
plt.legend(loc='lower left')
plt.ylabel('SNR (dB)')
plt.xlabel('Time Step')
plt.ylim(-20, 50) # SNR range (y axis)
plt.grid(True, alpha=0.3)
plt.show()

# polar plot
theta_plot = np.linspace(-90, 90, 360) 
fig, axs = plt.subplots(2, 2, figsize=(10, 10), subplot_kw={'projection': 'polar'})
axs = axs.flatten()

for i, t in enumerate(snapshot_times):
    data = snapshots[t]
    w_curr = data['w_mm']
    state = data['state']
    
    
    resp = get_array_response(w_curr, N, d, lambd, theta_plot)
    resp_db = 20 * np.log10(np.abs(resp) + 1e-9)
    resp_db = resp_db - np.max(resp_db) 
    
    ax = axs[i]
    ax.plot(np.deg2rad(theta_plot), resp_db, 'r-', linewidth=2)
    ax.set_title(f"T={t}: {state}", fontsize=8, pad=10)
    ax.set_ylim(-15, 0)
    ax.set_yticks([-30, -20, -10, 0])
    ax.set_theta_zero_location("N") 
    ax.set_theta_direction(-1)      # Clockwise
    
    
    current_angles = data['angles']
    current_gains = np.abs(data['gains'])
    

    for k in range(4):
        ang = current_angles[k]
        gain = current_gains[k]
        
        # Only mark significant paths (> -20dB relative to max possible)
        if gain > 1e-4: # If path is active
            ax.plot([np.deg2rad(ang), np.deg2rad(ang)], [-40, 0], 'g--', alpha=0.5)
        else:
            # Blocked path
            ax.plot([np.deg2rad(ang), np.deg2rad(ang)], [-40, 0], 'k:', alpha=0.3)

plt.tight_layout()
plt.show()
