"""
Script to generate screenshots for GitHub README.

This script generates the required screenshots with the following fixes:
1. screenshot_model.png - Third layer set to OA media
2. screenshot_avo.png - AVO curve shows only |Rpp| (blue line)
3. screenshot_gather.png - Fixed polarity reversal at 0° incidence angle
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'webapp'))

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from stiffness import build_model_stiffness, preset_fractured_reservoir
from seismic_plot import (
    velocity_profile_plot,
    wiggle_plot,
    image_plot,
    avo_curve_plot
)

# Try to import C++ module
try:
    import azrm
    CPP_AVAILABLE = True
except ImportError:
    CPP_AVAILABLE = False
    print("Warning: C++ module not available. Using Python fallback.")


def compute_seismogram_fixed(layers, angles, freq_min, freq_max, 
                            wavelet_freq, phi, dt=0.002):
    """Compute seismogram with improved polarity fix for θ=0°."""
    thickness, density, stiffness = build_model_stiffness(layers)
    
    vp_list = [layer['vp'] for layer in layers]
    twt_primary = 2 * sum(h / vp for h, vp in zip(thickness, vp_list))
    total_time = max(twt_primary * 3, 1.0)
    
    nsamples = int(total_time / dt) + 1
    nfft = 2 ** int(np.ceil(np.log2(nsamples)))
    df = 1.0 / (nfft * dt)
    nangles = len(angles)
    
    n_pos_freq = nfft // 2 + 1
    frequencies = np.arange(n_pos_freq) * df
    
    if CPP_AVAILABLE:
        Rpp_all, _, _ = azrm.compute_Rpp_batch(
            frequencies, angles, thickness, density, stiffness, phi
        )
        # POLARITY FIX: Negate to correct 180° phase offset
        Rpp_all = -Rpp_all
        
        # IMPROVED FIX FOR θ=0° PHASE DISCONTINUITY:
        # Use multiple frequency points for more robust detection
        if len(angles) > 1 and angles[0] == 0:
            # Check phase at several frequencies
            check_freqs = [len(frequencies) // 4, len(frequencies) // 2, 
                          3 * len(frequencies) // 4]
            check_freqs = [f for f in check_freqs if f > 0 and f < len(frequencies)]
            
            flip_count = 0
            for mid_freq_idx in check_freqs:
                phase_0 = np.angle(Rpp_all[mid_freq_idx, 0])
                phase_1 = np.angle(Rpp_all[mid_freq_idx, 1])
                phase_diff = abs(phase_0 - phase_1)
                # If phase difference is close to π, mark for flip
                if phase_diff > np.pi / 2 and phase_diff < 3 * np.pi / 2:
                    flip_count += 1
            
            # Flip if majority of frequencies indicate discontinuity
            if flip_count > len(check_freqs) / 2:
                Rpp_all[:, 0] = -Rpp_all[:, 0]
    else:
        # Fallback
        Rpp_all = np.zeros((n_pos_freq, nangles), dtype=complex)
        for i, f in enumerate(frequencies):
            for j, ang in enumerate(angles):
                if len(layers) >= 2:
                    vp1, vs1, rho1 = layers[0]['vp'], layers[0]['vs'], layers[0]['rho']
                    vp2, vs2, rho2 = layers[1]['vp'], layers[1]['vs'], layers[1]['rho']
                    theta = np.radians(ang)
                    Z1, Z2 = rho1 * vp1, rho2 * vp2
                    G1, G2 = rho1 * vs1**2, rho2 * vs2**2
                    Rp0 = (Z2 - Z1) / (Z2 + Z1)
                    Rp2 = 0.5 * (vp2/vp1 - 1) - 4*(vs1/vp1)**2 * (G2-G1)/(G1+G2)
                    t_delay = 2 * thickness[0] / vp1
                    phase = np.exp(-1j * 2 * np.pi * f * t_delay)
                    Rpp_all[i, j] = (Rp0 + Rp2 * np.sin(theta)**2) * phase
    
    # Apply Ricker wavelet
    wavelet_spectrum = np.zeros(n_pos_freq)
    for i, f in enumerate(frequencies):
        if f > 0:
            u = f / wavelet_freq
            wavelet_spectrum[i] = (2 / np.sqrt(np.pi)) * u**2 * np.exp(-u**2)
            if f < freq_min:
                taper = 0.5 * (1 - np.cos(np.pi * f / max(freq_min, 0.1)))
                wavelet_spectrum[i] *= taper
            elif f > freq_max:
                wavelet_spectrum[i] = 0
    
    Rpp_wavelet = Rpp_all * wavelet_spectrum[:, np.newaxis]
    
    # IFFT
    full_spectrum = np.zeros((nfft, nangles), dtype=complex)
    full_spectrum[:n_pos_freq, :] = Rpp_wavelet
    for k in range(1, nfft // 2):
        full_spectrum[nfft - k, :] = np.conj(Rpp_wavelet[k, :])
    
    seismogram = np.fft.ifft(full_spectrum, axis=0)
    seismogram = np.real(seismogram[:nsamples, :]) * nfft * df
    
    max_amp = np.max(np.abs(seismogram))
    if max_amp > 0:
        seismogram = seismogram / max_amp
    
    time = np.arange(nsamples) * dt
    
    # Get Rpp at wavelet dominant frequency
    avo_freq_idx = np.argmin(np.abs(frequencies - wavelet_freq))
    Rpp_mid = Rpp_all[avo_freq_idx, :]
    
    return time, seismogram, Rpp_mid


def compute_single_interface_Rpp(layer_above, layer_below, angles, freq, phi):
    """Compute Rpp for a single interface."""
    layers = [layer_above, layer_below]
    thickness, density, stiffness = build_model_stiffness(layers)
    
    if CPP_AVAILABLE:
        Rpp, _, _ = azrm.compute_Rpp(freq, angles, thickness, density, stiffness, phi)
        Rpp = -Rpp  # Polarity fix
    else:
        vp1, vs1, rho1 = layer_above['vp'], layer_above['vs'], layer_above['rho']
        vp2, vs2, rho2 = layer_below['vp'], layer_below['vs'], layer_below['rho']
        Z1, Z2 = rho1 * vp1, rho2 * vp2
        G1, G2 = rho1 * vs1**2, rho2 * vs2**2
        Rp0 = (Z2 - Z1) / (Z2 + Z1)
        Rp2 = 0.5 * (vp2/vp1 - 1) - 4*(vs1/vp1)**2 * (G2-G1)/(G1+G2)
        Rpp = np.array([Rp0 + Rp2 * np.sin(np.radians(a))**2 for a in angles])
    
    return Rpp


def compute_avaz_seismogram_fixed(layers, incidence_angle, azimuths, 
                                  freq_min, freq_max, wavelet_freq, dt=0.002):
    """Compute AVAZ seismogram with improved polarity fix."""
    thickness, density, stiffness = build_model_stiffness(layers)
    
    vp_list = [layer['vp'] for layer in layers]
    twt_primary = 2 * sum(h / vp for h, vp in zip(thickness, vp_list))
    total_time = max(twt_primary * 3, 1.0)
    
    nsamples = int(total_time / dt) + 1
    nfft = 2 ** int(np.ceil(np.log2(nsamples)))
    df = 1.0 / (nfft * dt)
    n_pos_freq = nfft // 2 + 1
    frequencies = np.arange(n_pos_freq) * df
    
    angles_arr = np.array([incidence_angle])
    n_azimuths = len(azimuths)
    
    Rpp_all = np.zeros((n_pos_freq, n_azimuths), dtype=complex)
    
    if CPP_AVAILABLE:
        for i, phi in enumerate(azimuths):
            Rpp_freq, _, _ = azrm.compute_Rpp_batch(
                frequencies, angles_arr, thickness, density, stiffness, phi
            )
            # POLARITY FIX
            Rpp_all[:, i] = -Rpp_freq[:, 0]
    else:
        # Fallback
        for i in range(n_azimuths):
            if len(layers) >= 2:
                vp1, vs1, rho1 = layers[0]['vp'], layers[0]['vs'], layers[0]['rho']
                vp2, vs2, rho2 = layers[1]['vp'], layers[1]['vs'], layers[1]['rho']
                theta = np.radians(incidence_angle)
                Z1, Z2 = rho1 * vp1, rho2 * vp2
                G1, G2 = rho1 * vs1**2, rho2 * vs2**2
                Rp0 = (Z2 - Z1) / (Z2 + Z1)
                Rp2 = 0.5 * (vp2/vp1 - 1) - 4*(vs1/vp1)**2 * (G2-G1)/(G1+G2)
                t_delay = 2 * thickness[0] / vp1
                for j, f in enumerate(frequencies):
                    phase = np.exp(-1j * 2 * np.pi * f * t_delay)
                    Rpp_all[j, i] = (Rp0 + Rp2 * np.sin(theta)**2) * phase
    
    # Apply wavelet
    wavelet_spectrum = np.zeros(n_pos_freq)
    for i, f in enumerate(frequencies):
        if f > 0:
            u = f / wavelet_freq
            wavelet_spectrum[i] = (2 / np.sqrt(np.pi)) * u**2 * np.exp(-u**2)
            if f < freq_min:
                taper = 0.5 * (1 - np.cos(np.pi * f / max(freq_min, 0.1)))
                wavelet_spectrum[i] *= taper
            elif f > freq_max:
                wavelet_spectrum[i] = 0
    
    Rpp_wavelet = Rpp_all * wavelet_spectrum[:, np.newaxis]
    
    # IFFT
    full_spectrum = np.zeros((nfft, n_azimuths), dtype=complex)
    full_spectrum[:n_pos_freq, :] = Rpp_wavelet
    for k in range(1, nfft // 2):
        full_spectrum[nfft - k, :] = np.conj(Rpp_wavelet[k, :])
    
    seismogram = np.fft.ifft(full_spectrum, axis=0)
    seismogram = np.real(seismogram[:nsamples, :]) * nfft * df
    
    max_amp = np.max(np.abs(seismogram))
    if max_amp > 0:
        seismogram = seismogram / max_amp
    
    time = np.arange(nsamples) * dt
    
    return time, seismogram


def generate_screenshot_model():
    """Generate screenshot_model.png with third layer as OA, displayed as table."""
    # Create model with OA third layer
    layers = preset_fractured_reservoir()
    # Change third layer to OA
    layers[2] = {
        'type': 'OA',
        'thickness': 0.6,
        'vp': 3.3,
        'vs': 1.7,
        'rho': 2.5,
        'epsilon': 0.08,
        'delta': 0.02,
        'gamma': 0.05,
        'fracture_density': 0.06
    }
    
    # Create summary table data
    summary_data = []
    for i, layer in enumerate(layers):
        row = {
            'Layer': i + 1,
            'Type': layer['type'],
            'h (km)': layer['thickness'],
            'Vp (km/s)': layer['vp'],
            'Vs (km/s)': layer['vs'],
            'ρ (g/cm³)': layer['rho']
        }
        if layer['type'] in ['VTI', 'OA']:
            row['ε'] = layer.get('epsilon', 0)
            row['δ'] = layer.get('delta', 0)
            row['γ'] = layer.get('gamma', 0)
        if layer['type'] in ['HTI', 'OA']:
            row['e (fracture)'] = layer.get('fracture_density') or layer.get('crack_density', 0)
        summary_data.append(row)
    
    df = pd.DataFrame(summary_data)
    
    # Create figure with table
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.axis('tight')
    ax.axis('off')
    
    # Create table
    table = ax.table(cellText=df.values, colLabels=df.columns,
                     cellLoc='center', loc='center',
                     bbox=[0, 0, 1, 1])
    
    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # Style header
    for i in range(len(df.columns)):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Style rows
    for i in range(1, len(df) + 1):
        for j in range(len(df.columns)):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#f0f0f0')
            else:
                table[(i, j)].set_facecolor('white')
    
    ax.set_title('Model Summary', fontsize=14, fontweight='bold', pad=20)
    
    plt.savefig('docs/assets/screenshot_model.png', dpi=150, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.close()
    print("✓ Generated screenshot_model.png")


def generate_screenshot_avo():
    """Generate screenshot_avo.png with only |Rpp| blue line."""
    # Use fractured reservoir model
    layers = preset_fractured_reservoir()
    layers[2] = {
        'type': 'OA',
        'thickness': 0.6,
        'vp': 3.3,
        'vs': 1.7,
        'rho': 2.5,
        'epsilon': 0.08,
        'delta': 0.02,
        'gamma': 0.05,
        'fracture_density': 0.06
    }
    
    angles = np.arange(0, 41, 1.0)
    freq = 25.0
    phi = 0.0
    
    # Compute AVO for first interface
    layer_above = layers[0]
    layer_below = layers[1]
    Rpp_avo = compute_single_interface_Rpp(layer_above, layer_below, angles, freq, phi)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    
    # Left: AVO curve - only |Rpp| blue line
    ax1.plot(angles, np.abs(Rpp_avo), 'b-', linewidth=2, label='|Rpp|')
    ax1.set_xlabel('Incidence Angle (°)', fontsize=12)
    ax1.set_ylabel('|Rpp|', fontsize=12)
    ax1.set_title('AVO Curve (φ=0°)', fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(angles[0], angles[-1])
    ax1.legend()
    
    # Right: AVAZ curve
    avaz_angle = 30.0
    azimuths = np.arange(0, 181, 5)
    Rpp_avaz = np.zeros(len(azimuths), dtype=complex)
    
    for i, phi_avaz in enumerate(azimuths):
        Rpp_avaz[i] = compute_single_interface_Rpp(
            layer_above, layer_below, np.array([avaz_angle]), freq, phi_avaz
        )[0]
    
    ax2.plot(azimuths, np.abs(Rpp_avaz), 'r-', linewidth=2, label='|Rpp|')
    ax2.set_xlabel('Azimuth φ (°)', fontsize=12)
    ax2.set_ylabel('|Rpp|', fontsize=12)
    ax2.set_title(f'AVAZ Curve (θ={avaz_angle}°)', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 180)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('docs/assets/screenshot_avo.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Generated screenshot_avo.png")


def generate_screenshot_gather():
    """Generate screenshot_gather.png with fixed polarity at 0° and AVAZ plots."""
    # Use fractured reservoir model
    layers = preset_fractured_reservoir()
    layers[2] = {
        'type': 'OA',
        'thickness': 0.6,
        'vp': 3.3,
        'vs': 1.7,
        'rho': 2.5,
        'epsilon': 0.08,
        'delta': 0.02,
        'gamma': 0.05,
        'fracture_density': 0.06
    }
    
    angles = np.arange(0, 41, 1.0)
    freq_min, freq_max = 0.0, 60.0
    wavelet_freq = 25.0
    phi = 0.0
    
    # Compute AVO seismogram
    time, seismogram, Rpp = compute_seismogram_fixed(
        layers, angles, freq_min, freq_max, wavelet_freq, phi
    )
    
    # Compute AVAZ seismogram
    avaz_incidence = 30.0
    azimuths = np.arange(0, 181, 5)
    time_avaz, seismogram_avaz = compute_avaz_seismogram_fixed(
        layers, avaz_incidence, azimuths, freq_min, freq_max, wavelet_freq
    )
    
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 2, 2],
                          height_ratios=[1, 1], hspace=0.35, wspace=0.35)
    
    # Velocity profile
    ax_vel = fig.add_subplot(gs[:, 0])
    velocity_profile_plot(layers, ax=ax_vel, title='Velocity Profile')
    
    # AVO Wiggle plot
    ax_avo_wig = fig.add_subplot(gs[0, 1])
    wiggle_plot(seismogram, time, angles, ax=ax_avo_wig,
                title='AVO Wiggle Display', scale=1.2)
    
    # AVO Image plot
    ax_avo_img = fig.add_subplot(gs[0, 2])
    image_plot(seismogram, time, angles, ax=ax_avo_img, title='AVO Image Display')
    
    # AVAZ Wiggle plot
    ax_avaz_wig = fig.add_subplot(gs[1, 1])
    wiggle_plot(seismogram_avaz, time_avaz, azimuths, ax=ax_avaz_wig,
                title=f'AVAZ Wiggle Display (θ={avaz_incidence}°)', scale=1.2)
    ax_avaz_wig.set_xlabel('Azimuth φ (°)', fontsize=12)
    
    # AVAZ Image plot
    ax_avaz_img = fig.add_subplot(gs[1, 2])
    image_plot(seismogram_avaz, time_avaz, azimuths, ax=ax_avaz_img,
               title=f'AVAZ Image Display (θ={avaz_incidence}°)')
    ax_avaz_img.set_xlabel('Azimuth φ (°)', fontsize=12)
    
    plt.savefig('docs/assets/screenshot_gather.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Generated screenshot_gather.png")


if __name__ == '__main__':
    print("Generating screenshots...")
    print()
    
    generate_screenshot_model()
    generate_screenshot_avo()
    generate_screenshot_gather()
    
    print()
    print("All screenshots generated successfully!")

