"""
AzRM Seismic Forward Modeling Web Application

A Streamlit-based GUI for anisotropic reflectivity method seismic forward modeling.
Supports Isotropic, VTI, HTI, and Orthorhombic media.

Author: Rui Yang, 2024
"""

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Dict, Optional
import json

# Local modules
from stiffness import (
    build_model_stiffness,
    build_layer_stiffness,
    preset_two_layer_vti,
    preset_shale_sand,
    preset_fractured_reservoir
)
from seismic_plot import (
    wiggle_plot,
    image_plot,
    velocity_profile_plot,
    create_seismic_figure,
    create_comprehensive_figure,
    fig_to_bytes
)

# Try to import C++ module, fall back to pure Python if not available
try:
    import azrm
    CPP_AVAILABLE = True
except ImportError:
    CPP_AVAILABLE = False
    st.warning("C++ module not found. Please run `./build_python.sh` to compile.")


# =============================================================================
# Page configuration
# =============================================================================
st.set_page_config(
    page_title="AzRM Seismic Forward Modeling",
    page_icon="🌊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        font-size: 3.5rem;
        font-weight: bold;
        color: #1E88E5;
        text-align: center;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        font-size: 1.4rem;
        color: #666;
        text-align: center;
        margin-bottom: 0.5rem;
    }
    .author-info {
        font-size: 1.1rem;
        color: #888;
        text-align: center;
        margin-bottom: 2rem;
        font-style: italic;
    }
    .layer-card {
        background-color: #f0f2f6;
        border-radius: 10px;
        padding: 1rem;
        margin-bottom: 1rem;
    }
    .stButton>button {
        width: 100%;
    }
    /* Increase sidebar width */
    section[data-testid="stSidebar"] {
        min-width: 300px !important;
    }
    section[data-testid="stSidebar"] > div {
        width: 300px !important;
    }
</style>
""", unsafe_allow_html=True)


# =============================================================================
# Session state initialization
# =============================================================================
if 'layers' not in st.session_state:
    st.session_state.layers = preset_two_layer_vti()

if 'results' not in st.session_state:
    st.session_state.results = None


# =============================================================================
# Helper functions
# =============================================================================
def ricker_wavelet(freq: float, dt: float, duration: float = 0.1) -> np.ndarray:
    """Generate Ricker wavelet."""
    t = np.arange(-duration/2, duration/2, dt)
    w = (1 - 2 * (np.pi * freq * t)**2) * np.exp(-(np.pi * freq * t)**2)
    return w


def compute_single_interface_Rpp(layer_above: Dict, layer_below: Dict,
                                  angles: np.ndarray, freq: float, phi: float) -> np.ndarray:
    """
    Compute Rpp for a single interface between two layers.

    Returns
    -------
    Rpp : ndarray
        Complex reflection coefficients at each angle
    """
    # Build two-layer model
    layers = [layer_above, layer_below]
    thickness, density, stiffness = build_model_stiffness(layers)

    if CPP_AVAILABLE:
        Rpp, _, _ = azrm.compute_Rpp(freq, angles, thickness, density, stiffness, phi)
        # POLARITY FIX: Negate to correct 180° phase offset
        Rpp = -Rpp
    else:
        # Simple Zoeppritz approximation
        vp1, vs1, rho1 = layer_above['vp'], layer_above['vs'], layer_above['rho']
        vp2, vs2, rho2 = layer_below['vp'], layer_below['vs'], layer_below['rho']
        Z1, Z2 = rho1 * vp1, rho2 * vp2
        G1, G2 = rho1 * vs1**2, rho2 * vs2**2
        Rp0 = (Z2 - Z1) / (Z2 + Z1)
        Rp2 = 0.5 * (vp2/vp1 - 1) - 4*(vs1/vp1)**2 * (G2-G1)/(G1+G2)
        Rpp = np.array([Rp0 + Rp2 * np.sin(np.radians(a))**2 for a in angles])

    return Rpp


def compute_avaz(layer_above: Dict, layer_below: Dict,
                 angle: float, freq: float, azimuths: np.ndarray) -> np.ndarray:
    """
    Compute Rpp as a function of azimuth (AVAZ) for a single interface.

    Parameters
    ----------
    layer_above, layer_below : dict
        Layer parameters
    angle : float
        Fixed incidence angle (degrees)
    freq : float
        Frequency (Hz)
    azimuths : ndarray
        Array of azimuth angles (degrees)

    Returns
    -------
    Rpp : ndarray
        Complex reflection coefficients at each azimuth
    """
    layers = [layer_above, layer_below]
    thickness, density, stiffness = build_model_stiffness(layers)
    angles_arr = np.array([angle])

    Rpp_avaz = np.zeros(len(azimuths), dtype=complex)

    if CPP_AVAILABLE:
        for i, phi in enumerate(azimuths):
            Rpp, _, _ = azrm.compute_Rpp(freq, angles_arr, thickness, density, stiffness, phi)
            # POLARITY FIX: Negate to correct 180° phase offset
            Rpp_avaz[i] = -Rpp[0]
    else:
        # For isotropic/VTI, Rpp doesn't vary with azimuth
        Rpp_avaz[:] = compute_single_interface_Rpp(layer_above, layer_below, angles_arr, freq, 0)[0]

    return Rpp_avaz


def compute_avaz_seismogram(layers: List[Dict],
                             incidence_angle: float,
                             azimuths: np.ndarray,
                             freq_min: float,
                             freq_max: float,
                             wavelet_freq: float,
                             dt: float = 0.002) -> tuple:
    """
    Compute AVAZ seismogram (fixed incidence angle, varying azimuth).

    Parameters
    ----------
    layers : list of dict
        Layer parameters
    incidence_angle : float
        Fixed incidence angle (degrees)
    azimuths : ndarray
        Array of azimuth angles (degrees)
    freq_min, freq_max : float
        Frequency range (Hz)
    wavelet_freq : float
        Wavelet dominant frequency (Hz)
    dt : float
        Time sampling interval (s)

    Returns
    -------
    time : ndarray
        Time array (s)
    seismogram : ndarray
        Synthetic seismogram (nTime x nAzimuths)
    """
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

    # Collect Rpp for all azimuths
    Rpp_all = np.zeros((n_pos_freq, n_azimuths), dtype=complex)

    if CPP_AVAILABLE:
        for i, phi in enumerate(azimuths):
            Rpp_freq, _, _ = azrm.compute_Rpp_batch(
                frequencies, angles_arr, thickness, density, stiffness, phi
            )
            # POLARITY FIX
            Rpp_all[:, i] = -Rpp_freq[:, 0]
    else:
        # Fallback: isotropic approximation (no azimuth variation)
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


def compute_seismogram(layers: List[Dict],
                       angles: np.ndarray,
                       freq_min: float,
                       freq_max: float,
                       wavelet_freq: float,
                       phi: float,
                       dt: float = 0.002) -> tuple:
    """
    Compute synthetic seismogram using reflectivity method.

    The reflectivity method computes Rpp(f) which contains ALL internal multiples.
    The phase in Rpp(f) encodes the arrival times of primaries and multiples.

    Returns
    -------
    time : ndarray
        Time array (s)
    seismogram : ndarray
        Synthetic seismogram (nTime x nAngles)
    Rpp_mid : ndarray
        Reflection coefficients at mid-frequency for AVO display
    """
    # Build model
    thickness, density, stiffness = build_model_stiffness(layers)

    # Calculate total two-way time (add extra for multiples)
    vp_list = [layer['vp'] for layer in layers]
    twt_primary = 2 * sum(h / vp for h, vp in zip(thickness, vp_list))
    # Allow 3x primary time for multiples
    total_time = max(twt_primary * 3, 1.0)

    # Time and frequency arrays
    nsamples = int(total_time / dt) + 1
    nfft = 2 ** int(np.ceil(np.log2(nsamples)))
    df = 1.0 / (nfft * dt)
    nangles = len(angles)

    # Frequency array for positive frequencies: 0, df, 2*df, ..., (nfft/2)*df
    n_pos_freq = nfft // 2 + 1
    frequencies = np.arange(n_pos_freq) * df

    # Compute reflection coefficients using C++
    if CPP_AVAILABLE:
        # Compute at all positive frequencies
        Rpp_all, Rpsv_all, Rpsh_all = azrm.compute_Rpp_batch(
            frequencies, angles, thickness, density, stiffness, phi
        )
        # POLARITY FIX: The C++ code has a 180° phase offset due to eigenvector
        # normalization convention. Negate to correct the polarity.
        Rpp_all = -Rpp_all

        # IMPROVED FIX FOR θ=0° PHASE DISCONTINUITY:
        # The C++ eigenvector sign convention differs between normal (θ=0°) and
        # oblique incidence, causing a ~180° phase jump. Detect and correct this.
        if len(angles) > 1 and angles[0] == 0:
            # Check phase at multiple frequencies for more robust detection
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
        # Fall back to simple approximation (for testing UI without C++)
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
                    # Add phase for travel time to first interface
                    t_delay = 2 * thickness[0] / vp1
                    phase = np.exp(-1j * 2 * np.pi * f * t_delay)
                    Rpp_all[i, j] = (Rp0 + Rp2 * np.sin(theta)**2) * phase

    # Apply Ricker wavelet in frequency domain with bandpass
    wavelet_spectrum = np.zeros(n_pos_freq)
    for i, f in enumerate(frequencies):
        if f > 0:
            # Ricker wavelet spectrum
            u = f / wavelet_freq
            wavelet_spectrum[i] = (2 / np.sqrt(np.pi)) * u**2 * np.exp(-u**2)

            # Apply bandpass (taper outside freq_min to freq_max)
            if f < freq_min:
                taper = 0.5 * (1 - np.cos(np.pi * f / max(freq_min, 0.1)))
                wavelet_spectrum[i] *= taper
            elif f > freq_max:
                wavelet_spectrum[i] = 0

    # Multiply by wavelet
    Rpp_wavelet = Rpp_all * wavelet_spectrum[:, np.newaxis]

    # Build full FFT spectrum (Hermitian symmetric for real output)
    full_spectrum = np.zeros((nfft, nangles), dtype=complex)

    # Positive frequencies (indices 0 to nfft/2)
    full_spectrum[:n_pos_freq, :] = Rpp_wavelet

    # Negative frequencies (indices nfft/2+1 to nfft-1)
    # f[-k] = conj(f[k]) for k = 1, 2, ..., nfft/2 - 1
    for k in range(1, nfft // 2):
        full_spectrum[nfft - k, :] = np.conj(Rpp_wavelet[k, :])

    # Inverse FFT to get time domain
    seismogram = np.fft.ifft(full_spectrum, axis=0)
    seismogram = np.real(seismogram[:nsamples, :]) * nfft * df

    # Normalize
    max_amp = np.max(np.abs(seismogram))
    if max_amp > 0:
        seismogram = seismogram / max_amp

    # Time array
    time = np.arange(nsamples) * dt

    # Get Rpp at wavelet dominant frequency for AVO display
    avo_freq_idx = np.argmin(np.abs(frequencies - wavelet_freq))
    Rpp_mid = Rpp_all[avo_freq_idx, :]

    return time, seismogram, Rpp_mid


# =============================================================================
# Sidebar - Model Parameters
# =============================================================================
with st.sidebar:
    st.markdown("## ⚙️ Model Parameters")

    # Preset models
    st.markdown("### 📋 Preset Models")
    preset = st.selectbox(
        "Load preset model:",
        ["Custom", "Two-layer VTI", "Shale-Sand", "Fractured Reservoir"],
        key="preset_select"
    )

    if preset == "Two-layer VTI":
        st.session_state.layers = preset_two_layer_vti()
    elif preset == "Shale-Sand":
        st.session_state.layers = preset_shale_sand()
    elif preset == "Fractured Reservoir":
        st.session_state.layers = preset_fractured_reservoir()

    st.markdown("---")

    # Layer controls
    st.markdown("### 📊 Layer Configuration")

    col1, col2 = st.columns(2)
    with col1:
        if st.button("➕ Add Layer"):
            new_layer = {
                'type': 'VTI',
                'thickness': 0.3,
                'vp': 3.0,
                'vs': 1.5,
                'rho': 2.4,
                'epsilon': 0.05,
                'delta': 0.02,
                'gamma': 0.03,
                'crack_density': 0.0
            }
            st.session_state.layers.append(new_layer)
            st.rerun()

    with col2:
        if len(st.session_state.layers) > 1:
            if st.button("➖ Remove"):
                st.session_state.layers.pop()
                st.rerun()

    st.markdown(f"**Number of layers: {len(st.session_state.layers)}**")

    st.markdown("---")

    # Forward modeling parameters
    st.markdown("### 🔧 Forward Modeling Parameters")

    col1, col2 = st.columns(2)
    with col1:
        angle_min = st.number_input("θ min (°)", value=0, min_value=0, max_value=60)
        freq_min = st.number_input("f min (Hz)", value=0.0, min_value=0.0)
    with col2:
        angle_max = st.number_input("θ max (°)", value=40, min_value=1, max_value=60)
        freq_max = st.number_input("f max (Hz)", value=60.0, min_value=1.0)

    angle_step = st.slider("Angle step (°)", 0.5, 5.0, 1.0, 0.5)
    wavelet_freq = st.slider("Wavelet frequency (Hz)", 5, 60, 25)
    phi = st.slider("Azimuth φ (°)", 0, 180, 0)

    st.markdown("---")

    # Display options
    st.markdown("### 🎨 Display Options")
    display_mode = st.radio(
        "Display mode:",
        ["Wiggle", "Image", "Both"],
        horizontal=True
    )

    st.markdown("---")

    # Run button
    run_button = st.button("🚀 Run Forward Modeling", type="primary", use_container_width=True)


# =============================================================================
# Main area - Layer editing
# =============================================================================
st.markdown('<p class="main-header">🌊 AzRM Seismic Forward Modeling</p>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Anisotropic Reflectivity Method for Layered Media</p>', unsafe_allow_html=True)
st.markdown('<p class="author-info">Rui Yang, Tongji University/ Stanford University, 2024-2025</p>', unsafe_allow_html=True)

# Create tabs
tab1, tab2, tab3 = st.tabs(["📝 Layer Editor", "📊 Results", "ℹ️ About"])

with tab1:
    st.markdown("### Edit Layer Properties")

    # Create columns for layer cards
    for i, layer in enumerate(st.session_state.layers):
        with st.expander(f"**Layer {i+1}** - {layer['type']}", expanded=(i < 3)):
            col1, col2, col3 = st.columns(3)

            with col1:
                layer['type'] = st.selectbox(
                    "Media Type",
                    ["ISO", "VTI", "HTI", "OA"],
                    index=["ISO", "VTI", "HTI", "OA"].index(layer.get('type', 'VTI')),
                    key=f"type_{i}"
                )
                layer['thickness'] = st.number_input(
                    "Thickness (km)",
                    value=float(layer.get('thickness', 0.3)),
                    min_value=0.01,
                    step=0.1,
                    key=f"h_{i}"
                )

            with col2:
                layer['vp'] = st.number_input(
                    "Vp (km/s)",
                    value=float(layer.get('vp', 3.0)),
                    min_value=0.5,
                    max_value=8.0,
                    step=0.1,
                    key=f"vp_{i}"
                )
                layer['vs'] = st.number_input(
                    "Vs (km/s)",
                    value=float(layer.get('vs', 1.5)),
                    min_value=0.1,
                    max_value=5.0,
                    step=0.1,
                    key=f"vs_{i}"
                )
                layer['rho'] = st.number_input(
                    "Density (g/cm³)",
                    value=float(layer.get('rho', 2.4)),
                    min_value=1.0,
                    max_value=4.0,
                    step=0.1,
                    key=f"rho_{i}"
                )

            with col3:
                # Show Thomsen parameters for VTI and OA
                if layer['type'] in ['VTI', 'OA']:
                    layer['epsilon'] = st.number_input(
                        "ε (epsilon)",
                        value=float(layer.get('epsilon', 0.05)),
                        min_value=-0.3,
                        max_value=0.5,
                        step=0.01,
                        format="%.3f",
                        key=f"eps_{i}"
                    )
                    layer['delta'] = st.number_input(
                        "δ (delta)",
                        value=float(layer.get('delta', 0.02)),
                        min_value=-0.3,
                        max_value=0.3,
                        step=0.01,
                        format="%.3f",
                        key=f"del_{i}"
                    )
                    layer['gamma'] = st.number_input(
                        "γ (gamma)",
                        value=float(layer.get('gamma', 0.03)),
                        min_value=-0.3,
                        max_value=0.5,
                        step=0.01,
                        format="%.3f",
                        key=f"gam_{i}"
                    )

                # Show crack density for HTI and OA
                if layer['type'] in ['HTI', 'OA']:
                    layer['crack_density'] = st.number_input(
                        "Crack density",
                        value=float(layer.get('crack_density', 0.05)),
                        min_value=0.0,
                        max_value=0.15,
                        step=0.01,
                        format="%.3f",
                        key=f"cd_{i}"
                    )

    # Model summary table
    st.markdown("### Model Summary")
    summary_data = []
    total_thickness = 0
    for i, layer in enumerate(st.session_state.layers):
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
            row['e (crack)'] = layer.get('crack_density', 0)
        summary_data.append(row)
        total_thickness += layer['thickness']

    df = pd.DataFrame(summary_data)
    st.dataframe(df, use_container_width=True, hide_index=True)
    st.info(f"Total model thickness: **{total_thickness:.2f} km**")

with tab2:
    # Run forward modeling if button clicked
    if run_button:
        with st.spinner("Computing synthetic seismogram..."):
            angles = np.arange(angle_min, angle_max + angle_step, angle_step)

            try:
                time, seismogram, Rpp = compute_seismogram(
                    st.session_state.layers,
                    angles,
                    freq_min,
                    freq_max,
                    wavelet_freq,
                    phi
                )

                st.session_state.results = {
                    'time': time,
                    'angles': angles,
                    'seismogram': seismogram,
                    'Rpp': Rpp,
                    'phi': phi,
                    'wavelet_freq': wavelet_freq
                }

                st.success("✅ Computation completed!")

            except Exception as e:
                st.error(f"❌ Error: {str(e)}")
                st.session_state.results = None

    # Display results
    if st.session_state.results is not None:
        results = st.session_state.results
        time = results['time']
        angles = results['angles']
        seismogram = results['seismogram'].copy()
        Rpp = results['Rpp']

        # Create visualization
        st.markdown("### Synthetic Seismic Gather")

        if display_mode == "Both":
            col1, col2, col3 = st.columns([1, 2, 2])

            with col1:
                st.markdown("#### Velocity Profile")
                fig_vel, ax_vel = plt.subplots(figsize=(4, 6))
                velocity_profile_plot(st.session_state.layers, ax=ax_vel)
                st.pyplot(fig_vel)
                plt.close(fig_vel)

            with col2:
                st.markdown("#### Wiggle Display")
                fig_wig, ax_wig = plt.subplots(figsize=(8, 6))
                wiggle_plot(seismogram, time, angles, ax=ax_wig, scale=1.5)
                st.pyplot(fig_wig)
                plt.close(fig_wig)

            with col3:
                st.markdown("#### Image Display")
                fig_img, ax_img = plt.subplots(figsize=(8, 6))
                image_plot(seismogram, time, angles, ax=ax_img)
                st.pyplot(fig_img)
                plt.close(fig_img)

        elif display_mode == "Wiggle":
            col1, col2 = st.columns([1, 3])

            with col1:
                st.markdown("#### Velocity Profile")
                fig_vel, ax_vel = plt.subplots(figsize=(4, 6))
                velocity_profile_plot(st.session_state.layers, ax=ax_vel)
                st.pyplot(fig_vel)
                plt.close(fig_vel)

            with col2:
                st.markdown("#### Wiggle Display")
                fig_wig, ax_wig = plt.subplots(figsize=(12, 8))
                wiggle_plot(seismogram, time, angles, ax=ax_wig, scale=1.5)
                st.pyplot(fig_wig)
                plt.close(fig_wig)

        else:  # Image mode
            col1, col2 = st.columns([1, 3])

            with col1:
                st.markdown("#### Velocity Profile")
                fig_vel, ax_vel = plt.subplots(figsize=(4, 6))
                velocity_profile_plot(st.session_state.layers, ax=ax_vel)
                st.pyplot(fig_vel)
                plt.close(fig_vel)

            with col2:
                st.markdown("#### Image Display")
                fig_img, ax_img = plt.subplots(figsize=(12, 8))
                image_plot(seismogram, time, angles, ax=ax_img)
                st.pyplot(fig_img)
                plt.close(fig_img)

        # =================================================================
        # AVAZ Seismic Gather Section
        # =================================================================
        st.markdown("---")
        st.markdown("### AVAZ Seismic Gather")
        st.markdown("*Fixed incidence angle, varying azimuth*")

        col_avaz_ctrl1, col_avaz_ctrl2 = st.columns(2)
        with col_avaz_ctrl1:
            avaz_incidence = st.slider(
                "Fixed incidence angle θ (°)",
                min_value=0, max_value=60, value=30,
                key="avaz_incidence_angle"
            )
        with col_avaz_ctrl2:
            avaz_step = st.slider(
                "Azimuth step (°)",
                min_value=1, max_value=10, value=5,
                key="avaz_step"
            )

        # Compute AVAZ seismogram
        azimuths_avaz = np.arange(0, 181, avaz_step)
        with st.spinner("Computing AVAZ seismogram..."):
            time_avaz, seismogram_avaz = compute_avaz_seismogram(
                st.session_state.layers,
                avaz_incidence,
                azimuths_avaz,
                freq_min,
                freq_max,
                results['wavelet_freq']
            )

        # Display AVAZ gather
        col_avaz_wig, col_avaz_img = st.columns(2)

        with col_avaz_wig:
            st.markdown(f"#### AVAZ Wiggle Display (θ={avaz_incidence}°)")
            fig_avaz_wig, ax_avaz_wig = plt.subplots(figsize=(8, 6))
            wiggle_plot(seismogram_avaz, time_avaz, azimuths_avaz, ax=ax_avaz_wig, scale=1.5)
            ax_avaz_wig.set_xlabel('Azimuth φ (°)', fontsize=12)
            ax_avaz_wig.set_title(f'AVAZ Gather (θ={avaz_incidence}°)', fontsize=12)
            st.pyplot(fig_avaz_wig)
            plt.close(fig_avaz_wig)

        with col_avaz_img:
            st.markdown(f"#### AVAZ Image Display (θ={avaz_incidence}°)")
            fig_avaz_img, ax_avaz_img = plt.subplots(figsize=(8, 6))
            image_plot(seismogram_avaz, time_avaz, azimuths_avaz, ax=ax_avaz_img)
            ax_avaz_img.set_xlabel('Azimuth φ (°)', fontsize=12)
            ax_avaz_img.set_title(f'AVAZ Gather (θ={avaz_incidence}°)', fontsize=12)
            st.pyplot(fig_avaz_img)
            plt.close(fig_avaz_img)

        # Check for azimuthal variation
        max_trace_amp = np.max(np.abs(seismogram_avaz), axis=0)
        if len(max_trace_amp) > 1:
            amp_variation = (np.max(max_trace_amp) - np.min(max_trace_amp)) / (np.mean(max_trace_amp) + 1e-10)
            if amp_variation > 0.1:
                st.warning(f"⚠️ Significant azimuthal amplitude variation detected: {amp_variation*100:.1f}%")
            else:
                st.success(f"✓ Low azimuthal amplitude variation: {amp_variation*100:.1f}%")

        # =================================================================
        # AVO and AVAZ Analysis Section
        # =================================================================
        st.markdown("---")
        st.markdown("### AVO & AVAZ Analysis")

        # Interface selection
        n_layers = len(st.session_state.layers)
        if n_layers >= 2:
            interface_options = [f"Interface {i+1}: Layer {i+1} / Layer {i+2}"
                                for i in range(n_layers - 1)]
            selected_interface = st.selectbox(
                "Select interface for AVO/AVAZ analysis:",
                interface_options,
                key="interface_select"
            )
            interface_idx = interface_options.index(selected_interface)

            layer_above = st.session_state.layers[interface_idx]
            layer_below = st.session_state.layers[interface_idx + 1]

            # Display interface info
            st.info(f"**Layer {interface_idx+1}** ({layer_above['type']}): "
                   f"Vp={layer_above['vp']:.2f}, Vs={layer_above['vs']:.2f}, ρ={layer_above['rho']:.2f}  →  "
                   f"**Layer {interface_idx+2}** ({layer_below['type']}): "
                   f"Vp={layer_below['vp']:.2f}, Vs={layer_below['vs']:.2f}, ρ={layer_below['rho']:.2f}")

            # Two columns for AVO and AVAZ
            col_avo, col_avaz = st.columns(2)

            with col_avo:
                st.markdown("#### AVO Curve (|Rpp| vs θ)")

                # AVO parameters - azimuth angle slider
                avo_phi = st.slider("Azimuth angle φ (°)",
                                   min_value=0, max_value=180, value=int(results.get('phi', 0)),
                                   key="avo_phi_slider")

                # Compute AVO for selected interface
                avo_freq = results.get('wavelet_freq', 25)
                Rpp_avo = compute_single_interface_Rpp(
                    layer_above, layer_below, angles, avo_freq, avo_phi
                )

                fig_avo, ax_avo = plt.subplots(figsize=(6, 4))
                ax_avo.plot(angles, np.abs(Rpp_avo), 'b-', linewidth=2, label='|Rpp|')
                ax_avo.set_xlabel('Incidence Angle (°)', fontsize=12)
                ax_avo.set_ylabel('|Rpp|', fontsize=12)
                ax_avo.set_title(f'AVO at Interface {interface_idx+1} (φ={avo_phi}°)', fontsize=12)
                ax_avo.grid(True, alpha=0.3)
                ax_avo.set_xlim(angles[0], angles[-1])
                ax_avo.legend()
                st.pyplot(fig_avo)
                plt.close(fig_avo)

            with col_avaz:
                st.markdown("#### AVAZ Curve (|Rpp| vs φ)")

                # AVAZ parameters
                avaz_angle = st.slider("Fixed incidence angle for AVAZ (°)",
                                       int(angles[0]), int(angles[-1]), 20,
                                       key="avaz_angle")

                # Compute AVAZ for selected interface
                azimuths = np.arange(0, 181, 5)
                Rpp_avaz = compute_avaz(
                    layer_above, layer_below, avaz_angle, avo_freq, azimuths
                )

                fig_avaz, ax_avaz = plt.subplots(figsize=(6, 4))
                ax_avaz.plot(azimuths, np.abs(Rpp_avaz), 'r-', linewidth=2, label='|Rpp|')
                ax_avaz.set_xlabel('Azimuth φ (°)', fontsize=12)
                ax_avaz.set_ylabel('|Rpp|', fontsize=12)
                ax_avaz.set_title(f'AVAZ at Interface {interface_idx+1} (θ={avaz_angle}°)', fontsize=12)
                ax_avaz.grid(True, alpha=0.3)
                ax_avaz.set_xlim(0, 180)
                ax_avaz.legend()
                st.pyplot(fig_avaz)
                plt.close(fig_avaz)

            # Show anisotropy indicator
            if len(Rpp_avaz) > 1:
                avaz_variation = (np.max(np.abs(Rpp_avaz)) - np.min(np.abs(Rpp_avaz))) / (np.mean(np.abs(Rpp_avaz)) + 1e-10)
                if avaz_variation > 0.05:
                    st.warning(f"⚠️ Significant azimuthal anisotropy detected! "
                              f"AVAZ variation: {avaz_variation*100:.1f}%")
                else:
                    st.success(f"✓ Low azimuthal anisotropy. AVAZ variation: {avaz_variation*100:.1f}%")

        else:
            st.warning("Need at least 2 layers for AVO/AVAZ analysis.")

        # Export options
        st.markdown("---")
        st.markdown("### 💾 Export")
        col1, col2, col3 = st.columns(3)

        with col1:
            # Create comprehensive figure with 5 subplots for export
            # Get AVAZ parameters from current sliders
            avaz_incidence_export = st.session_state.get('avaz_incidence_angle', 30)
            avaz_step_export = st.session_state.get('avaz_step', 5)
            azimuths_avaz_export = np.arange(0, 181, avaz_step_export)
            
            # Compute AVAZ seismogram for export
            time_avaz_export, seismogram_avaz_export = compute_avaz_seismogram(
                st.session_state.layers,
                avaz_incidence_export,
                azimuths_avaz_export,
                freq_min,
                freq_max,
                results['wavelet_freq']
            )
            
            # Create comprehensive figure with 5 subplots
            fig_export = create_comprehensive_figure(
                seismogram, time, angles,
                seismogram_avaz_export, time_avaz_export, azimuths_avaz_export,
                st.session_state.layers,
                avaz_incidence_export
            )
            png_bytes = fig_to_bytes(fig_export, 'png', dpi=300)
            st.download_button(
                "📥 Download PNG",
                data=png_bytes,
                file_name="seismic_gather.png",
                mime="image/png"
            )
            plt.close(fig_export)

        with col2:
            # Export model as JSON
            model_json = json.dumps(st.session_state.layers, indent=2)
            st.download_button(
                "📥 Download Model (JSON)",
                data=model_json,
                file_name="model.json",
                mime="application/json"
            )

        with col3:
            # Export data as MAT format
            from scipy.io import savemat
            import io
            
            # Get AVAZ data for export
            avaz_incidence_export = st.session_state.get('avaz_incidence_angle', 30)
            avaz_step_export = st.session_state.get('avaz_step', 5)
            azimuths_avaz_export = np.arange(0, 181, avaz_step_export)
            
            # Compute AVAZ seismogram for export
            time_avaz_export, seismogram_avaz_export = compute_avaz_seismogram(
                st.session_state.layers,
                avaz_incidence_export,
                azimuths_avaz_export,
                freq_min,
                freq_max,
                results['wavelet_freq']
            )
            
            # Prepare data dictionary for MAT file
            # Convert layers to JSON string for better MATLAB compatibility
            layers_json_str = json.dumps(st.session_state.layers)
            
            mat_data = {
                'time_avo': time,
                'angles': angles,
                'seismogram_avo': seismogram,
                'Rpp_real': np.real(Rpp),
                'Rpp_imag': np.imag(Rpp),
                'time_avaz': time_avaz_export,
                'azimuths': azimuths_avaz_export,
                'seismogram_avaz': seismogram_avaz_export,
                'avaz_incidence_angle': np.array([avaz_incidence_export]),
                'layers_json': layers_json_str,
                'freq_min': np.array([freq_min]),
                'freq_max': np.array([freq_max]),
                'wavelet_freq': np.array([results['wavelet_freq']]),
                'phi': np.array([results['phi']])
            }
            
            # Save to BytesIO buffer
            mat_buffer = io.BytesIO()
            savemat(mat_buffer, mat_data, format='5', oned_as='column')
            mat_buffer.seek(0)

            st.download_button(
                "📥 Download Data (MAT)",
                data=mat_buffer.getvalue(),
                file_name="seismic_data.mat",
                mime="application/octet-stream"
            )

    else:
        st.info("👆 Configure the model in the sidebar and click **Run Forward Modeling** to see results.")

with tab3:
    st.markdown("""
    ## About AzRM

    **AzRM (Azimuthal Reflectivity Method)** is a seismic forward modeling tool
    that computes synthetic seismograms for anisotropic layered media using
    the reflectivity method (Fryer & Frazer, 1984, 1987).

    ### Supported Media Types

    | Type | Description | Parameters |
    |------|-------------|------------|
    | **ISO** | Isotropic | Vp, Vs, ρ |
    | **VTI** | Vertical Transverse Isotropy | Vp, Vs, ρ, ε, δ, γ |
    | **HTI** | Horizontal Transverse Isotropy | Vp, Vs, ρ, crack density |
    | **OA** | Orthorhombic | Vp, Vs, ρ, ε, δ, γ, crack density |

    ### Algorithm

    1. Build 6×6 stiffness tensors for each layer
    2. Solve Christoffel equation using numerical eigenvalue decomposition
    3. Apply recursive reflectivity algorithm (Kennett, 1983)
    4. Compute frequency-domain reflection coefficients
    5. Apply Ricker wavelet and inverse FFT

    ### References
    - Yang R, Chen H, Guo Z, et al. An effective azimuthal reflectivity modeling (AzRM) tool for generating seismic data in anisotropic shale reservoirs[J]. Geophysics, 2025, 90(5): 1-76.
    - Fryer G J, Frazer L N. Seismic waves in stratified anisotropic media[J]. Geophysical Journal International, 1984, 78(3): 691-710.
    - Fryer G J, Frazer L N. Seismic waves in stratified anisotropic media—II. Elastodynamic eigensolutions for some anisotropic systems[J]. Geophysical Journal International, 1987, 91(1): 73-101.
    - Schoenberg M, Helbig K. Orthorhombic media: Modeling elastic wave behavior in a vertically fractured earth[J]. Geophysics, 1997, 62(6): 1954-1974.
    - Kennett B L. Seismic wave propagation in stratified anisotropic media[J]. Geophysical Journal International, 1983, 72(1): 1-37.

    ### Author

    Rui Yang, Tongji University, 2025

    ### Technical Info
    """)

    if CPP_AVAILABLE:
        st.success(f"✅ C++ module loaded (version {azrm.__version__})")
        st.info(f"OpenMP threads: {azrm.get_num_threads()}")
    else:
        st.warning("⚠️ C++ module not available. Using Python fallback (slower).")
        st.markdown("""
        To compile the C++ module:
        ```bash
        cd /path/to/project
        chmod +x build_python.sh
        ./build_python.sh
        ```
        """)
