"""
seismic_plot.py - Professional Seismic Visualization Module

This module provides functions for professional seismic data visualization,
including:
  - Wiggle plots with variable area fill
  - Image plots with seismic colormaps
  - Velocity profile plots
  - AVO/AVA curve plots

Author: Rui Yang, 2024
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import PolyCollection
from typing import Optional, Tuple, List
import io


def seismic_colormap() -> mcolors.LinearSegmentedColormap:
    """
    Create a professional red-white-blue seismic colormap.

    Returns
    -------
    cmap : LinearSegmentedColormap
        Seismic colormap suitable for seismic amplitude display
    """
    colors = [
        (0.0, 'blue'),
        (0.25, 'lightblue'),
        (0.5, 'white'),
        (0.75, 'lightsalmon'),
        (1.0, 'red')
    ]
    positions = [c[0] for c in colors]
    color_values = [c[1] for c in colors]

    return mcolors.LinearSegmentedColormap.from_list('seismic_rwb',
                                                      list(zip(positions, color_values)))


def wiggle_plot(data: np.ndarray,
                time: np.ndarray,
                angles: np.ndarray,
                ax: Optional[plt.Axes] = None,
                fill_positive: bool = True,
                fill_negative: bool = False,
                positive_color: str = 'black',
                negative_color: str = 'red',
                line_color: str = 'black',
                line_width: float = 0.5,
                scale: float = 1.0,
                clip: float = 1.0,
                title: str = 'Seismic Gather',
                xlabel: str = 'Angle (deg)',
                ylabel: str = 'Time (s)') -> plt.Axes:
    """
    Create a professional wiggle plot with variable area fill.

    Parameters
    ----------
    data : ndarray
        2D seismic data array with shape (nTime, nAngles)
    time : ndarray
        1D time array (s)
    angles : ndarray
        1D angle array (degrees)
    ax : Axes, optional
        Matplotlib axes to plot on. If None, creates new figure.
    fill_positive : bool
        Fill positive amplitudes
    fill_negative : bool
        Fill negative amplitudes
    positive_color : str
        Color for positive fill
    negative_color : str
        Color for negative fill
    line_color : str
        Color for trace lines
    line_width : float
        Width of trace lines
    scale : float
        Amplitude scaling factor
    clip : float
        Clipping value (0 to 1)
    title : str
        Plot title
    xlabel : str
        X-axis label
    ylabel : str
        Y-axis label

    Returns
    -------
    ax : Axes
        Matplotlib axes with the plot
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))

    ntime, nangles = data.shape

    # Normalize data
    max_amp = np.max(np.abs(data)) + 1e-10
    data_norm = data / max_amp * scale

    # Clip data
    data_norm = np.clip(data_norm, -clip, clip)

    # Calculate trace spacing
    if len(angles) > 1:
        trace_spacing = np.mean(np.diff(angles))
    else:
        trace_spacing = 1.0

    # Plot each trace
    for i, angle in enumerate(angles):
        trace = data_norm[:, i] * trace_spacing * 0.8
        x_base = angle
        x_trace = x_base + trace

        # Plot trace line
        ax.plot(x_trace, time, color=line_color, linewidth=line_width)

        # Fill positive amplitudes
        if fill_positive:
            x_fill = np.copy(x_trace)
            x_fill[trace < 0] = x_base
            verts = list(zip(x_fill, time))
            verts.append((x_base, time[-1]))
            verts.append((x_base, time[0]))
            poly = plt.Polygon(verts, facecolor=positive_color,
                              edgecolor='none', alpha=0.8)
            ax.add_patch(poly)

        # Fill negative amplitudes
        if fill_negative:
            x_fill = np.copy(x_trace)
            x_fill[trace > 0] = x_base
            verts = list(zip(x_fill, time))
            verts.append((x_base, time[-1]))
            verts.append((x_base, time[0]))
            poly = plt.Polygon(verts, facecolor=negative_color,
                              edgecolor='none', alpha=0.8)
            ax.add_patch(poly)

    # Set axes properties
    ax.set_xlim(angles[0] - trace_spacing, angles[-1] + trace_spacing)
    ax.set_ylim(time[-1], time[0])  # Reverse y-axis (time increases downward)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    return ax


def image_plot(data: np.ndarray,
               time: np.ndarray,
               angles: np.ndarray,
               ax: Optional[plt.Axes] = None,
               cmap: str = 'seismic',
               vmin: Optional[float] = None,
               vmax: Optional[float] = None,
               aspect: str = 'auto',
               title: str = 'Seismic Gather',
               xlabel: str = 'Angle (deg)',
               ylabel: str = 'Time (s)',
               colorbar: bool = True,
               colorbar_label: str = 'Amplitude') -> Tuple[plt.Axes, Optional[plt.colorbar]]:
    """
    Create a seismic image plot with color mapping.

    Parameters
    ----------
    data : ndarray
        2D seismic data array with shape (nTime, nAngles)
    time : ndarray
        1D time array (s)
    angles : ndarray
        1D angle array (degrees)
    ax : Axes, optional
        Matplotlib axes to plot on
    cmap : str
        Colormap name ('seismic', 'RdBu_r', 'gray', etc.)
    vmin, vmax : float, optional
        Color scale limits
    aspect : str
        Aspect ratio ('auto', 'equal', or float)
    title : str
        Plot title
    xlabel, ylabel : str
        Axis labels
    colorbar : bool
        Whether to add colorbar
    colorbar_label : str
        Colorbar label

    Returns
    -------
    ax : Axes
        Matplotlib axes
    cbar : colorbar or None
        Colorbar object if colorbar=True
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))

    # Set symmetric color limits if not specified
    if vmin is None and vmax is None:
        vmax = np.max(np.abs(data))
        vmin = -vmax

    # Create image
    extent = [angles[0], angles[-1], time[-1], time[0]]
    im = ax.imshow(data, aspect=aspect, cmap=cmap,
                   extent=extent, vmin=vmin, vmax=vmax)

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')

    cbar = None
    if colorbar:
        cbar = plt.colorbar(im, ax=ax, label=colorbar_label, shrink=0.8)

    return ax, cbar


def velocity_profile_plot(layers: list,
                          ax: Optional[plt.Axes] = None,
                          show_vp: bool = True,
                          show_vs: bool = True,
                          show_rho: bool = True,
                          title: str = 'Velocity Profile') -> plt.Axes:
    """
    Plot velocity and density profiles for a layered model.

    Parameters
    ----------
    layers : list
        List of layer parameter dictionaries
    ax : Axes, optional
        Matplotlib axes to plot on
    show_vp, show_vs, show_rho : bool
        Whether to show each property
    title : str
        Plot title

    Returns
    -------
    ax : Axes
        Matplotlib axes with the plot
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 8))

    # Calculate depth
    depths = [0]
    vp_vals = []
    vs_vals = []
    rho_vals = []

    for layer in layers:
        d = depths[-1]
        h = layer['thickness']
        depths.extend([d, d + h])
        vp_vals.extend([layer['vp'], layer['vp']])
        vs_vals.extend([layer['vs'], layer['vs']])
        rho_vals.extend([layer['rho'], layer['rho']])

    depths = depths[1:]  # Remove initial 0

    # Plot
    if show_vp:
        ax.plot(vp_vals, depths, 'b-', linewidth=2, label='Vp (km/s)')
    if show_vs:
        ax.plot(vs_vals, depths, 'r-', linewidth=2, label='Vs (km/s)')
    if show_rho:
        ax.plot(rho_vals, depths, 'g-', linewidth=2, label=r'$\rho$ (g/cm³)')

    ax.set_ylim(max(depths), 0)  # Reverse y-axis
    ax.set_xlabel('Value', fontsize=12)
    ax.set_ylabel('Depth (km)', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    return ax


def avo_curve_plot(angles: np.ndarray,
                   Rpp: np.ndarray,
                   ax: Optional[plt.Axes] = None,
                   show_real: bool = True,
                   show_imag: bool = False,
                   show_abs: bool = True,
                   title: str = 'AVO Curve',
                   xlabel: str = 'Angle (deg)',
                   ylabel: str = 'Reflection Coefficient') -> plt.Axes:
    """
    Plot AVO (Amplitude vs Offset/Angle) curves.

    Parameters
    ----------
    angles : ndarray
        1D angle array (degrees)
    Rpp : ndarray
        Complex PP reflection coefficients
    ax : Axes, optional
        Matplotlib axes to plot on
    show_real, show_imag, show_abs : bool
        Whether to show each component
    title : str
        Plot title
    xlabel, ylabel : str
        Axis labels

    Returns
    -------
    ax : Axes
        Matplotlib axes with the plot
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    if show_real:
        ax.plot(angles, np.real(Rpp), 'b-', linewidth=2, label='Real(Rpp)')
    if show_imag:
        ax.plot(angles, np.imag(Rpp), 'r--', linewidth=2, label='Imag(Rpp)')
    if show_abs:
        ax.plot(angles, np.abs(Rpp), 'k-', linewidth=2, label='|Rpp|')

    ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    return ax


def create_seismic_figure(data: np.ndarray,
                          time: np.ndarray,
                          angles: np.ndarray,
                          layers: list,
                          Rpp: Optional[np.ndarray] = None,
                          display_mode: str = 'wiggle') -> plt.Figure:
    """
    Create a comprehensive seismic figure with multiple panels.

    Parameters
    ----------
    data : ndarray
        2D seismic data array with shape (nTime, nAngles)
    time : ndarray
        1D time array (s)
    angles : ndarray
        1D angle array (degrees)
    layers : list
        List of layer parameter dictionaries
    Rpp : ndarray, optional
        Complex reflection coefficients for AVO plot
    display_mode : str
        'wiggle', 'image', or 'both'

    Returns
    -------
    fig : Figure
        Matplotlib figure object
    """
    if display_mode == 'both':
        fig = plt.figure(figsize=(16, 10))
        gs = fig.add_gridspec(2, 3, width_ratios=[1, 2, 2],
                              height_ratios=[1, 1], hspace=0.3, wspace=0.3)

        # Velocity profile (left, spans both rows)
        ax_vel = fig.add_subplot(gs[:, 0])
        velocity_profile_plot(layers, ax=ax_vel, title='Model')

        # Wiggle plot (top right)
        ax_wiggle = fig.add_subplot(gs[0, 1])
        wiggle_plot(data, time, angles, ax=ax_wiggle,
                    title='Wiggle Display', scale=1.2)

        # Image plot (top far right)
        ax_image = fig.add_subplot(gs[0, 2])
        image_plot(data, time, angles, ax=ax_image, title='Image Display')

        # AVO curve (bottom, spans two columns)
        if Rpp is not None:
            ax_avo = fig.add_subplot(gs[1, 1:])
            avo_curve_plot(angles, Rpp, ax=ax_avo,
                          show_real=True, show_abs=True)

    elif display_mode == 'wiggle':
        fig = plt.figure(figsize=(14, 8))
        gs = fig.add_gridspec(1, 3, width_ratios=[1, 2, 1], wspace=0.3)

        ax_vel = fig.add_subplot(gs[0])
        velocity_profile_plot(layers, ax=ax_vel, title='Model')

        ax_wiggle = fig.add_subplot(gs[1])
        wiggle_plot(data, time, angles, ax=ax_wiggle, scale=1.2)

        if Rpp is not None:
            ax_avo = fig.add_subplot(gs[2])
            avo_curve_plot(angles, Rpp, ax=ax_avo)

    else:  # image mode
        fig = plt.figure(figsize=(14, 8))
        gs = fig.add_gridspec(1, 3, width_ratios=[1, 2, 1], wspace=0.3)

        ax_vel = fig.add_subplot(gs[0])
        velocity_profile_plot(layers, ax=ax_vel, title='Model')

        ax_image = fig.add_subplot(gs[1])
        image_plot(data, time, angles, ax=ax_image)

        if Rpp is not None:
            ax_avo = fig.add_subplot(gs[2])
            avo_curve_plot(angles, Rpp, ax=ax_avo)

    fig.tight_layout()
    return fig


def create_comprehensive_figure(seismogram_avo: np.ndarray,
                                time_avo: np.ndarray,
                                angles_avo: np.ndarray,
                                seismogram_avaz: np.ndarray,
                                time_avaz: np.ndarray,
                                azimuths_avaz: np.ndarray,
                                layers: list,
                                avaz_incidence: float = 30.0) -> plt.Figure:
    """
    Create a comprehensive figure with 5 subplots:
    1. Velocity profile
    2. AVO wiggle display
    3. AVO image display
    4. AVAZ wiggle display
    5. AVAZ image display

    Parameters
    ----------
    seismogram_avo : ndarray
        2D AVO seismic data array with shape (nTime, nAngles)
    time_avo : ndarray
        1D time array for AVO (s)
    angles_avo : ndarray
        1D angle array for AVO (degrees)
    seismogram_avaz : ndarray
        2D AVAZ seismic data array with shape (nTime, nAzimuths)
    time_avaz : ndarray
        1D time array for AVAZ (s)
    azimuths_avaz : ndarray
        1D azimuth array for AVAZ (degrees)
    layers : list
        List of layer parameter dictionaries
    avaz_incidence : float
        Fixed incidence angle for AVAZ (degrees)

    Returns
    -------
    fig : Figure
        Matplotlib figure object with 5 subplots
    """
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 2, 2],
                          height_ratios=[1, 1], hspace=0.35, wspace=0.35)

    # 1. Velocity profile (left, spans both rows)
    ax_vel = fig.add_subplot(gs[:, 0])
    velocity_profile_plot(layers, ax=ax_vel, title='Velocity Profile')

    # 2. AVO wiggle display (top middle)
    ax_avo_wig = fig.add_subplot(gs[0, 1])
    wiggle_plot(seismogram_avo, time_avo, angles_avo, ax=ax_avo_wig,
                title='AVO Wiggle Display', scale=1.2)

    # 3. AVO image display (top right)
    ax_avo_img = fig.add_subplot(gs[0, 2])
    image_plot(seismogram_avo, time_avo, angles_avo, ax=ax_avo_img,
               title='AVO Image Display')

    # 4. AVAZ wiggle display (bottom middle)
    ax_avaz_wig = fig.add_subplot(gs[1, 1])
    wiggle_plot(seismogram_avaz, time_avaz, azimuths_avaz, ax=ax_avaz_wig,
                title=f'AVAZ Wiggle Display (θ={avaz_incidence}°)', scale=1.2)
    ax_avaz_wig.set_xlabel('Azimuth φ (°)', fontsize=12)

    # 5. AVAZ image display (bottom right)
    ax_avaz_img = fig.add_subplot(gs[1, 2])
    image_plot(seismogram_avaz, time_avaz, azimuths_avaz, ax=ax_avaz_img,
               title=f'AVAZ Image Display (θ={avaz_incidence}°)')
    ax_avaz_img.set_xlabel('Azimuth φ (°)', fontsize=12)

    fig.tight_layout()
    return fig


def fig_to_bytes(fig: plt.Figure, format: str = 'png', dpi: int = 150) -> bytes:
    """
    Convert matplotlib figure to bytes for download.

    Parameters
    ----------
    fig : Figure
        Matplotlib figure
    format : str
        Output format ('png', 'pdf', 'svg')
    dpi : int
        Resolution for raster formats

    Returns
    -------
    bytes
        Figure data as bytes
    """
    buf = io.BytesIO()
    fig.savefig(buf, format=format, dpi=dpi, bbox_inches='tight')
    buf.seek(0)
    return buf.getvalue()
