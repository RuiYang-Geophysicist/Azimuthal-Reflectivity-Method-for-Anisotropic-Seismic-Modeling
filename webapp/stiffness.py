"""
stiffness.py - Stiffness Matrix Construction for Anisotropic Media

This module provides functions to construct 6x6 stiffness tensors (Voigt notation)
for various types of anisotropic media:
  - Isotropic (ISO)
  - VTI (Vertical Transverse Isotropy)
  - HTI (Horizontal Transverse Isotropy)
  - Orthorhombic (OA)

Author: Rui Yang, 2024
"""

import numpy as np
from typing import Tuple, Optional


def stiffness_isotropic(vp: float, vs: float, rho: float) -> np.ndarray:
    """
    Construct stiffness matrix for isotropic medium.

    Parameters
    ----------
    vp : float
        P-wave velocity (km/s)
    vs : float
        S-wave velocity (km/s)
    rho : float
        Density (g/cm^3)

    Returns
    -------
    C : ndarray
        6x6 stiffness matrix (GPa)
    """
    # Lame parameters (in GPa, since rho in g/cm^3 and v in km/s)
    mu = rho * vs**2
    lam = rho * vp**2 - 2 * mu

    C = np.zeros((6, 6))

    # C11 = C22 = C33 = lambda + 2*mu
    C[0, 0] = C[1, 1] = C[2, 2] = lam + 2 * mu

    # C12 = C13 = C23 = lambda
    C[0, 1] = C[1, 0] = lam
    C[0, 2] = C[2, 0] = lam
    C[1, 2] = C[2, 1] = lam

    # C44 = C55 = C66 = mu
    C[3, 3] = C[4, 4] = C[5, 5] = mu

    return C


def stiffness_vti(vp: float, vs: float, rho: float,
                   epsilon: float, delta: float, gamma: float) -> np.ndarray:
    """
    Construct stiffness matrix for VTI medium using Thomsen parameters.

    Parameters
    ----------
    vp : float
        Vertical P-wave velocity (km/s)
    vs : float
        Vertical S-wave velocity (km/s)
    rho : float
        Density (g/cm^3)
    epsilon : float
        Thomsen parameter epsilon (P-wave anisotropy)
    delta : float
        Thomsen parameter delta (near-vertical P-wave anisotropy)
    gamma : float
        Thomsen parameter gamma (SH-wave anisotropy)

    Returns
    -------
    C : ndarray
        6x6 stiffness matrix (GPa)

    Reference
    ---------
    Thomsen, L. (1986). Weak elastic anisotropy. Geophysics, 51(10), 1954-1966.
    """
    # Vertical stiffnesses
    C33 = rho * vp**2
    C44 = rho * vs**2
    C55 = C44  # VTI symmetry

    # Horizontal stiffnesses from Thomsen parameters
    C11 = C33 * (1 + 2 * epsilon)
    C66 = C44 * (1 + 2 * gamma)

    # C13 from delta (Thomsen's formula)
    # delta = (C13 + C44)^2 - (C33 - C44)^2 / (2*C33*(C33-C44))
    # Solving for C13:
    C13 = np.sqrt(2 * delta * C33 * (C33 - C55) + (C33 - C55)**2) - C55

    # C12 from VTI symmetry: C12 = C11 - 2*C66
    C12 = C11 - 2 * C66

    # Build 6x6 matrix
    C = np.zeros((6, 6))
    C[0, 0] = C11
    C[1, 1] = C11  # VTI: C22 = C11
    C[2, 2] = C33
    C[0, 1] = C[1, 0] = C12
    C[0, 2] = C[2, 0] = C13
    C[1, 2] = C[2, 1] = C13  # VTI: C23 = C13
    C[3, 3] = C44
    C[4, 4] = C44  # VTI: C55 = C44
    C[5, 5] = C66

    return C


def stiffness_hti_hudson(vp: float, vs: float, rho: float,
                          crack_density: float, aspect_ratio: float = 0.01) -> np.ndarray:
    """
    Construct stiffness matrix for HTI medium using Hudson's crack theory.

    The HTI medium is modeled as an isotropic background with aligned
    vertical penny-shaped cracks. The crack normal is along the x1 axis.

    Parameters
    ----------
    vp : float
        Background P-wave velocity (km/s)
    vs : float
        Background S-wave velocity (km/s)
    rho : float
        Density (g/cm^3)
    crack_density : float
        Crack density parameter e = N * a^3 / V, typically 0 to 0.1
    aspect_ratio : float, optional
        Crack aspect ratio (thickness/diameter), default 0.01

    Returns
    -------
    C : ndarray
        6x6 stiffness matrix (GPa)

    Reference
    ---------
    Hudson, J.A. (1981). Wave speeds and attenuation of elastic waves in
    material containing cracks. Geophys. J. R. Astr. Soc., 64, 133-150.
    """
    # Background isotropic moduli
    mu = rho * vs**2
    lam = rho * vp**2 - 2 * mu
    K = lam + 2 * mu / 3  # Bulk modulus

    # Poisson's ratio
    nu = (vp**2 - 2 * vs**2) / (2 * (vp**2 - vs**2))

    # Hudson's U parameters (for dry cracks)
    # For penny-shaped cracks with small aspect ratio
    U1 = 16 * (1 - nu) / (3 * (2 - nu))
    U3 = 4 * (1 - nu) / (3 * (2 - nu) / (1 - nu / 2))

    # First-order corrections
    e = crack_density
    dC11 = -lam**2 * e * U1 / (mu * (lam + 2 * mu))
    dC13 = -lam * e * U1 / (lam + 2 * mu)
    dC33 = -(lam + 2 * mu) * e * U1
    dC44 = 0  # For vertical cracks with x1 normal
    dC55 = -mu * e * U3
    dC66 = -mu * e * U3

    # Build HTI stiffness matrix
    # HTI with x1 as symmetry axis (crack normal along x1)
    C = np.zeros((6, 6))

    # Background isotropic
    C[0, 0] = lam + 2 * mu + dC33  # C11 in HTI = C33 in VTI
    C[1, 1] = lam + 2 * mu + dC11  # C22
    C[2, 2] = lam + 2 * mu + dC11  # C33
    C[0, 1] = C[1, 0] = lam + dC13
    C[0, 2] = C[2, 0] = lam + dC13
    C[1, 2] = C[2, 1] = lam + dC11 - 2 * (mu + dC66)
    C[3, 3] = mu + dC66  # C44
    C[4, 4] = mu + dC55  # C55
    C[5, 5] = mu + dC55  # C66

    return C


def stiffness_orthorhombic(vp: float, vs: float, rho: float,
                            epsilon_v: float, delta_v: float, gamma_v: float,
                            crack_density: float) -> np.ndarray:
    """
    Construct stiffness matrix for orthorhombic medium.

    The orthorhombic medium is modeled as a VTI background with
    a single set of vertical aligned cracks (HTI perturbation).

    Parameters
    ----------
    vp : float
        Vertical P-wave velocity (km/s)
    vs : float
        Vertical S-wave velocity (km/s)
    rho : float
        Density (g/cm^3)
    epsilon_v : float
        VTI Thomsen parameter epsilon
    delta_v : float
        VTI Thomsen parameter delta
    gamma_v : float
        VTI Thomsen parameter gamma
    crack_density : float
        Crack density for vertical cracks (0 to 0.1)

    Returns
    -------
    C : ndarray
        6x6 stiffness matrix (GPa)

    Reference
    ---------
    Tsvankin, I. (1997). Anisotropic parameters and P-wave velocity for
    orthorhombic media. Geophysics, 62(4), 1292-1309.
    """
    # Start with VTI background
    C_vti = stiffness_vti(vp, vs, rho, epsilon_v, delta_v, gamma_v)

    # Add crack perturbation using Hudson theory
    # Crack normal along x1 direction

    mu = rho * vs**2
    lam = rho * vp**2 - 2 * mu

    # Poisson's ratio
    nu = (vp**2 - 2 * vs**2) / (2 * (vp**2 - vs**2) + 1e-10)

    # Hudson U parameters
    U1 = 16 * (1 - nu) / (3 * (2 - nu) + 1e-10)
    U3 = 4 * (1 - nu) / (3 * (2 - nu) / (1 - nu / 2 + 1e-10) + 1e-10)

    e = crack_density

    # Crack perturbations (simplified)
    dC11 = -(lam + 2 * mu) * e * U1
    dC12 = -lam * e * U1
    dC13 = -lam * e * U1
    dC55 = -mu * e * U3
    dC66 = -mu * e * U3

    # Apply perturbations (cracks with normal along x1)
    C = C_vti.copy()
    C[0, 0] += dC11
    C[0, 1] += dC12
    C[1, 0] += dC12
    C[0, 2] += dC13
    C[2, 0] += dC13
    C[4, 4] += dC55  # C55
    C[5, 5] += dC66  # C66

    return C


def build_layer_stiffness(layer_params: dict) -> np.ndarray:
    """
    Build stiffness matrix from layer parameters dictionary.

    Parameters
    ----------
    layer_params : dict
        Dictionary containing:
        - 'type': 'ISO', 'VTI', 'HTI', or 'OA'
        - 'vp': P-wave velocity (km/s)
        - 'vs': S-wave velocity (km/s)
        - 'rho': Density (g/cm^3)
        For VTI/OA:
        - 'epsilon': Thomsen epsilon
        - 'delta': Thomsen delta
        - 'gamma': Thomsen gamma
        For HTI/OA:
        - 'crack_density': Crack density parameter

    Returns
    -------
    C : ndarray
        6x6 stiffness matrix (GPa)
    """
    media_type = layer_params.get('type', 'ISO').upper()
    vp = layer_params['vp']
    vs = layer_params['vs']
    rho = layer_params['rho']

    if media_type == 'ISO':
        return stiffness_isotropic(vp, vs, rho)

    elif media_type == 'VTI':
        epsilon = layer_params.get('epsilon', 0.0)
        delta = layer_params.get('delta', 0.0)
        gamma = layer_params.get('gamma', 0.0)
        return stiffness_vti(vp, vs, rho, epsilon, delta, gamma)

    elif media_type == 'HTI':
        crack_density = layer_params.get('crack_density', 0.05)
        return stiffness_hti_hudson(vp, vs, rho, crack_density)

    elif media_type == 'OA':
        epsilon = layer_params.get('epsilon', 0.0)
        delta = layer_params.get('delta', 0.0)
        gamma = layer_params.get('gamma', 0.0)
        crack_density = layer_params.get('crack_density', 0.05)
        return stiffness_orthorhombic(vp, vs, rho, epsilon, delta, gamma, crack_density)

    else:
        raise ValueError(f"Unknown media type: {media_type}")


def build_model_stiffness(layers: list) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build stiffness matrices for a layered model.

    Parameters
    ----------
    layers : list
        List of layer parameter dictionaries

    Returns
    -------
    thickness : ndarray
        Layer thicknesses (km)
    density : ndarray
        Layer densities (g/cm^3)
    stiffness : ndarray
        Stiffness tensors with shape (nLayers, 6, 6)
    """
    n_layers = len(layers)
    thickness = np.zeros(n_layers)
    density = np.zeros(n_layers)
    stiffness = np.zeros((n_layers, 6, 6))

    for i, layer in enumerate(layers):
        thickness[i] = layer['thickness']
        density[i] = layer['rho']
        stiffness[i] = build_layer_stiffness(layer)

    return thickness, density, stiffness


# =============================================================================
# Preset models for quick testing
# =============================================================================

def preset_two_layer_vti() -> list:
    """Two-layer VTI model for quick testing."""
    return [
        {
            'type': 'VTI',
            'thickness': 0.5,
            'vp': 3.0,
            'vs': 1.5,
            'rho': 2.3,
            'epsilon': 0.05,
            'delta': 0.02,
            'gamma': 0.03
        },
        {
            'type': 'VTI',
            'thickness': 1.0,
            'vp': 3.5,
            'vs': 1.8,
            'rho': 2.5,
            'epsilon': 0.08,
            'delta': 0.04,
            'gamma': 0.05
        }
    ]


def preset_shale_sand() -> list:
    """Typical shale-sand model with VTI shale over isotropic sand."""
    return [
        {
            'type': 'VTI',
            'thickness': 0.3,
            'vp': 3.2,
            'vs': 1.6,
            'rho': 2.4,
            'epsilon': 0.15,
            'delta': 0.05,
            'gamma': 0.10
        },
        {
            'type': 'ISO',
            'thickness': 0.5,
            'vp': 3.8,
            'vs': 2.2,
            'rho': 2.3,
            'epsilon': 0.0,
            'delta': 0.0,
            'gamma': 0.0
        },
        {
            'type': 'VTI',
            'thickness': 0.8,
            'vp': 3.4,
            'vs': 1.7,
            'rho': 2.5,
            'epsilon': 0.12,
            'delta': 0.04,
            'gamma': 0.08
        }
    ]


def preset_fractured_reservoir() -> list:
    """Fractured reservoir model with HTI layer."""
    return [
        {
            'type': 'VTI',
            'thickness': 0.4,
            'vp': 3.0,
            'vs': 1.5,
            'rho': 2.4,
            'epsilon': 0.10,
            'delta': 0.03,
            'gamma': 0.06
        },
        {
            'type': 'HTI',
            'thickness': 0.2,
            'vp': 3.5,
            'vs': 2.0,
            'rho': 2.3,
            'crack_density': 0.08
        },
        {
            'type': 'VTI',
            'thickness': 0.6,
            'vp': 3.3,
            'vs': 1.7,
            'rho': 2.5,
            'epsilon': 0.08,
            'delta': 0.02,
            'gamma': 0.05
        }
    ]
