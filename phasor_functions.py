"""
Functions for spectral phasor analysis.

This module contains functions for calculating spectral phasor
coordinates from fluorescence emission spectra, generating
spectral phasor plots, visualizing fluorescence spectra, plotting
selected thermodynamic states, and saving figures.

The mathematical implementation of the spectral phasor
coordinates is kept consistent with the analysis used in the
associated publication.
"""


import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator


# -------------------------------------------------------------------
# Spectral phasor calculations
# -------------------------------------------------------------------

def G(lambdas, I):
    """
    Calculate the Gλ spectral phasor coordinate.

    Parameters
    ----------
    lambdas : array-like
        Wavelength values.

    I : array-like
        Fluorescence intensity values.

    Returns
    -------
    float
        Gλ spectral phasor coordinate.
    """

    min_lambda = np.min(lambdas)
    max_lambda = np.max(lambdas)

    return np.sum(
        (I / np.sum(I))
        * np.cos(
            2 * np.pi
            * (lambdas - min_lambda)
            / (max_lambda - min_lambda)
        )
    )


def S(lambdas, I):
    """
    Calculate the Sλ spectral phasor coordinate.

    Parameters
    ----------
    lambdas : array-like
        Wavelength values.

    I : array-like
        Fluorescence intensity values.

    Returns
    -------
    float
        Sλ spectral phasor coordinate.
    """

    min_lambda = np.min(lambdas)
    max_lambda = np.max(lambdas)

    return np.sum(
        (I / np.sum(I))
        * np.sin(
            2 * np.pi
            * (lambdas - min_lambda)
            / (max_lambda - min_lambda)
        )
    )


def theta(G_value, S_value):
    """
    Calculate the angular coordinate of a spectral phasor.

    Parameters
    ----------
    G_value : float
        Gλ spectral phasor coordinate.

    S_value : float
        Sλ spectral phasor coordinate.

    Returns
    -------
    float
        Phasor angle in degrees.
    """

    return np.degrees(
        np.arctan2(S_value, G_value) % (2 * np.pi)
    )


def M(G_value, S_value):
    """
    Calculate the magnitude of a spectral phasor.

    Parameters
    ----------
    G_value : float
        Gλ spectral phasor coordinate.

    S_value : float
        Sλ spectral phasor coordinate.

    Returns
    -------
    float
        Phasor magnitude.
    """

    return np.sqrt(
        G_value**2 + S_value**2
    )


# -------------------------------------------------------------------
# Plot formatting
# -------------------------------------------------------------------

def _format_axes(ax):
    """
    Apply the standard figure formatting used in the analysis.
    """

    for side in [
        'top',
        'right',
        'bottom',
        'left'
    ]:

        ax.spines[side].set_linewidth(2)

    ax.tick_params(
        axis='both',
        which='major',
        direction='in',
        length=6,
        width=1.5
    )

    ax.tick_params(
        axis='both',
        which='minor',
        direction='in',
        length=3,
        width=1.5
    )

    ax.tick_params(
        axis='both',
        labelsize=16,
        pad=10
    )

    for label in (
        ax.get_xticklabels()
        + ax.get_yticklabels()
    ):

        label.set_fontfamily(
            'sans-serif'
        )


# -------------------------------------------------------------------
# Fluorescence spectra
# -------------------------------------------------------------------

def create_spectra_plot(
    ax,
    wavelengths,
    spectra,
    labels=None,
    colormap='winter'
):
    """
    Plot fluorescence emission spectra.

    The first spectrum is displayed with a lighter color and the
    last spectrum with a more saturated color. The same color
    progression is used for the corresponding phasor points.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Matplotlib axis on which the spectra will be plotted.

    wavelengths : array-like
        Wavelength values.

    spectra : array-like
        Fluorescence spectra.

    labels : list, optional
        Labels associated with each spectrum.

    colormap : str, optional
        Matplotlib colormap used to represent the spectral progression.
        Default is 'Blues'.

    Returns
    -------
    matplotlib.axes.Axes
        The configured spectra plot axis.

    numpy.ndarray
        Colors assigned to each spectrum.
    """

    spectra = np.asarray(spectra)

    colors = plt.get_cmap(colormap)(
        np.linspace(0.25, 1.0, len(spectra))
    )

    for i, spectrum in enumerate(spectra):

        label = labels[i] if labels is not None else None

        ax.plot(
            wavelengths,
            spectrum,
            color=colors[i],
            linewidth=1.5,
            label=label
        )

    _format_axes(ax)

    ax.set_xlabel(
        'Wavelength (nm)',
        fontsize=20,
        labelpad=10
    )

    ax.set_ylabel(
        'Fluorescence intensity',
        fontsize=20,
        labelpad=10
    )

    ax.xaxis.set_minor_locator(
        AutoMinorLocator(2)
    )

    ax.yaxis.set_minor_locator(
        AutoMinorLocator(2)
    )

    return ax, colors


# -------------------------------------------------------------------
# Spectral phasor background
# -------------------------------------------------------------------

def create_phasor_plot(
    ax,
    universal_circle=False,
    grid=True
):
    """
    Configure an axis for spectral phasor visualization.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Matplotlib axis.

    universal_circle : bool, optional
        If True, display the universal semicircle.

    grid : bool, optional
        If True, display the phasor grid.

    Returns
    -------
    matplotlib.axes.Axes
        Configured phasor plot axis.
    """

    ax.clear()

    _format_axes(
        ax
    )

    t = np.linspace(
        0,
        2 * np.pi,
        201
    )

    if grid:

        radii = [
            1.0,
            0.8,
            0.6,
            0.4,
            0.2
        ]

        for r in radii:

            x = r * np.cos(t)
            y = r * np.sin(t)

            linewidth = (
                2
                if r == 1.0
                else 1
            )

            ax.plot(
                x,
                y,
                '-',
                color='gray',
                alpha=0.5,
                linewidth=linewidth
            )

        n_lines = 24

        angles = np.linspace(
            0,
            2 * np.pi,
            n_lines,
            endpoint=False
        )

        for angle in angles:

            x = np.cos(angle)
            y = np.sin(angle)

            ax.plot(
                [0, x],
                [0, y],
                '--',
                color='gray',
                alpha=0.4,
                linewidth=1
            )

    if universal_circle:

        phi = np.linspace(
            0,
            np.pi / 2,
            201
        )

        x_semi = (
            np.cos(phi)**2
        )

        y_semi = (
            np.cos(phi)
            * np.sin(phi)
        )

        ax.plot(
            x_semi,
            y_semi,
            color='orange',
            linewidth=2
        )

    ax.set_xlabel(
        r'$G_{\lambda}$',
        fontsize=22
    )

    ax.set_ylabel(
        r'$S_{\lambda}$',
        fontsize=22
    )

    ax.xaxis.set_minor_locator(
        AutoMinorLocator(1)
    )

    ax.yaxis.set_minor_locator(
        AutoMinorLocator(1)
    )

    ax.set_aspect(
        'equal',
        adjustable='box'
    )

    return ax


# -------------------------------------------------------------------
# Phasor plot limits
# -------------------------------------------------------------------

def set_phasor_limits(
    ax,
    G_values,
    S_values,
    padding=0.15
):
    """
    Set the phasor plot limits around the experimental points.
    """

    G_values = np.asarray(
        G_values
    )

    S_values = np.asarray(
        S_values
    )

    G_min = np.min(
        G_values
    )

    G_max = np.max(
        G_values
    )

    S_min = np.min(
        S_values
    )

    S_max = np.max(
        S_values
    )

    G_range = (
        G_max - G_min
    )

    S_range = (
        S_max - S_min
    )

    data_range = max(
        G_range,
        S_range
    )

    if data_range == 0:
        data_range = 0.1

    margin = (
        data_range
        * padding
    )

    center_G = (
        G_max + G_min
    ) / 2

    center_S = (
        S_max + S_min
    ) / 2

    half_range = (
        data_range / 2
        + margin
    )

    ax.set_xlim(
        center_G - half_range,
        center_G + half_range
    )

    ax.set_ylim(
        center_S - half_range,
        center_S + half_range
    )


# -------------------------------------------------------------------
# Phasor trajectory
# -------------------------------------------------------------------

def plot_phasor_trajectory(
    ax,
    G_values,
    S_values,
    colors,
    marker_size=100,
    line_width=2
):
    """
    Plot experimental spectral phasors and their trajectory.
    """

    G_values = np.asarray(
        G_values
    )

    S_values = np.asarray(
        S_values
    )

    ax.plot(
        G_values,
        S_values,
        color=colors[-1],
        linewidth=line_width,
        alpha=0.7
    )

    ax.scatter(
        G_values,
        S_values,
        c=colors,
        s=marker_size,
        edgecolors='none'
    )

    set_phasor_limits(
        ax,
        G_values,
        S_values
    )


# -------------------------------------------------------------------
# Selected thermodynamic states
# -------------------------------------------------------------------

def plot_state_connections(
    ax,
    G_values,
    S_values,
    state_indices,
    color='black',
    line_width=2.5
):
    """
    Connect selected thermodynamic states sequentially.

    The states are connected in the order provided by
    state_indices.

    For n selected states, n-1 lines are drawn.

    Examples
    --------
    Two states:

        initial -------- final

    Three states:

        initial ---- intermediate ---- final

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Matplotlib axis.

    G_values : array-like
        Gλ phasor coordinates.

    S_values : array-like
        Sλ phasor coordinates.

    state_indices : list of int
        Indices of the selected states in experimental order.

    color : str, optional
        Line color.

    line_width : float, optional
        Line width.
    """

    G_values = np.asarray(
        G_values
    )

    S_values = np.asarray(
        S_values
    )

    for i in range(
        len(state_indices) - 1
    ):

        idx1 = state_indices[i]
        idx2 = state_indices[i + 1]

        ax.plot(
            [
                G_values[idx1],
                G_values[idx2]
            ],
            [
                S_values[idx1],
                S_values[idx2]
            ],
            color=color,
            linewidth=line_width,
            zorder=3
        )


# -------------------------------------------------------------------
# Thermodynamic analysis
# -------------------------------------------------------------------


def distance(
    G1,
    S1,
    G2,
    S2
):
    """
    Calculate the Euclidean distance between two spectral phasors.
    """

    return np.sqrt(
        (G2 - G1)**2
        +
        (S2 - S1)**2
    )


def relative_distance(
    G1,
    G2,
    Gi,
    S1,
    S2,
    Si
):
    """
    Calculate the relative distance of an experimental phasor
    along a selected thermodynamic transition.

    The distance is normalized to the distance between the
    initial and final thermodynamic states.
    """

    return (
        distance(
            G1,
            S1,
            Gi,
            Si
        )
        /
        distance(
            G1,
            S1,
            G2,
            S2
        )
    )


def exp_phasor(
    x,
    I_initial,
    I_final,
    dG,
    m,
    T,
    R
):
    """
    Two-state thermodynamic model for the relative spectral
    phasor distance.

    Parameters
    ----------
    x : array-like
        Experimental variable.

    I_initial : float
        Integral of the initial-state spectrum.

    I_final : float
        Integral of the final-state spectrum.

    dG : float
        Standard free energy difference.

    m : float
        Dependence of the free energy on the experimental variable.

    T : float
        Temperature in Kelvin.

    R : float
        Gas constant in units consistent with dG.

    Returns
    -------
    numpy.ndarray
        Calculated relative distance.
    """

    exp_term = np.exp(
        (dG - m * x)
        /
        (R * T)
    )

    return (
        I_final
        /
        (
            I_final
            +
            I_initial * exp_term
        )
    )


# -------------------------------------------------------------------
# Thermodynamic plot
# -------------------------------------------------------------------


def plot_thermodynamic_graph(
    x,
    y,
    ax,
    color='#0084a5',
    xlabel='x',
    ylabel='Relative distance',
    tipo='line'
):
    """
    Plot relative spectral phasor distance as a function
    of the experimental variable.
    """

    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)

    if tipo == 'line':

        ax.plot(
            x,
            y,
            color=color,
            linewidth=2.5
        )

    else:

        ax.scatter(
            x,
            y,
            color=color,
            s=120
        )

    ax.set_xlabel(
        xlabel,
        fontsize=20,
        fontname='sans-serif',
        labelpad=10
    )

    ax.set_ylabel(
        ylabel,
        fontsize=20,
        fontname='sans-serif',
        labelpad=10
    )

    ax.tick_params(
        axis='both',
        labelsize=16,
        pad=10
    )

    ax.xaxis.set_minor_locator(
        AutoMinorLocator(2)
    )

    ax.yaxis.set_minor_locator(
        AutoMinorLocator(2)
    )

    ax.tick_params(
        axis='both',
        which='major',
        direction='in',
        length=6,
        width=1.5
    )

    ax.tick_params(
        axis='both',
        which='minor',
        direction='in',
        length=3,
        width=1.5
    )

    return ax
# -------------------------------------------------------------------
# Save figures
# -------------------------------------------------------------------

def save_figures(
    fig_spectra,
    fig_phasor,
    output_directory,
    base_name='spectral_phasor'
):
    """
    Save the fluorescence spectra and spectral phasor figures.

    Parameters
    ----------
    fig_spectra : matplotlib.figure.Figure
        Fluorescence spectra figure.

    fig_phasor : matplotlib.figure.Figure
        Fluorescence phasor figure.

    output_directory : pathlib.Path
        Output directory.

    base_name : str, optional
        Base name for the figures.
    """

    output_directory.mkdir(
        parents=True,
        exist_ok=True
    )

    fig_spectra.savefig(
        output_directory
        / f'{base_name}_spectra.png',
        dpi=300,
        bbox_inches='tight'
    )

    fig_phasor.savefig(
        output_directory
        / f'{base_name}_phasor.png',
        dpi=300,
        bbox_inches='tight'
    )



