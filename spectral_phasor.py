import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


def G(lambdas, I):
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
    return np.degrees(
        np.arctan2(S_value, G_value) % (2 * np.pi)
    )


def M(G_value, S_value):
    return np.sqrt(G_value**2 + S_value**2)


def create_phasor_plot(ax, universal_circle=False, grid=True):

    ax.clear()

    for side in ['top', 'right', 'bottom', 'left']:
        ax.spines[side].set_linewidth(2)

    t = np.linspace(0, 2 * np.pi, 201)

    if grid:

        radii = [1.0, 0.8, 0.6, 0.4, 0.2]

        for r in radii:

            x = r * np.cos(t)
            y = r * np.sin(t)

            linewidth = 2 if r == 1.0 else 1

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

        phi = np.linspace(0, np.pi / 2, 201)

        x_semi = np.cos(phi)**2
        y_semi = np.cos(phi) * np.sin(phi)

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

    ax.xaxis.set_minor_locator(AutoMinorLocator(1))
    ax.yaxis.set_minor_locator(AutoMinorLocator(1))

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

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontfamily('sans-serif')

    ax.set_aspect('equal', adjustable='box')

    return ax