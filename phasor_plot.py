"""
Generate fluorescence spectra and spectral phasor plots.
"""


from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt


# Project directory
PROJECT_DIR = Path(__file__).resolve().parent

sys.path.insert(0, str(PROJECT_DIR))


# Load data
from data_reader import load_data


# Phasor functions
from phasor_functions import (
    G,
    S,
    create_spectra_plot,
    create_phasor_plot,
    plot_phasor_trajectory,
    save_figures
)


# -------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------

wavelengths, spectra, labels = load_data()


# -------------------------------------------------------------------
# Calculate spectral phasors
# -------------------------------------------------------------------

G_values = np.array([
    G(wavelengths, spectrum)
    for spectrum in spectra
])

S_values = np.array([
    S(wavelengths, spectrum)
    for spectrum in spectra
])


# -------------------------------------------------------------------
# Create fluorescence spectra plot
# -------------------------------------------------------------------

fig_spectra, ax_spectra = plt.subplots(
    figsize=(8, 6)
)

_, colors = create_spectra_plot(
    ax_spectra,
    wavelengths,
    spectra,
    labels=labels
)

fig_spectra.tight_layout()


# -------------------------------------------------------------------
# Create spectral phasor plot
# -------------------------------------------------------------------

fig_phasor, ax_phasor = plt.subplots(
    figsize=(7, 7)
)

create_phasor_plot(
    ax_phasor,
    universal_circle=False,
    grid=True
)

plot_phasor_trajectory(
    ax_phasor,
    G_values,
    S_values,
    colors
)

fig_phasor.tight_layout()

# -------------------------------------------------------------------
# Display figures
# -------------------------------------------------------------------

plt.show(block=False)


# -------------------------------------------------------------------
# Optionally save figures
# -------------------------------------------------------------------

save = input(
    "\nSave figures to the 'Figures' folder? [y/n]: "
).strip().lower()

if save == 'y':

    base_name = input(
        "Enter a name for the figures: "
    ).strip()

    if not base_name:
        base_name = 'spectral_phasor'

    save_figures(
        fig_spectra,
        fig_phasor,
        PROJECT_DIR / 'Figures',
        base_name=base_name
    )

    print(
        f"\nFigures saved in the 'Figures' folder as "
        f"'{base_name}_spectra.png' and '{base_name}_phasor.png'."
    )

plt.show()