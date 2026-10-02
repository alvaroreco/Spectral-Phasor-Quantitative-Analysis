"""
Spectral phasor analysis and plotting.

This script:
    - Loads fluorescence spectra.
    - Calculates spectral phasor coordinates.
    - Displays the spectra and phasor trajectory.
    - Optionally saves G/S data as CSV or TXT.
    - Optionally saves figures.

All generated files are stored in the 'output' directory.

Spectrum labels are used only as identifiers in this analysis.
They do not need to be numerical.
"""

from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# PROJECT PATH
# =============================================================================

PROJECT_DIR = Path(__file__).resolve().parent

sys.path.insert(
    0,
    str(PROJECT_DIR)
)


# =============================================================================
# IMPORTS
# =============================================================================

from data_reader import load_data_with_condition

from phasor_functions import (
    G,
    S,
    create_spectra_plot,
    create_phasor_plot,
    plot_phasor_trajectory,
)


# =============================================================================
# LOAD DATA
# =============================================================================

(
    wavelengths,
    spectra,
    labels,
    condition
) = load_data_with_condition()


# =============================================================================
# CALCULATE SPECTRAL PHASORS
# =============================================================================

G_values = np.array(
    [
        G(wavelengths, spectrum)
        for spectrum in spectra
    ]
)

S_values = np.array(
    [
        S(wavelengths, spectrum)
        for spectrum in spectra
    ]
)


# =============================================================================
# SPECTRA FIGURE
# =============================================================================

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


# =============================================================================
# PHASOR FIGURE
# =============================================================================

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


# =============================================================================
# DISPLAY FIGURES
# =============================================================================

plt.show(block=False)
plt.pause(0.1)


# =============================================================================
# OUTPUT DIRECTORY
# =============================================================================

OUTPUT_DIR = PROJECT_DIR / "output"

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# =============================================================================
# SAVE OPTIONS
# =============================================================================

save_data_option = input(
    "\nSave spectral phasor data (G and S)? [y/n]: "
).strip().lower()

save_figures_option = input(
    "\nSave figures to the 'output' folder? [y/n]: "
).strip().lower()


# =============================================================================
# OUTPUT FILE NAME
# =============================================================================

save_data = save_data_option == "y"
save_figures = save_figures_option == "y"

if save_data or save_figures:

    file_name = input(
        "\nEnter a name for the output files: "
    ).strip()

    if not file_name:

        file_name = "spectral_phasor_data"

    file_name = Path(file_name).stem


# =============================================================================
# SAVE PHASOR DATA
# =============================================================================

if save_data:

    output_format = input(
        "\nSelect output format: "
        "[1] CSV  [2] TXT (tab-delimited): "
    ).strip()

    # The first column contains the spectrum label.
    # No numerical x-column is included here because this script
    # does not perform quantitative thermodynamic analysis.
    header = [
        condition,
        "G",
        "S"
    ]

    if output_format == "1":

        output_file = (
            OUTPUT_DIR
            / f"{file_name}_phasor.csv"
        )

        with open(
            output_file,
            "w",
            encoding="utf-8"
        ) as file:

            file.write(
                ",".join(header)
                + "\n"
            )

            for label, g, s in zip(
                labels,
                G_values,
                S_values
            ):

                file.write(
                    f"{label},{g},{s}\n"
                )

        print(
            "\nPhasor data saved as:"
        )

        print(output_file)

    elif output_format == "2":

        output_file = (
            OUTPUT_DIR
            / f"{file_name}_phasor.txt"
        )

        with open(
            output_file,
            "w",
            encoding="utf-8"
        ) as file:

            file.write(
                "\t".join(header)
                + "\n"
            )

            for label, g, s in zip(
                labels,
                G_values,
                S_values
            ):

                file.write(
                    f"{label}\t{g}\t{s}\n"
                )

        print(
            "\nPhasor data saved as:"
        )

        print(output_file)

    else:

        print(
            "\nInvalid output format. "
            "Phasor data were not saved."
        )

elif save_data_option == "n":

    print(
        "\nPhasor data were not saved."
    )

else:

    print(
        "\nInvalid option selected. "
        "Phasor data were not saved."
    )


# =============================================================================
# SAVE FIGURES
# =============================================================================

if save_figures:

    spectra_file = (
        OUTPUT_DIR
        / f"{file_name}_spectra.png"
    )

    phasor_file = (
        OUTPUT_DIR
        / f"{file_name}_phasor.png"
    )

    fig_spectra.savefig(
        spectra_file,
        dpi=300,
        bbox_inches="tight"
    )

    fig_phasor.savefig(
        phasor_file,
        dpi=300,
        bbox_inches="tight"
    )

    print(
        "\nFigures saved in the 'output' folder:"
    )

    print(spectra_file)
    print(phasor_file)

elif save_figures_option == "n":

    print(
        "\nFigures were not saved."
    )

else:

    print(
        "\nInvalid option selected. "
        "Figures were not saved."
    )


# =============================================================================
# KEEP FIGURES OPEN
# =============================================================================

plt.show()

