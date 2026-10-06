"""
Spectral phasor analysis and plotting.

This script:
    - Loads fluorescence spectra.
    - Calculates spectral phasor coordinates.
    - Displays the spectra and phasor trajectory.
    - Optionally saves G/S data.
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
# ERROR REPORTING
# =============================================================================

def report_error(error):
    """Display a user-friendly error message."""

    print("\n" + "=" * 70)
    print(
        "ERROR: The spectral phasor analysis "
        "could not be completed."
    )
    print("=" * 70)

    print("\n" + str(error))

    print("\nPlease check the following:")
    print("  - The selected data folder is correct.")
    print("  - The selected file format matches the data.")
    print("  - The selected files contain spectral data.")
    print("  - The number of header rows is correct.")
    print("  - Wavelength values are numerical.")
    print("  - Fluorescence intensity values are numerical.")
    print("  - All selected spectra use the same wavelength axis.")
    print("  - At least two spectra are available.")
    print("  - The spectra contain finite intensity values.")
    print("  - The output folder can be written if saving is requested.")

    print("\nTechnical information:")
    print(f"  Error type: {type(error).__name__}")

    print("\nNo further analysis was performed.")
    print("=" * 70)


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def main():

    # -------------------------------------------------------------------------
    # Load data
    # -------------------------------------------------------------------------

    (
        wavelengths,
        spectra,
        labels,
        condition
    ) = load_data_with_condition()

    # -------------------------------------------------------------------------
    # Calculate spectral phasors
    # -------------------------------------------------------------------------

    G_values = np.asarray(
        [
            G(
                wavelengths,
                spectrum
            )
            for spectrum in spectra
        ],
        dtype=float
    )

    S_values = np.asarray(
        [
            S(
                wavelengths,
                spectrum
            )
            for spectrum in spectra
        ],
        dtype=float
    )

    if not np.all(np.isfinite(G_values)):

        raise ValueError(
            "The calculated Gλ phasor coordinates contain "
            "non-finite values.\n\n"
            "Possible causes:\n"
            "  - A spectrum has a total intensity of zero.\n"
            "  - The wavelength range is zero.\n"
            "  - The spectral data contain invalid values."
        )

    if not np.all(np.isfinite(S_values)):

        raise ValueError(
            "The calculated Sλ phasor coordinates contain "
            "non-finite values.\n\n"
            "Possible causes:\n"
            "  - A spectrum has a total intensity of zero.\n"
            "  - The wavelength range is zero.\n"
            "  - The spectral data contain invalid values."
        )

    # -------------------------------------------------------------------------
    # Spectra figure
    # -------------------------------------------------------------------------

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

    # -------------------------------------------------------------------------
    # Phasor figure
    # -------------------------------------------------------------------------

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

    # -------------------------------------------------------------------------
    # Display figures
    # -------------------------------------------------------------------------

    plt.show(
        block=False
    )

    plt.pause(
        0.1
    )

    # -------------------------------------------------------------------------
    # Output directory
    # -------------------------------------------------------------------------

    output_directory = (
        PROJECT_DIR
        / "output"
    )

    try:

        output_directory.mkdir(
            parents=True,
            exist_ok=True
        )

    except OSError as error:

        raise OSError(
            "The 'output' directory could not be created.\n\n"
            f"Path:\n{output_directory}\n\n"
            "Please check that you have write permission "
            "for the project directory."
        ) from error

    # -------------------------------------------------------------------------
    # Save options
    # -------------------------------------------------------------------------

    save_data_option = input(
        "\nSave spectral phasor data (G and S)? [y/n]: "
    ).strip().lower()

    save_figures_option = input(
        "\nSave figures to the 'output' folder? [y/n]: "
    ).strip().lower()

    save_data = save_data_option == "y"
    save_figures = save_figures_option == "y"

    # -------------------------------------------------------------------------
    # Output file name
    # -------------------------------------------------------------------------

    if save_data or save_figures:

        file_name = input(
            "\nEnter a name for the output files: "
        ).strip()

        if not file_name:
            file_name = "spectral_phasor_data"

        file_name = Path(
            file_name
        ).stem

    # -------------------------------------------------------------------------
    # Save phasor data
    # -------------------------------------------------------------------------

    if save_data:

        output_format = input(
            "\nSelect output format: "
            "[1] CSV  [2] TXT (tab-delimited): "
        ).strip()

        header = [
            condition,
            "G",
            "S"
        ]

        try:

            if output_format == "1":

                output_file = (
                    output_directory
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

            elif output_format == "2":

                output_file = (
                    output_directory
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

            else:

                print(
                    "\nInvalid output format. "
                    "Phasor data were not saved."
                )

                output_file = None

            if output_file is not None:

                print(
                    "\nPhasor data saved as:"
                )

                print(
                    output_file
                )

        except OSError as error:

            raise OSError(
                "The phasor data could not be saved.\n\n"
                f"Output file:\n{output_file}\n\n"
                "Possible causes:\n"
                "  - The file is open in another program.\n"
                "  - You do not have write permission.\n"
                "  - The output directory cannot be accessed."
            ) from error

    elif save_data_option == "n":

        print(
            "\nPhasor data were not saved."
        )

    else:

        print(
            "\nInvalid option selected. "
            "Phasor data were not saved."
        )

    # -------------------------------------------------------------------------
    # Save figures
    # -------------------------------------------------------------------------

    if save_figures:

        spectra_file = (
            output_directory
            / f"{file_name}_spectra.png"
        )

        phasor_file = (
            output_directory
            / f"{file_name}_phasor.png"
        )

        try:

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

        except OSError as error:

            raise OSError(
                "The figures could not be saved.\n\n"
                f"Output directory:\n{output_directory}\n\n"
                "Possible causes:\n"
                "  - The output file is open in another program.\n"
                "  - You do not have write permission.\n"
                "  - The output directory cannot be accessed."
            ) from error

        print(
            "\nFigures saved in the 'output' folder:"
        )

        print(
            spectra_file
        )

        print(
            phasor_file
        )

    elif save_figures_option == "n":

        print(
            "\nFigures were not saved."
        )

    else:

        print(
            "\nInvalid option selected. "
            "Figures were not saved."
        )

    # -------------------------------------------------------------------------
    # Keep figures open
    # -------------------------------------------------------------------------

    plt.show()


# =============================================================================
# MAIN GUARD
# =============================================================================

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        print(
            "\n\nAnalysis interrupted by the user."
        )

    except (
        FileNotFoundError,
        NotADirectoryError,
        ValueError,
        RuntimeError,
        OSError
    ) as error:

        report_error(
            error
        )

    except Exception as error:

        print("\n" + "=" * 70)
        print("UNEXPECTED ERROR")
        print("=" * 70)

        print(
            f"\nError type: {type(error).__name__}"
        )

        print(
            f"Message: {error}"
        )

        print(
            "\nPlease report this message if the problem persists."
        )

        print("=" * 70)