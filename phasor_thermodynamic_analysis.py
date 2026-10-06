"""
Spectral phasor analysis for thermodynamic denaturation.

This script:

    - Loads fluorescence spectra using data_reader.py.
    - Supports spectra stored either:
        1. one spectrum per file, or
        2. multiple spectra in a single file.
    - Uses numerical spectrum labels as experimental x-values.
    - Requires numerical spectrum labels for thermodynamic analysis.
    - Calculates spectral phasor coordinates.
    - Allows the user to select thermodynamic states.
    - Calculates relative spectral phasor distances between
      consecutive thermodynamic states.
    - Performs thermodynamic fitting of each transition.
    - Uses the experimental variable name supplied by data_reader
      as the x-axis label.
    - Optionally saves figures in the 'output' directory.
    - Figures displayed on screen have titles.
    - Saved figures do not have titles.
    - Optionally exports analysis data and thermodynamic results.
    - Does not export smooth fitted curves as separate files.

Data loading is handled exclusively by data_reader.py.

General spectral phasor and thermodynamic functions are handled
by phasor_functions.py.

The thermodynamic model is:

    r(x) = I_final /
           (I_final + I_initial * exp((dG - m*x)/(R*T)))

where:

    dG = standard free energy difference
    m  = dependence of the free energy on the experimental variable
"""


from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# ======================================================================
# PROJECT PATH
# ======================================================================

PROJECT_DIR = Path(__file__).resolve().parent

if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))


# ======================================================================
# IMPORT REPOSITORY MODULES
# ======================================================================

from data_reader import load_data_with_x_values

from phasor_functions import (
    G,
    S,
    distance,
    relative_distance,
    exp_phasor,
    create_spectra_plot,
    create_phasor_plot,
    plot_phasor_trajectory,
    plot_state_connections,
    plot_thermodynamic_graph,
)


# ======================================================================
# GENERAL SETTINGS
# ======================================================================

UNIVERSAL_CIRCLE = False
GRID = True


# ======================================================================
# ERROR REPORTING
# ======================================================================

def report_error(error):
    """Display a user-friendly error message."""

    print("\n" + "=" * 70)
    print("ERROR: The thermodynamic analysis could not be completed.")
    print("=" * 70)

    print("\n" + str(error))

    print("\nPlease check the following:")
    print("  - The selected folder and files are correct.")
    print("  - The selected file format matches the data.")
    print("  - The number of header rows is correct.")
    print("  - Thermodynamic spectrum labels are numerical.")
    print("    Examples: 0.000, 0.026, 0.051")
    print("  - At least two spectra are available.")
    print("  - Initial and final thermodynamic states are different.")
    print("  - States were selected in increasing experimental order.")
    print("  - Each thermodynamic transition contains enough")
    print("    experimental points for fitting.")
    print("  - Temperature is greater than zero.")
    print("  - The selected energy units are correct.")
    print("  - The output folder can be written.")

    print("\nTechnical information:")
    print(f"  Error type: {type(error).__name__}")

    print("\nNo further analysis was performed.")
    print("=" * 70)


# ======================================================================
# THERMODYNAMIC STATE SELECTION
# ======================================================================

def select_thermodynamic_states(number_of_spectra):
    """
    Select the thermodynamic states used for the analysis.

    The first and last spectra are always selected.
    Intermediate states can be selected by spectrum number.
    """

    if number_of_spectra < 2:
        raise ValueError(
            "At least two spectra are required for thermodynamic analysis."
        )

    print("\nSelect the thermodynamic states.")
    print("The first and last spectra are automatically selected.")
    print(
        "You can select intermediate states by entering their "
        "spectrum numbers in increasing order.\n"
    )

    print(
        f"Number of spectra available: "
        f"{number_of_spectra}"
    )

    if number_of_spectra == 2:
        return [0, 1]

    print(
        f"\nEnter the intermediate state spectrum numbers "
        f"(2 to {number_of_spectra - 1})."
    )

    print("Separate multiple numbers with commas.")
    print("Press Enter if there are no intermediate states.")

    user_input = input(
        "\nIntermediate states: "
    ).strip()

    if user_input == "":
        return [
            0,
            number_of_spectra - 1
        ]

    try:
        intermediate_numbers = [
            int(value.strip())
            for value in user_input.split(",")
        ]

    except ValueError as error:
        raise ValueError(
            "Invalid state selection.\n"
            "Please enter spectrum numbers separated by commas."
        ) from error

    if len(set(intermediate_numbers)) != len(
        intermediate_numbers
    ):
        raise ValueError(
            "Intermediate states cannot be repeated."
        )

    for number in intermediate_numbers:

        if (
            number < 2
            or number > number_of_spectra - 1
        ):
            raise ValueError(
                f"Intermediate state {number} "
                "is outside the allowed range."
            )

    if intermediate_numbers != sorted(
        intermediate_numbers
    ):
        raise ValueError(
            "Intermediate states must be entered "
            "in increasing order."
        )

    # Convert from spectrum numbers (1-based)
    # to Python indices (0-based).

    state_indices = [0]

    for number in intermediate_numbers:
        state_indices.append(number - 1)

    state_indices.append(
        number_of_spectra - 1
    )

    return state_indices


# ======================================================================
# TEMPERATURE AND ENERGY UNITS
# ======================================================================

def select_temperature_and_units():
    """Select temperature and energy units."""

    print("\nThermodynamic parameters")
    print("------------------------")

    temperature_input = input(
        "Temperature (K): "
    ).strip()

    try:
        temperature = float(
            temperature_input
        )

    except ValueError as error:
        raise ValueError(
            "Temperature must be a numerical value.\n"
            "Example: 298.15"
        ) from error

    if not np.isfinite(temperature):
        raise ValueError(
            "Temperature must be a finite numerical value."
        )

    if temperature <= 0:
        raise ValueError(
            "Temperature must be greater than zero."
        )

    print("\nSelect energy units for ΔG:")
    print("1 - J/mol")
    print("2 - kJ/mol")
    print("3 - cal/mol")
    print("4 - kcal/mol")

    unit_choice = input(
        "\nChoice: "
    ).strip()

    if unit_choice == "1":

        gas_constant = 8.314462618
        energy_unit = "J/mol"

    elif unit_choice == "2":

        gas_constant = 0.008314462618
        energy_unit = "kJ/mol"

    elif unit_choice == "3":

        gas_constant = 1.987204258
        energy_unit = "cal/mol"

    elif unit_choice == "4":

        gas_constant = 0.001987204258
        energy_unit = "kcal/mol"

    else:

        raise ValueError(
            "Invalid energy-unit choice.\n"
            "Please select 1, 2, 3, or 4."
        )

    print(
        f"\nTemperature: {temperature:.2f} K"
    )

    print(
        f"Energy units: {energy_unit}"
    )

    print(
        f"Gas constant: {gas_constant}"
    )

    return (
        temperature,
        gas_constant,
        energy_unit
    )


# ======================================================================
# SORT DATA BY EXPERIMENTAL X-VALUE
# ======================================================================

def sort_spectral_data_by_x(
    wavelengths,
    spectra,
    labels,
    x_values
):
    """
    Sort spectra according to increasing experimental x-value.

    This is necessary because file names or file-column order may
    not correspond to the numerical order of the experimental variable.
    """

    x_values = np.asarray(
        x_values,
        dtype=float
    )

    if x_values.ndim != 1:
        raise ValueError(
            "Experimental x-values must form a one-dimensional array."
        )

    if not np.all(
        np.isfinite(x_values)
    ):
        raise ValueError(
            "The experimental x-values contain non-finite values."
        )

    if len(x_values) != len(spectra):
        raise ValueError(
            "The number of experimental x-values does not match "
            "the number of spectra."
        )

    order = np.argsort(
        x_values
    )

    x_values = x_values[order]

    spectra = np.asarray(
        spectra
    )[order]

    labels = [
        labels[index]
        for index in order
    ]

    return (
        wavelengths,
        spectra,
        labels,
        x_values
    )


# ======================================================================
# ANALYZE ONE TRANSITION
# ======================================================================

def calculate_transition(
    transition_number,
    initial_index,
    final_index,
    wavelengths,
    spectra,
    G_values,
    S_values,
    x_values,
    temperature,
    gas_constant
):
    """
    Calculate relative phasor distances and fit one thermodynamic
    transition.
    """

    # ------------------------------------------------------------------
    # Experimental points belonging to this transition
    # ------------------------------------------------------------------

    # The final state is excluded from the experimental points,
    # preserving the original analysis convention.

    indices = np.arange(
        initial_index,
        final_index
    )

    x_transition = np.asarray(
        x_values[indices],
        dtype=float
    )

    # ------------------------------------------------------------------
    # Validate experimental x-values
    # ------------------------------------------------------------------

    if len(x_transition) < 2:

        raise ValueError(
            f"Transition {transition_number} does not contain "
            "enough experimental points for thermodynamic fitting.\n"
            f"Initial state: spectrum {initial_index + 1} "
            f"(x = {x_values[initial_index]})\n"
            f"Final state: spectrum {final_index + 1} "
            f"(x = {x_values[final_index]})\n\n"
            "At least two experimental points are required "
            "to fit ΔG and m.\n"
            "Select an intermediate spectrum between these "
            "states or select thermodynamic states farther apart."
        )

    if not np.all(
        np.isfinite(x_transition)
    ):
        raise ValueError(
            f"Transition {transition_number} contains "
            "non-finite experimental x-values."
        )

    if len(
        np.unique(x_transition)
    ) < 2:
        raise ValueError(
            f"Transition {transition_number} contains "
            "repeated experimental x-values.\n"
            "The fitting procedure requires at least "
            "two different x-values."
        )

    # ------------------------------------------------------------------
    # Initial and final phasor coordinates
    # ------------------------------------------------------------------

    G_initial = G_values[
        initial_index
    ]

    S_initial = S_values[
        initial_index
    ]

    G_final = G_values[
        final_index
    ]

    S_final = S_values[
        final_index
    ]

    # ------------------------------------------------------------------
    # Validate initial/final phasor separation
    # ------------------------------------------------------------------

    state_distance = distance(
        G_initial,
        S_initial,
        G_final,
        S_final
    )

    if not np.isfinite(state_distance):

        raise ValueError(
            f"Transition {transition_number}: the distance between "
            "the initial and final phasor states is not finite."
        )

    if state_distance == 0:

        raise ValueError(
            f"Transition {transition_number}: the initial and final "
            "thermodynamic states have identical phasor coordinates.\n"
            "The relative phasor distance cannot be calculated "
            "for this transition."
        )

    # ------------------------------------------------------------------
    # Relative phasor distances
    # ------------------------------------------------------------------

    r_values = np.array(
        [
            relative_distance(
                G_initial,
                G_final,
                G_values[index],
                S_initial,
                S_final,
                S_values[index]
            )
            for index in indices
        ],
        dtype=float
    )

    if not np.all(
        np.isfinite(r_values)
    ):
        raise ValueError(
            f"Transition {transition_number} produced "
            "non-finite relative phasor distances.\n"
            "Check the selected initial and final states "
            "and the spectral data."
        )

    # ------------------------------------------------------------------
    # Integrated fluorescence
    # ------------------------------------------------------------------

    try:

        I_initial = np.trapz(
            spectra[initial_index],
            wavelengths
        )

        I_final = np.trapz(
            spectra[final_index],
            wavelengths
        )

    except (TypeError, ValueError) as error:

        raise ValueError(
            f"Transition {transition_number}: integrated fluorescence "
            "could not be calculated from the selected spectra."
        ) from error

    if not np.isfinite(I_initial):

        raise ValueError(
            f"Transition {transition_number}: the integrated "
            "fluorescence of the initial spectrum is not finite."
        )

    if not np.isfinite(I_final):

        raise ValueError(
            f"Transition {transition_number}: the integrated "
            "fluorescence of the final spectrum is not finite."
        )

    if I_initial == 0:

        raise ValueError(
            f"Transition {transition_number}: the integrated "
            "fluorescence of the initial spectrum is zero.\n"
            "The thermodynamic model cannot be fitted."
        )

    if I_final == 0:

        raise ValueError(
            f"Transition {transition_number}: the integrated "
            "fluorescence of the final spectrum is zero.\n"
            "The thermodynamic model cannot be fitted."
        )

    # ------------------------------------------------------------------
    # Initial guesses according to energy units
    # ------------------------------------------------------------------

    if np.isclose(
        gas_constant,
        8.314462618
    ):

        p0 = [
            10000,
            20000
        ]

    elif np.isclose(
        gas_constant,
        0.008314462618
    ):

        p0 = [
            10,
            20
        ]

    elif np.isclose(
        gas_constant,
        1.987204258
    ):

        p0 = [
            2390,
            4780
        ]

    elif np.isclose(
        gas_constant,
        0.001987204258
    ):

        p0 = [
            2.39,
            4.78
        ]

    else:

        raise ValueError(
            "The gas constant is not recognized."
        )

    # ------------------------------------------------------------------
    # Thermodynamic fit
    # ------------------------------------------------------------------

    def fit_function(
        x,
        dG,
        m
    ):
        """Wrapper around the thermodynamic model for curve_fit."""

        return exp_phasor(
            x,
            I_initial,
            I_final,
            dG,
            m,
            temperature,
            gas_constant
        )

    try:

        popt, pcov = curve_fit(
            fit_function,
            x_transition,
            r_values,
            p0=p0,
            maxfev=100000
        )

    except RuntimeError as error:

        raise RuntimeError(
            f"Thermodynamic fit failed for transition "
            f"{transition_number}.\n\n"
            f"Experimental x-values:\n"
            f"{x_transition}\n\n"
            f"Relative phasor distances:\n"
            f"{r_values}\n\n"
            "Possible causes:\n"
            "  - The selected transition does not contain "
            "enough experimental points.\n"
            "  - The experimental data do not follow the "
            "thermodynamic model.\n"
            "  - The initial and final states are too similar.\n"
            "  - The fitting parameters could not be determined "
            "from the available data."
        ) from error

    except (ValueError, TypeError) as error:

        raise ValueError(
            f"Invalid numerical values encountered while fitting "
            f"transition {transition_number}.\n"
            "Check the spectra, x-values, integrated fluorescence, "
            "and selected thermodynamic states."
        ) from error

    # ------------------------------------------------------------------
    # Validate fitted parameters
    # ------------------------------------------------------------------

    if not np.all(
        np.isfinite(popt)
    ):
        raise RuntimeError(
            f"Thermodynamic fit for transition "
            f"{transition_number} returned non-finite parameters."
        )

    if pcov.shape != (2, 2):

        raise RuntimeError(
            f"The covariance matrix for transition "
            f"{transition_number} has an unexpected shape."
        )

    with np.errstate(
        invalid="ignore"
    ):
        perr = np.sqrt(
            np.diag(pcov)
        )

    if not np.all(
        np.isfinite(perr)
    ):
        raise RuntimeError(
            f"The parameter uncertainties for transition "
            f"{transition_number} could not be determined."
        )

    dG_fit = popt[0]
    m_fit = popt[1]

    dG_error = perr[0]
    m_error = perr[1]

    # ------------------------------------------------------------------
    # Smooth fitted curve
    #
    # Kept in memory for plotting.
    # It is NOT exported as a separate file.
    # ------------------------------------------------------------------

    x_fit = np.linspace(
        np.min(x_transition),
        np.max(x_transition),
        500
    )

    y_fit = exp_phasor(
        x_fit,
        I_initial,
        I_final,
        dG_fit,
        m_fit,
        temperature,
        gas_constant
    )

    if not np.all(
        np.isfinite(y_fit)
    ):
        raise RuntimeError(
            f"The fitted curve for transition "
            f"{transition_number} contains non-finite values."
        )

    return {
        "transition_number":
            transition_number,

        "initial_index":
            initial_index,

        "final_index":
            final_index,

        "x":
            x_transition,

        "r":
            r_values,

        "x_fit":
            x_fit,

        "y_fit":
            y_fit,

        "I_initial":
            I_initial,

        "I_final":
            I_final,

        "dG":
            dG_fit,

        "dG_error":
            dG_error,

        "m":
            m_fit,

        "m_error":
            m_error,

        "p0":
            p0
    }


# ======================================================================
# FIGURE SAVING
# ======================================================================

def save_figure_without_title(
    figure,
    axis,
    output_file,
    figure_format
):
    """
    Save a figure without its title.

    The current title is temporarily removed before saving and
    restored immediately afterwards.
    """

    current_title = axis.get_title()

    axis.set_title("")

    try:

        if figure_format == "png":

            figure.savefig(
                output_file,
                format="png",
                dpi=300,
                bbox_inches="tight"
            )

        elif figure_format == "tiff":

            figure.savefig(
                output_file,
                format="tiff",
                dpi=300,
                bbox_inches="tight",
                pil_kwargs={
                    "compression": "tiff_lzw"
                }
            )

        else:

            raise ValueError(
                "Unsupported figure format."
            )

    except (OSError, ValueError) as error:

        raise OSError(
            f"Could not save figure:\n"
            f"{output_file}\n\n"
            f"Reason: {error}"
        ) from error

    finally:

        axis.set_title(
            current_title
        )


# ======================================================================
# EXPORT ANALYSIS DATA
# ======================================================================

def export_all_analysis_data(
    output_directory,
    base_name,
    condition,
    labels,
    x_values,
    G_values,
    S_values,
    state_indices,
    transition_results,
    temperature,
    gas_constant,
    energy_unit
):
    """
    Export numerical results from the thermodynamic analysis.

    Generated files:

        1. Human-readable summary text file.
        2. Spectral phasor data text file.
        3. One experimental rd-vs-x text file per transition.

    Smooth fitted curves are intentionally NOT exported.
    """

    output_directory.mkdir(
        parents=True,
        exist_ok=True
    )

    # ------------------------------------------------------------------
    # SUMMARY FILE
    # ------------------------------------------------------------------

    summary_file = (
        output_directory
        /
        f"{base_name}_summary.txt"
    )

    try:

        with open(
            summary_file,
            "w",
            encoding="utf-8"
        ) as file:

            file.write(
                "============================================================\n"
            )
            file.write(
                "SPECTRAL PHASOR THERMODYNAMIC ANALYSIS\n"
            )
            file.write(
                "============================================================\n\n"
            )

            file.write(
                f"Experimental variable:\n{condition}\n\n"
            )

            file.write(
                f"Temperature:\n{temperature:.6g} K\n\n"
            )

            file.write(
                f"Gas constant:\n{gas_constant:.12g} "
                f"{energy_unit}/(mol·K)\n\n"
            )

            file.write(
                f"Energy units:\n{energy_unit}\n\n"
            )

            # ----------------------------------------------------------
            # Loaded spectra
            # ----------------------------------------------------------

            file.write(
                "============================================================\n"
            )
            file.write(
                "SPECTRA LOADED\n"
            )
            file.write(
                "============================================================\n\n"
            )

            file.write(
                f"{'Spectrum':<12}"
                f"{'Label':<20}"
                f"{condition:<20}\n"
            )

            file.write(
                "-" * 52
                +
                "\n"
            )

            for index, (
                label,
                x_value
            ) in enumerate(
                zip(
                    labels,
                    x_values
                ),
                start=1
            ):

                file.write(
                    f"{index:<12}"
                    f"{str(label):<20}"
                    f"{x_value:<20.8g}\n"
                )

            file.write("\n")

            # ----------------------------------------------------------
            # Selected thermodynamic states
            # ----------------------------------------------------------

            file.write(
                "============================================================\n"
            )
            file.write(
                "SELECTED THERMODYNAMIC STATES\n"
            )
            file.write(
                "============================================================\n\n"
            )

            for state_number, index in enumerate(
                state_indices,
                start=1
            ):

                file.write(
                    f"State {state_number}: "
                    f"spectrum {index + 1} "
                    f"({labels[index]}, "
                    f"{condition} = {x_values[index]:.8g})\n"
                )

            file.write("\n")

            # ----------------------------------------------------------
            # Phasor table
            # ----------------------------------------------------------

            file.write(
                "============================================================\n"
            )
            file.write(
                "SPECTRAL PHASOR DATA\n"
            )
            file.write(
                "============================================================\n\n"
            )

            file.write(
                f"{condition:<20}"
                f"{'G':>16}"
                f"{'S':>16}\n"
            )

            file.write(
                "-" * 52
                +
                "\n"
            )

            for x_value, g_value, s_value in zip(
                x_values,
                G_values,
                S_values
            ):

                file.write(
                    f"{x_value:<20.8g}"
                    f"{g_value:>16.8g}"
                    f"{s_value:>16.8g}\n"
                )

            file.write("\n")

            # ----------------------------------------------------------
            # Transition results
            # ----------------------------------------------------------

            file.write(
                "============================================================\n"
            )
            file.write(
                "THERMODYNAMIC FITTING RESULTS\n"
            )
            file.write(
                "============================================================\n\n"
            )

            for result in transition_results:

                transition_number = (
                    result["transition_number"]
                )

                initial_index = (
                    result["initial_index"]
                )

                final_index = (
                    result["final_index"]
                )

                file.write(
                    f"Transition {transition_number}: "
                    f"State {transition_number} ↔ "
                    f"State {transition_number + 1}\n"
                )

                file.write(
                    "-" * 52
                    +
                    "\n"
                )

                file.write(
                    f"Initial state: "
                    f"{labels[initial_index]} "
                    f"({condition} = "
                    f"{x_values[initial_index]:.8g})\n"
                )

                file.write(
                    f"Final state: "
                    f"{labels[final_index]} "
                    f"({condition} = "
                    f"{x_values[final_index]:.8g})\n"
                )

                file.write(
                    f"I_initial: "
                    f"{result['I_initial']:.8g}\n"
                )

                file.write(
                    f"I_final: "
                    f"{result['I_final']:.8g}\n"
                )

                file.write(
                    f"Initial guess dG: "
                    f"{result['p0'][0]:.8g} "
                    f"{energy_unit}\n"
                )

                file.write(
                    f"Initial guess m: "
                    f"{result['p0'][1]:.8g} "
                    f"{energy_unit}/x\n"
                )

                file.write(
                    f"ΔG: "
                    f"{result['dG']:.8g} ± "
                    f"{result['dG_error']:.4g} "
                    f"{energy_unit}\n"
                )

                file.write(
                    f"m: "
                    f"{result['m']:.8g} ± "
                    f"{result['m_error']:.4g} "
                    f"{energy_unit}/x\n"
                )

                file.write("\n")

            # ----------------------------------------------------------
            # Final compact results table
            # ----------------------------------------------------------

            file.write(
                "============================================================\n"
            )
            file.write(
                "FINAL RESULTS TABLE\n"
            )
            file.write(
                "============================================================\n\n"
            )

            headers = ["Parameter"]

            for result in transition_results:

                headers.append(
                    f"State "
                    f"{result['transition_number']} ↔ "
                    f"State "
                    f"{result['transition_number'] + 1}"
                )

            rows = []

            dG_row = ["ΔG"]

            for result in transition_results:

                dG_row.append(
                    f"{result['dG']:.6g} ± "
                    f"{result['dG_error']:.3g} "
                    f"{energy_unit}"
                )

            rows.append(dG_row)

            m_row = ["m"]

            for result in transition_results:

                m_row.append(
                    f"{result['m']:.6g} ± "
                    f"{result['m_error']:.3g} "
                    f"{energy_unit}/x"
                )

            rows.append(m_row)

            column_widths = []

            for column in zip(
                headers,
                *rows
            ):

                width = max(
                    len(str(value))
                    for value in column
                )

                column_widths.append(
                    width + 4
                )

            def format_row(values):
                return " ".join(
                    str(value).center(width)
                    for value, width in zip(
                        values,
                        column_widths
                    )
                )

            file.write(
                format_row(headers)
                +
                "\n"
            )

            file.write(
                "-"
                *
                (
                    sum(column_widths)
                    +
                    len(column_widths)
                    -
                    1
                )
                +
                "\n"
            )

            for row in rows:

                file.write(
                    format_row(row)
                    +
                    "\n"
                )

    except OSError as error:

        raise OSError(
            f"Could not write the summary file:\n"
            f"{summary_file}\n\n"
            f"Reason: {error}"
        ) from error

    # ------------------------------------------------------------------
    # PHASOR DATA FILE
    # ------------------------------------------------------------------

    phasor_file = (
        output_directory
        /
        f"{base_name}_phasor.txt"
    )

    try:

        with open(
            phasor_file,
            "w",
            encoding="utf-8"
        ) as file:

            file.write(
                f"{condition}\tG\tS\n"
            )

            for x_value, g_value, s_value in zip(
                x_values,
                G_values,
                S_values
            ):

                file.write(
                    f"{x_value:.10g}\t"
                    f"{g_value:.10g}\t"
                    f"{s_value:.10g}\n"
                )

    except OSError as error:

        raise OSError(
            f"Could not write the phasor data file:\n"
            f"{phasor_file}\n\n"
            f"Reason: {error}"
        ) from error

    # ------------------------------------------------------------------
    # TRANSITION DATA FILES
    # ------------------------------------------------------------------

    transition_files = []

    for result in transition_results:

        transition_number = (
            result["transition_number"]
        )

        rd_file = (
            output_directory
            /
            f"{base_name}_transition_"
            f"{transition_number}_rd_vs_x.txt"
        )

        try:

            with open(
                rd_file,
                "w",
                encoding="utf-8"
            ) as file:

                file.write(
                    f"{condition}\trd\n"
                )

                for x_value, r_value in zip(
                    result["x"],
                    result["r"]
                ):

                    file.write(
                        f"{x_value:.10g}\t"
                        f"{r_value:.10g}\n"
                    )

        except OSError as error:

            raise OSError(
                f"Could not write the transition data file:\n"
                f"{rd_file}\n\n"
                f"Reason: {error}"
            ) from error

        transition_files.append(
            rd_file
        )

    # ------------------------------------------------------------------
    # Print exported files
    # ------------------------------------------------------------------

    print(
        "\nAll analysis data were exported to:"
    )

    print(
        output_directory
    )

    print(
        f"\nSummary:"
        f"\n  {summary_file.name}"
    )

    print(
        f"\nPhasor data:"
        f"\n  {phasor_file.name}"
    )

    print(
        "\nTransition data:"
    )

    for file_path in transition_files:

        print(
            f"  {file_path.name}"
        )


# ======================================================================
# PRINT FINAL RESULTS TABLE
# ======================================================================

def print_results_table(
    transition_results,
    energy_unit
):
    """Print the final thermodynamic fitting results."""

    print(
        "\n\n=============================================="
    )

    print(
        " Thermodynamic fitting results"
    )

    print(
        "==============================================\n"
    )

    headers = [
        "Parameter"
    ]

    for result in transition_results:

        headers.append(
            f"State "
            f"{result['transition_number']} ↔ "
            f"State "
            f"{result['transition_number'] + 1}"
        )

    rows = []

    dG_row = [
        "ΔG"
    ]

    for result in transition_results:

        dG_row.append(
            f"{result['dG']:.6g} ± "
            f"{result['dG_error']:.3g} "
            f"{energy_unit}"
        )

    rows.append(
        dG_row
    )

    m_row = [
        "m"
    ]

    for result in transition_results:

        m_row.append(
            f"{result['m']:.6g} ± "
            f"{result['m_error']:.3g} "
            f"{energy_unit}/x"
        )

    rows.append(
        m_row
    )

    column_widths = []

    for column in zip(
        headers,
        *rows
    ):

        width = max(
            len(str(value))
            for value in column
        )

        column_widths.append(
            width + 4
        )

    def format_row(values):

        return " ".join(
            str(value).center(width)
            for value, width in zip(
                values,
                column_widths
            )
        )

    print(
        format_row(headers)
    )

    print(
        "-"
        *
        (
            sum(column_widths)
            +
            len(column_widths)
            -
            1
        )
    )

    for row in rows:

        print(
            format_row(row)
        )

    print()


# ======================================================================
# MAIN PROGRAM
# ======================================================================

def main():

    print(
        "\n=============================================="
    )

    print(
        " Spectral Phasor Thermodynamic Analysis"
    )

    print(
        "==============================================\n"
    )

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------

    try:

        (
            wavelengths,
            spectra,
            labels,
            x_values,
            condition
        ) = load_data_with_x_values()

    except Exception as error:

        raise RuntimeError(
            "The spectral data could not be loaded.\n"
            "The problem occurred while reading the selected "
            "folder, files, file format, or header."
        ) from error

    # ------------------------------------------------------------------
    # Validate loaded data
    # ------------------------------------------------------------------

    if x_values is None:

        raise ValueError(
            "No experimental x-values were found.\n"
            "Thermodynamic analysis requires numerical "
            "x-values for each spectrum."
        )

    if len(spectra) < 2:

        raise ValueError(
            "At least two spectra are required "
            "for thermodynamic analysis."
        )

    if len(x_values) != len(spectra):

        raise ValueError(
            "The number of experimental x-values does not match "
            "the number of spectra."
        )

    if len(labels) != len(spectra):

        raise ValueError(
            "The number of spectrum labels does not match "
            "the number of spectra."
        )

    if not np.all(
        np.isfinite(
            np.asarray(
                x_values,
                dtype=float
            )
        )
    ):

        raise ValueError(
            "The experimental x-values contain non-finite values."
        )

    # ------------------------------------------------------------------
    # Sort spectra numerically
    # ------------------------------------------------------------------

    (
        wavelengths,
        spectra,
        labels,
        x_values
    ) = sort_spectral_data_by_x(
        wavelengths,
        spectra,
        labels,
        x_values
    )

    # ------------------------------------------------------------------
    # Display loaded spectra
    # ------------------------------------------------------------------

    print(
        "\nSpectra loaded for thermodynamic analysis:"
    )

    print(
        "--------------------------------------------"
    )

    for index, (
        label,
        x_value
    ) in enumerate(
        zip(
            labels,
            x_values
        ),
        start=1
    ):

        print(
            f"{index}. "
            f"{label} "
            f"(x = {x_value})"
        )

    # ------------------------------------------------------------------
    # Temperature and energy units
    # ------------------------------------------------------------------

    (
        temperature,
        gas_constant,
        energy_unit
    ) = select_temperature_and_units()

    # ------------------------------------------------------------------
    # Calculate spectral phasors
    # ------------------------------------------------------------------

    try:

        G_values = np.array(
            [
                G(
                    wavelengths,
                    spectrum
                )
                for spectrum in spectra
            ],
            dtype=float
        )

        S_values = np.array(
            [
                S(
                    wavelengths,
                    spectrum
                )
                for spectrum in spectra
            ],
            dtype=float
        )

    except (ValueError, TypeError, ZeroDivisionError) as error:

        raise ValueError(
            "Spectral phasor coordinates could not be calculated.\n"
            "Check the wavelength and fluorescence data."
        ) from error

    if not np.all(
        np.isfinite(G_values)
    ):

        raise ValueError(
            "The calculated G phasor coordinates contain "
            "non-finite values.\n"
            "Check the wavelength and fluorescence data."
        )

    if not np.all(
        np.isfinite(S_values)
    ):

        raise ValueError(
            "The calculated S phasor coordinates contain "
            "non-finite values.\n"
            "Check the wavelength and fluorescence data."
        )

    # ------------------------------------------------------------------
    # Select thermodynamic states
    # ------------------------------------------------------------------

    state_indices = (
        select_thermodynamic_states(
            len(spectra)
        )
    )

    print(
        "\nSelected thermodynamic states:"
    )

    print(
        "--------------------------------"
    )

    for state_number, index in enumerate(
        state_indices,
        start=1
    ):

        print(
            f"State {state_number}: "
            f"spectrum {index + 1} "
            f"({labels[index]}, "
            f"x = {x_values[index]})"
        )

    # ------------------------------------------------------------------
    # Spectra figure
    # ------------------------------------------------------------------

    try:

        fig_spectra, ax_spectra = (
            plt.subplots(
                figsize=(8, 6)
            )
        )

        _, colors = create_spectra_plot(
            ax_spectra,
            wavelengths,
            spectra,
            labels=labels
        )

        ax_spectra.set_title(
            "Fluorescence spectra",
            fontsize=18
        )

        fig_spectra.tight_layout()

    except Exception as error:

        raise RuntimeError(
            "The fluorescence spectra figure could not be created."
        ) from error

    # ------------------------------------------------------------------
    # Phasor figure
    # ------------------------------------------------------------------

    try:

        fig_phasor, ax_phasor = (
            plt.subplots(
                figsize=(7, 7)
            )
        )

        create_phasor_plot(
            ax_phasor,
            universal_circle=UNIVERSAL_CIRCLE,
            grid=GRID
        )

        plot_phasor_trajectory(
            ax_phasor,
            G_values,
            S_values,
            colors
        )

        plot_state_connections(
            ax_phasor,
            G_values,
            S_values,
            state_indices
        )

        ax_phasor.set_title(
            "Spectral phasor plot",
            fontsize=18
        )

        fig_phasor.tight_layout()

    except Exception as error:

        raise RuntimeError(
            "The spectral phasor figure could not be created."
        ) from error

    # ------------------------------------------------------------------
    # Thermodynamic analysis
    # ------------------------------------------------------------------

    transition_results = []

    number_of_transitions = (
        len(state_indices) - 1
    )

    for transition_number in range(
        number_of_transitions
    ):

        initial_index = (
            state_indices[
                transition_number
            ]
        )

        final_index = (
            state_indices[
                transition_number + 1
            ]
        )

        print(
            "\n----------------------------------------------"
        )

        print(
            f"Transition "
            f"{transition_number + 1}: "
            f"State "
            f"{transition_number + 1} ↔ "
            f"State "
            f"{transition_number + 2}"
        )

        print(
            "----------------------------------------------"
        )

        result = calculate_transition(
            transition_number + 1,
            initial_index,
            final_index,
            wavelengths,
            spectra,
            G_values,
            S_values,
            x_values,
            temperature,
            gas_constant
        )

        transition_results.append(
            result
        )

        print(
            f"ΔG = "
            f"{result['dG']:.6g} ± "
            f"{result['dG_error']:.3g} "
            f"{energy_unit}"
        )

        print(
            f"m = "
            f"{result['m']:.6g} ± "
            f"{result['m_error']:.3g} "
            f"{energy_unit}/x"
        )

    # ------------------------------------------------------------------
    # Thermodynamic figures
    # ------------------------------------------------------------------

    thermodynamic_figures = []
    thermodynamic_axes = []

    for result in transition_results:

        try:

            fig, ax = plt.subplots(
                figsize=(7, 7)
            )

            # Experimental data

            plot_thermodynamic_graph(
                result["x"],
                result["r"],
                ax,
                color="#0084a5",
                xlabel=condition,
                ylabel="Relative distance",
                tipo="scatter",
                markersize=120
            )

            # Thermodynamic fit

            plot_thermodynamic_graph(
                result["x_fit"],
                result["y_fit"],
                ax,
                color="red",
                xlabel=condition,
                ylabel="Relative distance",
                tipo="line",
                linewidth=2.5
            )

            # Title shown on screen

            ax.set_title(
                f"State "
                f"{result['transition_number']} ↔ "
                f"State "
                f"{result['transition_number'] + 1}",
                fontsize=18
            )

            ax.set_box_aspect(1)

            fig.tight_layout()

        except Exception as error:

            raise RuntimeError(
                f"The thermodynamic figure for transition "
                f"{result['transition_number']} could not be created."
            ) from error

        thermodynamic_figures.append(
            fig
        )

        thermodynamic_axes.append(
            ax
        )

    # ------------------------------------------------------------------
    # Final results table
    # ------------------------------------------------------------------

    print_results_table(
        transition_results,
        energy_unit
    )

    # ------------------------------------------------------------------
    # Display figures
    # ------------------------------------------------------------------

    # Figures are displayed WITH their titles.

    plt.show(
        block=False
    )

    plt.pause(
        0.1
    )

    # ------------------------------------------------------------------
    # Save figures
    # ------------------------------------------------------------------

    save = input(
        "\nSave figures to the 'output' folder? [y/n]: "
    ).strip().lower()

    base_name = None
    output_directory = None

    if save == "y":

        base_name = input(
            "Enter a name for the figures: "
        ).strip()

        if not base_name:

            base_name = (
                "thermodynamic_analysis"
            )

        base_name = Path(
            base_name
        ).stem

        if not base_name:

            raise ValueError(
                "The output file name is empty or invalid."
            )

        print(
            "\nSelect figure format:"
        )

        print(
            "1 - PNG"
        )

        print(
            "2 - TIFF (LZW compressed)"
        )

        figure_format_choice = input(
            "\nChoice: "
        ).strip()

        if figure_format_choice == "1":

            figure_format = "png"

        elif figure_format_choice == "2":

            figure_format = "tiff"

        else:

            raise ValueError(
                "Invalid figure-format choice.\n"
                "Please select 1 or 2."
            )

        output_directory = (
            PROJECT_DIR
            /
            "output"
        )

        try:

            output_directory.mkdir(
                parents=True,
                exist_ok=True
            )

        except OSError as error:

            raise OSError(
                "The 'output' directory could not be created.\n"
                f"Directory: {output_directory}\n"
                f"Reason: {error}"
            ) from error

        # --------------------------------------------------------------
        # Save spectra
        # --------------------------------------------------------------

        spectra_file = (
            output_directory
            /
            f"{base_name}_spectra.{figure_format}"
        )

        save_figure_without_title(
            fig_spectra,
            ax_spectra,
            spectra_file,
            figure_format
        )

        # --------------------------------------------------------------
        # Save phasor
        # --------------------------------------------------------------

        phasor_file = (
            output_directory
            /
            f"{base_name}_phasor.{figure_format}"
        )

        save_figure_without_title(
            fig_phasor,
            ax_phasor,
            phasor_file,
            figure_format
        )

        # --------------------------------------------------------------
        # Save thermodynamic figures
        # --------------------------------------------------------------

        transition_files = []

        for result, fig, ax in zip(
            transition_results,
            thermodynamic_figures,
            thermodynamic_axes
        ):

            transition_number = (
                result["transition_number"]
            )

            transition_file = (
                output_directory
                /
                f"{base_name}_transition_"
                f"{transition_number}."
                f"{figure_format}"
            )

            save_figure_without_title(
                fig,
                ax,
                transition_file,
                figure_format
            )

            transition_files.append(
                transition_file
            )

        print(
            "\nFigures saved in the 'output' folder:"
        )

        print(
            f"  {spectra_file.name}"
        )

        print(
            f"  {phasor_file.name}"
        )

        for file_path in transition_files:

            print(
                f"  {file_path.name}"
            )

    elif save == "n":

        print(
            "\nFigures were not saved."
        )

    else:

        print(
            "\nInvalid option selected. "
            "Figures were not saved."
        )

    # ------------------------------------------------------------------
    # Export analysis data
    # ------------------------------------------------------------------

    export = input(
        "\nExport all analysis data and results? [y/n]: "
    ).strip().lower()

    if export == "y":

        # If figures were not saved, ask for an export name.
        # If figures were saved, reuse the same base name.

        if base_name is None:

            base_name = input(
                "Enter a name for the exported files: "
            ).strip()

            if not base_name:

                base_name = (
                    "thermodynamic_analysis"
                )

            base_name = Path(
                base_name
            ).stem

            if not base_name:

                raise ValueError(
                    "The output file name is empty or invalid."
                )

        if output_directory is None:

            output_directory = (
                PROJECT_DIR
                /
                "output"
            )

            try:

                output_directory.mkdir(
                    parents=True,
                    exist_ok=True
                )

            except OSError as error:

                raise OSError(
                    "The 'output' directory could not be created.\n"
                    f"Directory: {output_directory}\n"
                    f"Reason: {error}"
                ) from error

        export_all_analysis_data(
            output_directory,
            base_name,
            condition,
            labels,
            x_values,
            G_values,
            S_values,
            state_indices,
            transition_results,
            temperature,
            gas_constant,
            energy_unit
        )

    elif export == "n":

        print(
            "\nAnalysis data were not exported."
        )

    else:

        print(
            "\nInvalid option selected. "
            "Analysis data were not exported."
        )

    # ------------------------------------------------------------------
    # Keep figures open
    # ------------------------------------------------------------------

    plt.show()


# ======================================================================
# PROGRAM ENTRY POINT
# ======================================================================

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        print(
            "\n\nAnalysis interrupted by the user."
        )

    except (
        FileNotFoundError,
        ValueError,
        RuntimeError,
        OSError
    ) as error:

        report_error(
            error
        )

    except Exception as error:

        print(
            "\n" + "=" * 70
        )

        print(
            "UNEXPECTED ERROR"
        )

        print(
            "=" * 70
        )

        print(
            f"\nError type: "
            f"{type(error).__name__}"
        )

        print(
            f"Message: {error}"
        )

        print(
            "\nPlease report this message if the problem persists."
        )

        print(
            "=" * 70
        )