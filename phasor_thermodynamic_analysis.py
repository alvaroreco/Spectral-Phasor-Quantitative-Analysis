"""
Spectral phasor analysis for thermodynamic denaturation.

This script loads fluorescence spectra and the associated
experimental variable, calculates spectral phasor coordinates,
allows the user to select thermodynamic states, calculates
relative spectral phasor distances between consecutive states,
and performs thermodynamic fitting of each transition.

The thermodynamic model is:

    r(x) = I_final /
           (I_final + I_initial * exp((dG - m*x)/(R*T)))

where dG is the standard free energy difference and m is the
dependence of the free energy on the experimental variable.
"""

from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import AutoMinorLocator


# ============================================================
# PROJECT PATH
# ============================================================

PROJECT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_DIR))


# ============================================================
# IMPORTS FROM REPOSITORY MODULES
# ============================================================

from data_reader import load_data_with_x_values

from phasor_functions import (
    G,
    S,
    create_spectra_plot,
    create_phasor_plot,
    plot_phasor_trajectory,
    plot_state_connections,
    save_figures
)


# ============================================================
# GENERAL SETTINGS
# ============================================================

UNIVERSAL_CIRCLE = False
GRID = True


# ============================================================
# SELECT THERMODYNAMIC STATES
# ============================================================

def select_thermodynamic_states(number_of_spectra):

    print("\nSelect the thermodynamic states.")
    print("The first and last spectra are automatically selected.")
    print("You can select intermediate states by entering their")
    print("spectrum numbers in increasing order.\n")

    print(f"Number of spectra available: {number_of_spectra}")

    if number_of_spectra < 2:
        raise ValueError(
            "At least two spectra are required for thermodynamic analysis."
        )

    if number_of_spectra == 2:
        return [0, 1]

    print(
        f"\nEnter the intermediate state spectrum numbers "
        f"(2 to {number_of_spectra - 1})."
    )
    print("Separate multiple numbers with commas.")
    print("Press Enter if there are no intermediate states.")

    user_input = input("\nIntermediate states: ").strip()

    if user_input == "":
        return [0, number_of_spectra - 1]

    try:
        intermediate_numbers = [
            int(value.strip())
            for value in user_input.split(",")
        ]
    except ValueError:
        raise ValueError(
            "Invalid state selection. Please enter spectrum numbers "
            "separated by commas."
        )

    if len(set(intermediate_numbers)) != len(intermediate_numbers):
        raise ValueError("Intermediate states cannot be repeated.")

    for number in intermediate_numbers:
        if number < 2 or number > number_of_spectra - 1:
            raise ValueError(
                f"Intermediate state {number} is outside the allowed range."
            )

    if intermediate_numbers != sorted(intermediate_numbers):
        raise ValueError(
            "Intermediate states must be entered in increasing order."
        )

    # Convert from spectrum numbers (1-based) to Python indices (0-based)
    state_indices = [0]

    for number in intermediate_numbers:
        state_indices.append(number - 1)

    state_indices.append(number_of_spectra - 1)

    return state_indices


# ============================================================
# TEMPERATURE AND ENERGY UNITS
# ============================================================

def select_temperature_and_units():

    print("\nThermodynamic parameters")
    print("------------------------")

    temperature_input = input(
        "Temperature (K): "
    ).strip()

    try:
        T = float(temperature_input)
    except ValueError:
        raise ValueError("Temperature must be a numerical value.")

    if T <= 0:
        raise ValueError("Temperature must be greater than zero.")

    print("\nSelect energy units for ΔG:")
    print("1 - J/mol")
    print("2 - kJ/mol")
    print("3 - cal/mol")
    print("4 - kcal/mol")

    unit_choice = input("\nChoice: ").strip()

    if unit_choice == "1":
        R = 8.314462618
        unit = "J/mol"

    elif unit_choice == "2":
        R = 0.008314462618
        unit = "kJ/mol"

    elif unit_choice == "3":
        R = 1.987204258
        unit = "cal/mol"

    elif unit_choice == "4":
        R = 0.001987204258
        unit = "kcal/mol"

    else:
        raise ValueError(
            "Invalid choice. Please select 1, 2, 3, or 4."
        )

    print(f"\nTemperature: {T:.2f} K")
    print(f"Energy units: {unit}")
    print(f"Gas constant: {R}")

    return T, R, unit


# ============================================================
# DISTANCE FUNCTIONS
# ============================================================

def distance(G1, S1, G2, S2):
    """
    Euclidean distance between two spectral phasor coordinates.
    """
    return np.sqrt(
        (G2 - G1) ** 2 +
        (S2 - S1) ** 2
    )


def rd(G1, G2, Gi, S1, S2, Si):
    """
    Relative distance from the initial state to an experimental
    phasor point, normalized by the distance between the initial
    and final thermodynamic states.
    """

    denominator = distance(G1, S1, G2, S2)

    if denominator == 0:
        raise ValueError(
            "The initial and final thermodynamic states have "
            "identical phasor coordinates."
        )

    return (
        distance(G1, S1, Gi, Si)
        / denominator
    )


# ============================================================
# THERMODYNAMIC MODEL
# ============================================================

def exp_phasor(x, I_initial, I_final, dG, m, T, R):
    """
    Thermodynamic model for the relative phasor distance.
    """

    exp_term = np.exp(
        (dG - m * x) / (R * T)
    )

    return I_final / (
        I_final +
        I_initial * exp_term
    )


# ============================================================
# PLOT FORMATTING
# ============================================================

def graph(
    x,
    y,
    ax,
    color='black',
    xlabel='x',
    ylabel='y',
    tipo='line',
    linewidth=2.5,
    markersize=120,
    fontsize=20
):
    """
    Generic plotting function for thermodynamic analysis.

    Parameters
    ----------
    x, y : array-like
        Data to plot.
    ax : matplotlib.axes.Axes
        Axis on which the data are plotted.
    color : str
        Plot color.
    xlabel, ylabel : str
        Axis labels.
    tipo : str
        'line' or 'scatter'.
    linewidth : float
        Line width for line plots.
    markersize : float
        Marker size for scatter plots.
    fontsize : float
        Axis-label font size.
    """

    # --------------------------------------------------------
    # Axis formatting
    # --------------------------------------------------------

    for side in ['top', 'right', 'bottom', 'left']:
        ax.spines[side].set_linewidth(2)

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    if tipo == 'line':

        ax.plot(
            x,
            y,
            color=color,
            linewidth=linewidth
        )

    elif tipo == 'scatter':

        ax.scatter(
            x,
            y,
            color=color,
            s=markersize
        )

    else:
        raise ValueError(
            "tipo must be either 'line' or 'scatter'."
        )

    # --------------------------------------------------------
    # Labels
    # --------------------------------------------------------

    ax.set_xlabel(
        xlabel,
        fontsize=fontsize,
        fontname='sans-serif',
        labelpad=10
    )

    ax.set_ylabel(
        ylabel,
        fontsize=fontsize,
        fontname='sans-serif',
        labelpad=10
    )

    # --------------------------------------------------------
    # Tick formatting
    # --------------------------------------------------------

    ax.tick_params(
        axis='both',
        labelsize=16,
        pad=10
    )

    for label in (
        ax.get_xticklabels()
        + ax.get_yticklabels()
    ):
        label.set_fontfamily('sans-serif')

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
# ============================================================
# ANALYZE ONE TRANSITION
# ============================================================

def calculate_transition(
    transition_number,
    initial_index,
    final_index,
    wavelengths,
    spectra,
    G_values,
    S_values,
    x_values,
    T,
    R
):

    # --------------------------------------------------------
    # Experimental points belonging to this transition
    # --------------------------------------------------------

    indices = np.arange(
        initial_index,
        final_index
    )

    x_transition = np.asarray(
        x_values[indices],
        dtype=float
    )

    # --------------------------------------------------------
    # Initial and final phasor coordinates
    # --------------------------------------------------------

    G_initial = G_values[initial_index]
    S_initial = S_values[initial_index]

    G_final = G_values[final_index]
    S_final = S_values[final_index]

    # --------------------------------------------------------
    # Relative phasor distances
    # --------------------------------------------------------

    r_values = np.array([
        rd(
            G_initial,
            G_final,
            G_values[index],
            S_initial,
            S_final,
            S_values[index]
        )
        for index in indices
    ])

    # --------------------------------------------------------
    # Integrated fluorescence of initial and final spectra
    # --------------------------------------------------------

    I_initial = np.trapz(
        spectra[initial_index],
        wavelengths
    )
    
    I_final = np.trapz(
        spectra[final_index],
        wavelengths
    )

    # --------------------------------------------------------
    # Initial guesses according to energy units
    #
    # Original values:
    # J/mol  -> [10000, 20000]
    # kJ/mol -> [10, 20]
    # cal/mol -> [2390, 4780]
    # kcal/mol -> [2.39, 4.78]
    # --------------------------------------------------------

    if np.isclose(R, 8.314462618):
        p0 = [10000, 20000]

    elif np.isclose(R, 0.008314462618):
        p0 = [10, 20]

    elif np.isclose(R, 1.987204258):
        p0 = [2390, 4780]

    elif np.isclose(R, 0.001987204258):
        p0 = [2.39, 4.78]

    else:
        p0 = [10000, 20000]

   # --------------------------------------------------------
   # Thermodynamic fit
   # --------------------------------------------------------

    def fit_function(x, dG, m):
        return exp_phasor(
            x,
            I_initial,
            I_final,
            dG,
            m,
            T,
            R
        )
    
    
    popt, pcov = curve_fit(
        fit_function,
        x_transition,
        r_values,
        p0=p0,
        maxfev=100000
    )
    
    dG_fit = popt[0]
    m_fit = popt[1]
    
    perr = np.sqrt(
        np.diag(pcov)
    )
    
    dG_error = perr[0]
    m_error = perr[1]

    dG_fit = popt[0]
    m_fit = popt[1]

    perr = np.sqrt(
        np.diag(pcov)
    )

    dG_error = perr[0]
    m_error = perr[1]

    # --------------------------------------------------------
    # Smooth fitted curve
    # --------------------------------------------------------

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
        T,
        R
    )

    return {
        'transition_number': transition_number,
        'initial_index': initial_index,
        'final_index': final_index,
        'x': x_transition,
        'r': r_values,
        'x_fit': x_fit,
        'y_fit': y_fit,
        'I_initial': I_initial,
        'I_final': I_final,
        'dG': dG_fit,
        'dG_error': dG_error,
        'm': m_fit,
        'm_error': m_error,
        'p0': p0
    }


# ============================================================
# MAIN PROGRAM
# ============================================================

print("\n==============================================")
print(" Spectral Phasor Thermodynamic Analysis")
print("==============================================\n")


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

(
    wavelengths,
    spectra,
    labels,
    x_values,
    data_folder
) = load_data_with_x_values()


if x_values is None:
    raise ValueError(
        "\nNo experimental x-values were found.\n"
        "Thermodynamic analysis requires x-values for each spectrum."
    )


# ------------------------------------------------------------
# Temperature and energy units
# ------------------------------------------------------------

T, R, energy_unit = select_temperature_and_units()


# ------------------------------------------------------------
# Calculate spectral phasors
# ------------------------------------------------------------

G_values = np.array([
    G(wavelengths, spectrum)
    for spectrum in spectra
])

S_values = np.array([
    S(wavelengths, spectrum)
    for spectrum in spectra
])


# ------------------------------------------------------------
# Select thermodynamic states
# ------------------------------------------------------------

state_indices = select_thermodynamic_states(
    len(spectra)
)


print("\nSelected thermodynamic states:")
print("--------------------------------")

for state_number, index in enumerate(
    state_indices,
    start=1
):

    print(
        f"State {state_number}: "
        f"spectrum {index + 1} "
        f"(x = {x_values[index]})"
    )


# ============================================================
# SPECTRA FIGURE
# ============================================================

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


# ============================================================
# PHASOR FIGURE
# ============================================================

fig_phasor, ax_phasor = plt.subplots(
    figsize=(7, 7)
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

# Connect selected thermodynamic states
plot_state_connections(
    ax_phasor,
    G_values,
    S_values,
    state_indices
)

fig_phasor.tight_layout()


# ============================================================
# THERMODYNAMIC ANALYSIS
# ============================================================

transition_results = []

number_of_transitions = (
    len(state_indices) - 1
)

for transition_number in range(
    number_of_transitions
):

    initial_index = state_indices[
        transition_number
    ]

    final_index = state_indices[
        transition_number + 1
    ]

    print("\n----------------------------------------------")
    print(
        f"Transition {transition_number + 1}: "
        f"State {transition_number + 1} ↔ "
        f"State {transition_number + 2}"
    )
    print("----------------------------------------------")

    result = calculate_transition(
        transition_number + 1,
        initial_index,
        final_index,
        wavelengths,
        spectra,
        G_values,
        S_values,
        x_values,
        T,
        R
    )

    transition_results.append(result)

    print(
        f"ΔG = {result['dG']:.6g} ± "
        f"{result['dG_error']:.3g} {energy_unit}"
    )

    print(
        f"m = {result['m']:.6g} ± "
        f"{result['m_error']:.3g} {energy_unit}/x"
    )


# ============================================================
# THERMODYNAMIC FIGURES
# ============================================================

thermodynamic_figures = []

for result in transition_results:

    fig, ax = plt.subplots(
        figsize=(7, 7)
    )

    # Experimental data
    graph(
        result['x'],
        result['r'],
        ax,
        color='#0084a5',
        xlabel='x',
        ylabel='Relative distance',
        tipo='scatter',
        markersize=120
    )

    # Thermodynamic fit
    graph(
        result['x_fit'],
        result['y_fit'],
        ax,
        color='red',
        xlabel='x',
        ylabel='Relative distance',
        tipo='line',
        linewidth=2.5
    )

    ax.set_title(
        f"State {result['transition_number']} ↔ "
        f"State {result['transition_number'] + 1}",
        fontsize=18
    )

    # Square aspect ratio
    ax.set_box_aspect(1)

    fig.tight_layout()

    thermodynamic_figures.append(fig)


# ============================================================
# FINAL RESULTS TABLE
# ============================================================

print("\n\n==============================================")
print(" Thermodynamic fitting results")
print("==============================================\n")


# ------------------------------------------------------------
# Create table columns
# ------------------------------------------------------------

headers = ['Parameter']

for result in transition_results:

    headers.append(
        f"State {result['transition_number']} ↔ "
        f"State {result['transition_number'] + 1}"
    )


rows = []


# ΔG row
dG_row = ['ΔG']

for result in transition_results:

    dG_row.append(
        f"{result['dG']:.6g} ± "
        f"{result['dG_error']:.3g} {energy_unit}"
    )

rows.append(dG_row)


# m row
m_row = ['m']

for result in transition_results:

    m_row.append(
        f"{result['m']:.6g} ± "
        f"{result['m_error']:.3g} {energy_unit}/x"
    )

rows.append(m_row)


# ------------------------------------------------------------
# Determine column widths
# ------------------------------------------------------------

column_widths = []

for column in zip(
    headers,
    *rows
):

    width = max(
        len(str(value))
        for value in column
    )

    # Add some horizontal padding
    column_widths.append(
        width + 4
    )


# ------------------------------------------------------------
# Print centered table
# ------------------------------------------------------------

def format_row(values):

    return " ".join(
        str(value).center(width)
        for value, width in zip(
            values,
            column_widths
        )
    )


print(format_row(headers))

print(
    "-" * (
        sum(column_widths)
        + len(column_widths) - 1
    )
)

for row in rows:
    print(format_row(row))

print()


# ============================================================
# SHOW FIGURES
# ============================================================

plt.show(block=False)
plt.pause(0.1)


# ============================================================
# SAVE FIGURES
# ============================================================

save = input(
    "\nSave figures to the 'Figures' folder? [y/n]: "
).strip().lower()


if save == 'y':

    base_name = input(
        "Enter a name for the figures: "
    ).strip()

    if not base_name:
        base_name = 'thermodynamic_analysis'

    # Spectra + phasor
    save_figures(
        fig_spectra,
        fig_phasor,
        PROJECT_DIR / 'Figures',
        base_name=base_name
    )

    # Thermodynamic transitions
    output_directory = (
        PROJECT_DIR / 'Figures'
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True
    )

    for result, fig in zip(
        transition_results,
        thermodynamic_figures
    ):

        transition_number = (
            result['transition_number']
        )

        fig.savefig(
            output_directory /
            f'{base_name}_transition_{transition_number}.png',
            dpi=300,
            bbox_inches='tight'
        )

    print(
        f"\nFigures saved in the 'Figures' folder "
        f"with base name '{base_name}'."
    )


# Keep figures open
plt.show()