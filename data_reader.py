"""
Functions for loading fluorescence spectral data.

The module allows the user to select:
    - the data folder,
    - the data organization,
    - the file format,
    - the data files,
    - the number of header rows to skip.

Two data organizations are supported:

1. One spectrum per file
   The filename is used as the spectrum label.
   For example:
       A1.txt
       A2.txt
       Control.txt

2. Multiple spectra in one file
   The first row contains:
       condition | spectrum_label_1 | spectrum_label_2 | ...

   The first column contains wavelengths and the remaining columns
   contain the corresponding spectra.

For thermodynamic analysis, spectrum labels must be numerical.
For example:
       0.000
       0.026
       0.051

The original load_data() function keeps its three-value return
structure for compatibility.
"""

from pathlib import Path
import re

import numpy as np
import pandas as pd


# =============================================================================
# DATA FOLDER SELECTION
# =============================================================================

def select_data_folder(data_directory):
    """Select the folder containing the spectral data."""

    folders = sorted(
        [folder for folder in data_directory.iterdir() if folder.is_dir()],
        key=lambda path: path.name.lower()
    )

    if not folders:
        return data_directory

    print("\nAvailable data folders:")

    for index, folder in enumerate(folders, start=1):
        print(f"{index}. {folder.name}")

    print(
        "\nSelect a data folder by number "
        "(press Enter for the main 'data' folder): ",
        end=""
    )

    selection = input().strip()

    if selection == "":
        return data_directory

    try:
        selection = int(selection)
    except ValueError:
        raise ValueError("Invalid folder selection.")

    if selection < 1 or selection > len(folders):
        raise ValueError("Invalid folder selection.")

    return folders[selection - 1]


# =============================================================================
# FILE FORMAT SELECTION
# =============================================================================

def select_file_format():
    """Select the spectral file format."""

    formats = {
        "1": ".csv",
        "2": ".txt",
        "3": ".dat",
        "4": ".xlsx",
        "5": ".xls",
    }

    print("\nAvailable file formats:")
    print("1. CSV")
    print("2. TXT")
    print("3. DAT")
    print("4. Excel (.xlsx)")
    print("5. Excel (.xls)")

    selection = input("\nSelect the file format by number: ").strip()

    if selection not in formats:
        raise ValueError("Invalid file format selection.")

    return formats[selection]


# =============================================================================
# FILE SELECTION
# =============================================================================

def select_data_file(folder, file_extension):
    """Select one spectral data file."""

    files = sorted(
        [
            file
            for file in folder.iterdir()
            if file.is_file() and file.suffix.lower() == file_extension
        ],
        key=lambda path: path.name.lower()
    )

    if not files:
        raise FileNotFoundError(
            f"No {file_extension} files were found in:\n{folder}"
        )

    print(f"\nAvailable {file_extension} files:")

    for index, file in enumerate(files, start=1):
        print(f"{index}. {file.name}")

    selection = input("\nSelect the file by number: ").strip()

    try:
        selection = int(selection)
    except ValueError:
        raise ValueError("Invalid file selection.")

    if selection < 1 or selection > len(files):
        raise ValueError("Invalid file selection.")

    return files[selection - 1]


def _parse_file_selection(selection, number_of_files):
    """
    Parse a file selection string.

    Accepted formats:
        Enter       -> all files
        1-4         -> files 1 through 4
        1,3,5       -> files 1, 3 and 5
        1-4,7,10-12 -> combinations
    """

    if selection == "":
        return list(range(number_of_files))

    selected_indices = set()

    parts = [part.strip() for part in selection.split(",")]

    for part in parts:

        if "-" in part:
            bounds = [value.strip() for value in part.split("-")]

            if len(bounds) != 2:
                raise ValueError(
                    f"Invalid range '{part}'."
                )

            try:
                start = int(bounds[0])
                end = int(bounds[1])
            except ValueError:
                raise ValueError(
                    f"Invalid range '{part}'."
                )

            if start > end:
                start, end = end, start

            if start < 1 or end > number_of_files:
                raise ValueError(
                    f"File range '{part}' is outside the available files."
                )

            for index in range(start, end + 1):
                selected_indices.add(index - 1)

        else:
            try:
                index = int(part)
            except ValueError:
                raise ValueError(
                    f"Invalid file selection '{part}'."
                )

            if index < 1 or index > number_of_files:
                raise ValueError(
                    f"File number {index} is outside the available files."
                )

            selected_indices.add(index - 1)

    return sorted(selected_indices)


def select_data_files(folder, file_extension):
    """Select multiple spectral data files."""

    files = sorted(
        [
            file
            for file in folder.iterdir()
            if file.is_file() and file.suffix.lower() == file_extension
        ],
        key=lambda path: path.name.lower()
    )

    if not files:
        raise FileNotFoundError(
            f"No {file_extension} files were found in:\n{folder}"
        )

    print(f"\nAvailable {file_extension} files:")

    for index, file in enumerate(files, start=1):
        print(f"{index}. {file.name}")

    selection = input(
        "\nSelect the spectra to load "
        "(press Enter to load all files, or enter e.g. 1-4 or 1,3,5): "
    ).strip()

    selected_indices = _parse_file_selection(
        selection,
        len(files)
    )

    return [files[index] for index in selected_indices]


# =============================================================================
# DATA ORGANIZATION SELECTION
# =============================================================================

def select_data_organization():
    """Select how the spectra are organized."""

    print("\nHow are the spectra organized?")
    print("1. One spectrum per file")
    print("2. Multiple spectra in one file")

    selection = input("\nSelect an option: ").strip()

    if selection not in {"1", "2"}:
        raise ValueError("Invalid data organization selection.")

    return int(selection)


# =============================================================================
# HEADER SELECTION
# =============================================================================

def select_header_rows():
    """Select the number of header rows to skip."""

    selection = input(
        "\nNumber of header rows to skip (press Enter for 0): "
    ).strip()

    if selection == "":
        return 0

    try:
        skiprows = int(selection)
    except ValueError:
        raise ValueError(
            "The number of header rows must be an integer."
        )

    if skiprows < 0:
        raise ValueError(
            "The number of header rows cannot be negative."
        )

    return skiprows


# =============================================================================
# FILE READING
# =============================================================================

def read_spectral_file(file_path, skiprows=0, header=None):
    """
    Read a spectral data file.

    TXT and DAT files are attempted with several common encodings.
    """

    suffix = file_path.suffix.lower()

    if suffix == ".csv":

        return pd.read_csv(
            file_path,
            skiprows=skiprows,
            header=header
        )

    elif suffix in [".xlsx", ".xls"]:

        return pd.read_excel(
            file_path,
            skiprows=skiprows,
            header=header
        )

    elif suffix in [".txt", ".dat"]:

        encodings = [
            "utf-8",
            "cp1252",
            "latin1",
        ]

        last_error = None

        for encoding in encodings:

            try:
                return pd.read_csv(
                    file_path,
                    sep=None,
                    engine="python",
                    header=header,
                    skiprows=skiprows,
                    encoding=encoding
                )

            except UnicodeDecodeError as error:
                last_error = error

        raise last_error

    else:

        raise ValueError(
            f"Unsupported file format: {suffix}"
        )


# =============================================================================
# NUMERICAL LABEL EXTRACTION
# =============================================================================

def extract_x_value(label):
    """
    Extract a numerical x-value from a spectrum label.

    This function is used only for quantitative thermodynamic analysis.

    Examples accepted:
        0.000
        0.026
        0.051
        -0.100
        1e-3
    """

    label = str(label).strip()

    try:
        return float(label)

    except ValueError as error:

        raise ValueError(
            f"Could not extract a numerical x-value from the "
            f"spectrum label '{label}'.\n\n"
            "For thermodynamic analysis, spectrum labels must be "
            "numerical, for example:\n"
            "    0.000\n"
            "    0.026\n"
            "    0.051"
        ) from error


def extract_x_values(labels):
    """Convert a sequence of numerical spectrum labels into x-values."""

    return np.asarray(
        [extract_x_value(label) for label in labels],
        dtype=float
    )


# =============================================================================
# MULTIPLE SPECTRA IN ONE FILE
# =============================================================================

def load_multiple_spectra_file(file_path, skiprows):
    """
    Load a file containing multiple spectra.

    Expected structure:

        condition    label_1    label_2    label_3
        wavelength   intensity  intensity  intensity
        wavelength   intensity  intensity  intensity
        ...

    The first cell of the first row is the experimental condition.
    The remaining cells of the first row are the spectrum labels.
    """

    data = read_spectral_file(
        file_path,
        skiprows=skiprows,
        header=None
    )

    if data.shape[1] < 3:

        raise ValueError(
            "A multiple-spectrum file must contain at least "
            "one wavelength column and two spectra."
        )

    # First cell of the first row:
    # experimental variable name, e.g. [Denaturant]
    condition = str(data.iloc[0, 0]).strip()

    # Remaining cells of first row:
    # spectrum labels, e.g. 0.000, 0.026, 0.051
    labels = [
        str(label).strip()
        for label in data.iloc[0, 1:]
    ]

    # First column, starting from second row:
    # wavelengths
    wavelengths = pd.to_numeric(
        data.iloc[1:, 0],
        errors="coerce"
    ).to_numpy(dtype=float)

    spectra = []

    for column in range(1, data.shape[1]):

        spectrum = pd.to_numeric(
            data.iloc[1:, column],
            errors="coerce"
        ).to_numpy(dtype=float)

        spectra.append(spectrum)

    spectra = np.asarray(
        spectra,
        dtype=float
    )

    # Remove rows with invalid wavelength values.
    valid = np.isfinite(wavelengths)

    wavelengths = wavelengths[valid]
    spectra = spectra[:, valid]

    if len(labels) != len(spectra):

        raise ValueError(
            "The number of spectrum labels does not match "
            "the number of spectra."
        )

    return (
        wavelengths,
        spectra,
        labels,
        condition
    )


# =============================================================================
# SINGLE SPECTRUM FILE
# =============================================================================

def load_single_spectrum_file(file_path, skiprows):
    """
    Load one spectrum from a file.

    Expected structure:

        wavelength    intensity
        wavelength    intensity
        ...
    """

    data = read_spectral_file(
        file_path,
        skiprows=skiprows,
        header=None
    )

    if data.shape[1] < 2:

        raise ValueError(
            "The file must contain at least two columns."
        )

    wavelengths = pd.to_numeric(
        data.iloc[:, 0],
        errors="coerce"
    ).to_numpy(dtype=float)

    spectrum = pd.to_numeric(
        data.iloc[:, 1],
        errors="coerce"
    ).to_numpy(dtype=float)

    valid = (
        np.isfinite(wavelengths)
        & np.isfinite(spectrum)
    )

    return (
        wavelengths[valid],
        spectrum[valid]
    )


# =============================================================================
# MULTIPLE FILES, ONE SPECTRUM PER FILE
# =============================================================================

def load_individual_spectra(file_paths, skiprows):
    """
    Load several spectra stored in separate files.

    The filename stem is used as the spectrum label.

    Examples:
        A1.txt  -> A1
        A2.txt  -> A2
        Control.txt -> Control

    No numerical x-value is required at this stage.
    """

    wavelengths = None
    spectra = []
    labels = []

    # Alphabetical sorting.
    # The filename is treated only as a label.
    file_paths = sorted(
        file_paths,
        key=lambda path: path.stem.lower()
    )

    condition = file_paths[0].parent.name

    for file_path in file_paths:

        current_wavelengths, spectrum = load_single_spectrum_file(
            file_path,
            skiprows
        )

        if wavelengths is None:

            wavelengths = current_wavelengths

        else:

            if not np.allclose(
                wavelengths,
                current_wavelengths
            ):

                raise ValueError(
                    f"The wavelength values in "
                    f"'{file_path.name}' do not match "
                    f"the other spectra."
                )

        spectra.append(spectrum)

        # Filename stem = spectrum label
        labels.append(file_path.stem)

    spectra = np.asarray(
        spectra,
        dtype=float
    )

    return (
        wavelengths,
        spectra,
        labels,
        condition
    )


# =============================================================================
# MAIN DATA LOADING FUNCTION
# =============================================================================

def load_spectral_data():
    """
    Load spectral data according to the user's selections.

    Returns:
        wavelengths
        spectra
        labels
        condition
        folder
    """

    project_directory = Path(__file__).resolve().parent
    data_directory = project_directory / "data"

    if not data_directory.exists():

        raise FileNotFoundError(
            f"The data directory was not found:\n"
            f"{data_directory}"
        )

    # Select data folder
    folder = select_data_folder(
        data_directory
    )

    # Select organization
    organization = select_data_organization()

    # Select file format
    file_extension = select_file_format()

    if organization == 1:

        # Multiple files, one spectrum per file
        file_paths = select_data_files(
            folder,
            file_extension
        )

    else:

        # Multiple spectra in one file
        file_path = select_data_file(
            folder,
            file_extension
        )

    # Select number of header rows
    skiprows = select_header_rows()

    # -------------------------------------------------------------------------
    # One spectrum per file
    # -------------------------------------------------------------------------

    if organization == 1:

        (
            wavelengths,
            spectra,
            labels,
            condition
        ) = load_individual_spectra(
            file_paths,
            skiprows
        )

    # -------------------------------------------------------------------------
    # Multiple spectra in one file
    # -------------------------------------------------------------------------

    else:

        (
            wavelengths,
            spectra,
            labels,
            condition
        ) = load_multiple_spectra_file(
            file_path,
            skiprows
        )

    # -------------------------------------------------------------------------
    # Validation
    # -------------------------------------------------------------------------

    if len(spectra) < 2:

        raise ValueError(
            "At least two spectra are required "
            "for spectral phasor analysis."
        )

    if len(wavelengths) != spectra.shape[1]:

        raise ValueError(
            "The wavelength array and spectra have "
            "incompatible dimensions."
        )

    return (
        wavelengths,
        spectra,
        labels,
        condition,
        folder
    )


# =============================================================================
# PUBLIC LOADING FUNCTIONS
# =============================================================================

def load_data():
    """
    Load spectral data.

    Returns:
        wavelengths
        spectra
        labels
    """

    (
        wavelengths,
        spectra,
        labels,
        _,
        _
    ) = load_spectral_data()

    return (
        wavelengths,
        spectra,
        labels
    )


def load_data_with_condition():
    """
    Load spectral data for spectral phasor plotting.

    Returns:
        wavelengths
        spectra
        labels
        condition
    """

    (
        wavelengths,
        spectra,
        labels,
        condition,
        _
    ) = load_spectral_data()

    return (
        wavelengths,
        spectra,
        labels,
        condition
    )


def load_data_with_x_values():
    """
    Load spectral data for quantitative thermodynamic analysis.

    Spectrum labels must be numerical because they are used as
    experimental x-values.

    The name of the selected data folder is returned as the
    experimental variable name.

    Returns:
        wavelengths
        spectra
        labels
        x_values
        condition
    """

    (
        wavelengths,
        spectra,
        labels,
        _,
        folder
    ) = load_spectral_data()

    x_values = extract_x_values(
        labels
    )

    # The selected folder name is used as the experimental
    # variable name for thermodynamic analysis.
    condition = folder.name

    return (
        wavelengths,
        spectra,
        labels,
        x_values,
        condition
    )