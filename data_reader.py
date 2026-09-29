"""
Functions for loading fluorescence spectral data.

The module allows the user to select the data folder, file format,
data organization, and data files.

For files containing multiple spectra, experimental x-values are
automatically extracted from column labels in the format:

    x=0.000
    x=0.150
    x=0.300
    ...

For files containing one spectrum each, all files are loaded by
default. A subset of files can also be selected manually.

The original load_data() function keeps its three-value return
structure for compatibility with phasor_plot.py.
"""


from pathlib import Path

import numpy as np
import pandas as pd


# -------------------------------------------------------------------
# Select data folder
# -------------------------------------------------------------------

def select_data_folder(data_directory):
    """
    Select the folder containing the experimental data.

    Press Enter to use the main data directory.
    """

    folders = sorted(
        [
            folder
            for folder in data_directory.iterdir()
            if folder.is_dir()
        ]
    )

    if not folders:
        return data_directory

    print("\nAvailable data folders:")

    for i, folder in enumerate(folders, start=1):
        print(f"{i}. {folder.name}")

    while True:

        selection = input(
            "\nSelect a data folder by number "
            "(press Enter for the main 'data' folder): "
        ).strip()

        if selection == "":
            return data_directory

        try:
            selection = int(selection)

            if 1 <= selection <= len(folders):
                return folders[selection - 1]

        except ValueError:
            pass

        print(
            "Invalid selection. Please enter a valid number."
        )


# -------------------------------------------------------------------
# Select file format
# -------------------------------------------------------------------

def select_file_format():
    """
    Select the format of the spectral data files.
    """

    print("\nAvailable file formats:")

    print("1. CSV")
    print("2. TXT")
    print("3. DAT")
    print("4. Excel (.xlsx)")
    print("5. Excel (.xls)")

    formats = {
        "1": ".csv",
        "2": ".txt",
        "3": ".dat",
        "4": ".xlsx",
        "5": ".xls"
    }

    while True:

        selection = input(
            "\nSelect the file format by number: "
        ).strip()

        if selection in formats:
            return formats[selection]

        print(
            "Invalid selection. Please enter a number "
            "between 1 and 5."
        )


# -------------------------------------------------------------------
# Select one data file
# -------------------------------------------------------------------

def select_data_file(folder, file_extension):
    """
    Select a single data file.
    """

    files = sorted(
        [
            file
            for file in folder.iterdir()
            if file.is_file()
            and file.suffix.lower() == file_extension
        ]
    )

    if not files:

        raise FileNotFoundError(
            f"No {file_extension} files were found "
            f"in '{folder}'."
        )

    print(
        f"\nAvailable {file_extension} files:"
    )

    for i, file in enumerate(files, start=1):
        print(f"{i}. {file.name}")

    while True:

        selection = input(
            "\nSelect the data file by number: "
        ).strip()

        try:

            selection = int(selection)

            if 1 <= selection <= len(files):
                return files[selection - 1]

        except ValueError:
            pass

        print(
            "Invalid selection. Please enter a valid number."
        )


# -------------------------------------------------------------------
# Select multiple data files
# -------------------------------------------------------------------

def select_data_files(folder, file_extension):
    """
    Select multiple data files.

    By default, pressing Enter loads all files.

    The user can also specify files using ranges or individual
    indices, for example:

        1-4
        1,3,5
        1-3,5
    """

    files = sorted(
        [
            file
            for file in folder.iterdir()
            if file.is_file()
            and file.suffix.lower() == file_extension
        ]
    )

    if not files:

        raise FileNotFoundError(
            f"No {file_extension} files were found "
            f"in '{folder}'."
        )

    print(
        f"\nAvailable {file_extension} files:"
    )

    for i, file in enumerate(files, start=1):
        print(f"{i}. {file.name}")

    while True:

        selection = input(
            "\nSelect the spectra to load "
            "(press Enter to load all files, "
            "or enter e.g. 1-4 or 1,3,5): "
        ).strip()

        # -----------------------------------------------------------
        # Load all files
        # -----------------------------------------------------------

        if selection == "":
            return files

        # -----------------------------------------------------------
        # Parse selected files
        # -----------------------------------------------------------

        try:

            selected_indices = set()

            parts = selection.split(",")

            for part in parts:

                part = part.strip()

                if "-" in part:

                    start, end = part.split("-", 1)

                    start = int(start.strip())
                    end = int(end.strip())

                    if start > end:
                        raise ValueError

                    selected_indices.update(
                        range(start, end + 1)
                    )

                else:

                    selected_indices.add(
                        int(part)
                    )

            # -------------------------------------------------------
            # Validate indices
            # -------------------------------------------------------

            if not selected_indices:
                raise ValueError

            if any(
                index < 1 or index > len(files)
                for index in selected_indices
            ):
                raise ValueError

            selected_indices = sorted(
                selected_indices
            )

            return [
                files[index - 1]
                for index in selected_indices
            ]

        except ValueError:

            print(
                "Invalid selection. Please use a format such as "
                "'1-4' or '1,3,5'."
            )


# -------------------------------------------------------------------
# Select data organization
# -------------------------------------------------------------------

def select_data_organization():
    """
    Select how the spectra are organized.
    """

    print("\nHow are the spectra organized?")

    print("1. One spectrum per file")
    print("2. Multiple spectra in one file")

    while True:

        selection = input(
            "\nSelect an option: "
        ).strip()

        if selection in ["1", "2"]:
            return int(selection)

        print(
            "Invalid selection. Please enter 1 or 2."
        )


# -------------------------------------------------------------------
# Number of header rows
# -------------------------------------------------------------------

def select_header_rows():
    """
    Ask how many rows should be skipped before reading the data.
    """

    while True:

        selection = input(
            "\nNumber of header rows to skip "
            "(press Enter for 0): "
        ).strip()

        if selection == "":
            return 0

        try:

            skiprows = int(selection)

            if skiprows >= 0:
                return skiprows

        except ValueError:
            pass

        print(
            "Please enter a non-negative integer."
        )


# -------------------------------------------------------------------
# Read spectral file
# -------------------------------------------------------------------

def read_spectral_file(file_path, skiprows=0):
    """
    Read a spectral data file.
    """

    suffix = file_path.suffix.lower()

    if suffix == ".csv":

        data = pd.read_csv(
            file_path,
            skiprows=skiprows
        )

    elif suffix in [".xlsx", ".xls"]:

        data = pd.read_excel(
            file_path,
            skiprows=skiprows
        )

    elif suffix in [".txt", ".dat"]:

        data = pd.read_csv(
            file_path,
            sep=None,
            engine="python",
            header=None,
            skiprows=skiprows
        )

    else:

        raise ValueError(
            f"Unsupported file format: {suffix}"
        )

    return data


# -------------------------------------------------------------------
# Extract x-values
# -------------------------------------------------------------------

def extract_x_values(labels):
    """
    Extract numerical x-values from spectrum labels.

    Accepted formats include:

        x=0.000
        x=0.150
        x=0.300

    or simply:

        0.000
        0.150
        0.300
    """

    x_values = []

    for label in labels:

        label = str(label).strip()

        # -----------------------------------------------------------
        # Remove optional 'x=' prefix
        # -----------------------------------------------------------

        if label.startswith("x="):

            value = label.split(
                "=",
                1
            )[1]

        else:

            value = label

        # -----------------------------------------------------------
        # Convert to numerical value
        # -----------------------------------------------------------

        try:

            x_value = float(value)

        except ValueError as error:

            raise ValueError(
                f"Could not extract a numerical x-value "
                f"from the spectrum label '{label}'.\n\n"
                "Accepted formats are, for example:\n"
                "    x=0.150\n"
                "    0.150"
            ) from error

        x_values.append(
            x_value
        )

    return np.asarray(
        x_values,
        dtype=float
    )

# -------------------------------------------------------------------
# Load multiple spectra from one file
# -------------------------------------------------------------------

def load_multiple_spectra_file(
    file_path,
    skiprows
):
    """
    Load multiple fluorescence spectra from one file.

    The first column contains wavelengths.
    Each subsequent column contains one spectrum.

    Spectrum labels are taken from the column names and x-values
    are extracted automatically.
    """

    data = read_spectral_file(
        file_path,
        skiprows
    )

    if data.shape[1] < 3:

        raise ValueError(
            "A multiple-spectrum file must contain at least "
            "one wavelength column and two spectra."
        )

    # Wavelengths

    wavelengths = np.asarray(
        data.iloc[:, 0],
        dtype=float
    )

    # Spectra

    spectra = np.asarray(
        data.iloc[:, 1:],
        dtype=float
    ).T

    # Labels

    labels = [
        str(label)
        for label in data.columns[1:]
    ]

    # x-values

    x_values = extract_x_values(
        labels
    )

    if len(x_values) != len(spectra):

        raise ValueError(
            "The number of x-values does not match "
            "the number of spectra."
        )

    return (
        wavelengths,
        spectra,
        labels,
        x_values
    )


# -------------------------------------------------------------------
# Load one spectrum
# -------------------------------------------------------------------

def load_single_spectrum_file(
    file_path,
    skiprows
):
    """
    Load one fluorescence spectrum.

    The first column contains wavelengths.
    The second column contains fluorescence intensity.
    """

    data = read_spectral_file(
        file_path,
        skiprows
    )

    if data.shape[1] < 2:

        raise ValueError(
            "The file must contain at least "
            "two columns."
        )

    wavelengths = np.asarray(
        data.iloc[:, 0],
        dtype=float
    )

    spectrum = np.asarray(
        data.iloc[:, 1],
        dtype=float
    )

    return (
        wavelengths,
        spectrum
    )


# -------------------------------------------------------------------
# Load multiple individual spectrum files
# -------------------------------------------------------------------

def load_individual_spectra(
    file_paths,
    skiprows
):
    """
    Load spectra stored in separate files.

    All files must contain the same wavelength values.
    """

    wavelengths = None
    spectra = []
    labels = []

    for file_path in file_paths:

        (
            current_wavelengths,
            spectrum
        ) = load_single_spectrum_file(
            file_path,
            skiprows
        )

        # -----------------------------------------------------------
        # Check wavelengths
        # -----------------------------------------------------------

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
                    "the other spectra."
                )

        spectra.append(
            spectrum
        )

        labels.append(
            file_path.stem
        )

    spectra = np.asarray(
        spectra
    )

    # No x-values are assumed for separate files.

    x_values = None

    return (
        wavelengths,
        spectra,
        labels,
        x_values
    )


# -------------------------------------------------------------------
# Main spectral-data loader
# -------------------------------------------------------------------

def load_spectral_data():
    """
    Load fluorescence spectral data.

    Returns
    -------
    wavelengths : numpy.ndarray

    spectra : numpy.ndarray

    labels : list

    x_values : numpy.ndarray or None

    folder : pathlib.Path
    """

    # ---------------------------------------------------------------
    # Project and data directory
    # ---------------------------------------------------------------

    project_directory = (
        Path(__file__).resolve().parent
    )

    data_directory = (
        project_directory / "data"
    )

    if not data_directory.exists():

        raise FileNotFoundError(
            f"The data directory was not found:\n"
            f"{data_directory}"
        )

    # ---------------------------------------------------------------
    # Select folder
    # ---------------------------------------------------------------

    folder = select_data_folder(
        data_directory
    )

    # ---------------------------------------------------------------
    # Select organization
    # ---------------------------------------------------------------

    organization = select_data_organization()

    # ---------------------------------------------------------------
    # Select format
    # ---------------------------------------------------------------

    file_extension = select_file_format()

    # ---------------------------------------------------------------
    # Select files
    # ---------------------------------------------------------------

    if organization == 1:

        file_paths = select_data_files(
            folder,
            file_extension
        )

    else:

        file_path = select_data_file(
            folder,
            file_extension
        )

    # ---------------------------------------------------------------
    # Header rows
    # ---------------------------------------------------------------

    skiprows = select_header_rows()

    # ---------------------------------------------------------------
    # One spectrum per file
    # ---------------------------------------------------------------

    if organization == 1:

        (
            wavelengths,
            spectra,
            labels,
            x_values
        ) = load_individual_spectra(
            file_paths,
            skiprows
        )

    # ---------------------------------------------------------------
    # Multiple spectra in one file
    # ---------------------------------------------------------------

    else:

        (
            wavelengths,
            spectra,
            labels,
            x_values
        ) = load_multiple_spectra_file(
            file_path,
            skiprows
        )

    # ---------------------------------------------------------------
    # Validation
    # ---------------------------------------------------------------

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
        x_values,
        folder
    )


# -------------------------------------------------------------------
# Original loader
# -------------------------------------------------------------------

def load_data():
    """
    Load fluorescence spectral data.

    This function preserves the original three-value return
    structure used by phasor_plot.py.
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


# -------------------------------------------------------------------
# Loader for thermodynamic analysis
# -------------------------------------------------------------------

def load_data_with_x_values():
    """
    Load fluorescence spectral data together with x-values.
    """

    return load_spectral_data()