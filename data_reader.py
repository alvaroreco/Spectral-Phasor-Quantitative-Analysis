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

   The selected data folder name is used as the
   experimental variable name (condition).

2. Multiple spectra in one file
   The first row contains:
       condition | spectrum_label_1 | spectrum_label_2 | ...

   The first cell of the first row is used as the
   experimental variable name.

   The first column contains wavelengths and the remaining
   columns contain the corresponding spectra.

For thermodynamic analysis, spectrum labels must be numerical.
For example:
       0.000
       0.026
       0.051

The original load_data() function keeps its three-value return
structure for compatibility.
"""

from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# DATA FOLDER SELECTION
# =============================================================================

def select_data_folder(data_directory):
    """Select the folder containing the spectral data."""

    try:
        folders = sorted(
            [
                folder
                for folder in data_directory.iterdir()
                if folder.is_dir()
            ],
            key=lambda path: path.name.lower()
        )

    except OSError as error:

        raise OSError(
            f"The data directory could not be accessed:\n"
            f"{data_directory}\n\n"
            "Possible causes:\n"
            "  - The directory does not exist.\n"
            "  - You do not have permission to access it.\n"
            "  - The directory is being used by another program."
        ) from error

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

    except ValueError as error:

        raise ValueError(
            f"Invalid folder selection: '{selection}'.\n\n"
            "Please enter the number of one of the folders "
            "listed above, or press Enter to use the main "
            "'data' folder."
        ) from error

    if selection < 1 or selection > len(folders):

        raise ValueError(
            f"Folder number {selection} is not available.\n\n"
            f"Please select a number between 1 and {len(folders)}."
        )

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

        raise ValueError(
            f"Invalid file format selection: '{selection}'.\n\n"
            "Please select one of the available options:\n"
            "  1 - CSV\n"
            "  2 - TXT\n"
            "  3 - DAT\n"
            "  4 - Excel (.xlsx)\n"
            "  5 - Excel (.xls)"
        )

    return formats[selection]


# =============================================================================
# FILE SELECTION
# =============================================================================

def select_data_file(folder, file_extension):
    """Select one spectral data file."""

    try:
        files = sorted(
            [
                file
                for file in folder.iterdir()
                if (
                    file.is_file()
                    and file.suffix.lower() == file_extension
                )
            ],
            key=lambda path: path.name.lower()
        )

    except OSError as error:

        raise OSError(
            f"The selected data folder could not be accessed:\n"
            f"{folder}\n\n"
            "Possible causes:\n"
            "  - The folder does not exist.\n"
            "  - You do not have permission to access it."
        ) from error

    if not files:

        raise FileNotFoundError(
            f"No {file_extension} files were found in:\n"
            f"{folder}\n\n"
            "Please check that:\n"
            "  - The files are located in the selected folder.\n"
            "  - The selected file format is correct.\n"
            "  - The file extension matches the selected format."
        )

    print(f"\nAvailable {file_extension} files:")

    for index, file in enumerate(files, start=1):
        print(f"{index}. {file.name}")

    selection = input("\nSelect the file by number: ").strip()

    try:
        selection = int(selection)

    except ValueError as error:

        raise ValueError(
            f"Invalid file selection: '{selection}'.\n\n"
            "Please enter the number of one of the files "
            "listed above."
        ) from error

    if selection < 1 or selection > len(files):

        raise ValueError(
            f"File number {selection} is not available.\n\n"
            f"Please select a number between 1 and {len(files)}."
        )

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

    parts = [
        part.strip()
        for part in selection.split(",")
    ]

    for part in parts:

        if not part:

            raise ValueError(
                "An empty file selection was found.\n\n"
                "Please use a format such as:\n"
                "  1-4\n"
                "  1,3,5\n"
                "  1-4,7,10-12"
            )

        if "-" in part:

            bounds = [
                value.strip()
                for value in part.split("-")
            ]

            if len(bounds) != 2:

                raise ValueError(
                    f"Invalid file range '{part}'.\n\n"
                    "Please use a range such as '1-4'."
                )

            try:
                start = int(bounds[0])
                end = int(bounds[1])

            except ValueError as error:

                raise ValueError(
                    f"Invalid file range '{part}'.\n\n"
                    "File ranges must contain two integers, "
                    "for example '1-4'."
                ) from error

            if start > end:
                start, end = end, start

            if start < 1 or end > number_of_files:

                raise ValueError(
                    f"File range '{part}' is outside the "
                    f"available files.\n\n"
                    f"Available file numbers: 1-{number_of_files}."
                )

            for index in range(start, end + 1):
                selected_indices.add(index - 1)

        else:

            try:
                index = int(part)

            except ValueError as error:

                raise ValueError(
                    f"Invalid file selection '{part}'.\n\n"
                    "Please enter file numbers such as "
                    "'1', '1,3,5', or '1-4'."
                ) from error

            if index < 1 or index > number_of_files:

                raise ValueError(
                    f"File number {index} is outside the "
                    f"available files.\n\n"
                    f"Available file numbers: 1-{number_of_files}."
                )

            selected_indices.add(index - 1)

    if not selected_indices:

        raise ValueError(
            "No files were selected.\n\n"
            "Please select at least one file."
        )

    return sorted(selected_indices)


def select_data_files(folder, file_extension):
    """Select multiple spectral data files."""

    try:
        files = sorted(
            [
                file
                for file in folder.iterdir()
                if (
                    file.is_file()
                    and file.suffix.lower() == file_extension
                )
            ],
            key=lambda path: path.name.lower()
        )

    except OSError as error:

        raise OSError(
            f"The selected data folder could not be accessed:\n"
            f"{folder}\n\n"
            "Possible causes:\n"
            "  - The folder does not exist.\n"
            "  - You do not have permission to access it."
        ) from error

    if not files:

        raise FileNotFoundError(
            f"No {file_extension} files were found in:\n"
            f"{folder}\n\n"
            "Please check that:\n"
            "  - The files are located in the selected folder.\n"
            "  - The selected file format is correct.\n"
            "  - The file extension matches the selected format."
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

    return [
        files[index]
        for index in selected_indices
    ]


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

        raise ValueError(
            f"Invalid data organization selection: '{selection}'.\n\n"
            "Please select:\n"
            "  1 - One spectrum per file\n"
            "  2 - Multiple spectra in one file"
        )

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

    except ValueError as error:

        raise ValueError(
            f"Invalid number of header rows: '{selection}'.\n\n"
            "The number of header rows must be a non-negative integer."
        ) from error

    if skiprows < 0:

        raise ValueError(
            "The number of header rows cannot be negative.\n\n"
            "Please enter 0 or a positive integer."
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

    try:

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

    except UnicodeDecodeError as error:

        raise ValueError(
            f"The file could not be decoded:\n"
            f"{file_path}\n\n"
            "The file may use an unsupported text encoding."
        ) from error

    except (
        pd.errors.EmptyDataError,
        pd.errors.ParserError,
        pd.errors.ParserWarning,
    ) as error:

        raise ValueError(
            f"The file could not be parsed:\n"
            f"{file_path}\n\n"
            "Possible causes:\n"
            "  - The file is empty or corrupted.\n"
            "  - The selected number of header rows is incorrect.\n"
            "  - The file delimiter or structure is not supported.\n"
            "  - The selected file format does not match the file contents."
        ) from error

    except (
        FileNotFoundError,
        PermissionError,
        OSError,
    ) as error:

        raise OSError(
            f"The spectral file could not be opened:\n"
            f"{file_path}\n\n"
            "Possible causes:\n"
            "  - The file does not exist.\n"
            "  - The file is open in another program.\n"
            "  - You do not have permission to access the file."
        ) from error


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
        value = float(label)

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

    if not np.isfinite(value):

        raise ValueError(
            f"The spectrum label '{label}' does not represent "
            "a finite numerical x-value.\n\n"
            "Please use numerical labels such as:\n"
            "    0.000\n"
            "    0.026\n"
            "    0.051"
        )

    return value


def extract_x_values(labels):
    """Convert a sequence of numerical spectrum labels into x-values."""

    return np.asarray(
        [
            extract_x_value(label)
            for label in labels
        ],
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

    The first cell of the first row is the experimental
    variable name.

    The remaining cells of the first row are the spectrum labels.
    """

    try:

        data = read_spectral_file(
            file_path,
            skiprows=skiprows,
            header=None
        )

    except Exception as error:

        raise type(error)(
            f"Could not load the multiple-spectrum file:\n"
            f"{file_path}\n\n"
            f"{error}"
        ) from error

    if data.empty:

        raise ValueError(
            f"The selected file is empty:\n"
            f"{file_path}\n\n"
            "Please check the file contents and the number "
            "of header rows to skip."
        )

    if data.shape[1] < 3:

        raise ValueError(
            f"Invalid multiple-spectrum file structure:\n"
            f"{file_path}\n\n"
            "A multiple-spectrum file must contain at least:\n"
            "  - one wavelength column\n"
            "  - two spectrum columns\n"
            "  - one label for each spectrum\n\n"
            "Expected structure:\n"
            "  condition | spectrum 1 | spectrum 2 | ...\n"
            "  wavelength | intensity | intensity | ..."
        )

    # First cell of the first row:
    # experimental variable name, e.g. [Denaturant]
    condition = str(data.iloc[0, 0]).strip()

    if not condition or condition.lower() == "nan":

        raise ValueError(
            f"The experimental variable name is missing in:\n"
            f"{file_path}\n\n"
            "The first cell of the first row must contain "
            "the condition name, for example '[Denaturant]'."
        )

    # Remaining cells of first row:
    # spectrum labels, e.g. 0.000, 0.026, 0.051
    labels = [
        str(label).strip()
        for label in data.iloc[0, 1:]
    ]

    invalid_labels = [
        label
        for label in labels
        if not label or label.lower() == "nan"
    ]

    if invalid_labels:

        raise ValueError(
            f"One or more spectrum labels are missing in:\n"
            f"{file_path}\n\n"
            "Every spectrum must have a label in the first row."
        )

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
    valid_wavelengths = np.isfinite(wavelengths)

    wavelengths = wavelengths[
        valid_wavelengths
    ]

    spectra = spectra[
        :,
        valid_wavelengths
    ]

    if len(wavelengths) == 0:

        raise ValueError(
            f"No valid numerical wavelength values were found in:\n"
            f"{file_path}\n\n"
            "Possible causes:\n"
            "  - The number of header rows is incorrect.\n"
            "  - The wavelength column contains text.\n"
            "  - The selected file does not have the expected format."
        )

    invalid_intensities = ~np.isfinite(spectra)

    if np.any(invalid_intensities):

        invalid_count = np.sum(
            invalid_intensities
        )

        raise ValueError(
            f"Invalid intensity values were found in:\n"
            f"{file_path}\n\n"
            f"Number of non-numerical or non-finite intensity "
            f"values: {invalid_count}\n\n"
            "Possible causes:\n"
            "  - The number of header rows is incorrect.\n"
            "  - One or more intensity values are missing.\n"
            "  - The file contains text in the intensity columns.\n"
            "  - The selected file does not have the expected format."
        )

    if len(labels) != len(spectra):

        raise ValueError(
            f"The number of spectrum labels does not match "
            f"the number of spectra in:\n"
            f"{file_path}\n\n"
            f"Number of labels: {len(labels)}\n"
            f"Number of spectra: {len(spectra)}"
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

    try:

        data = read_spectral_file(
            file_path,
            skiprows=skiprows,
            header=None
        )

    except Exception as error:

        raise type(error)(
            f"Could not load the spectral file:\n"
            f"{file_path}\n\n"
            f"{error}"
        ) from error

    if data.empty:

        raise ValueError(
            f"The spectral file is empty:\n"
            f"{file_path}\n\n"
            "Please check the file contents and the number "
            "of header rows to skip."
        )

    if data.shape[1] < 2:

        raise ValueError(
            f"Invalid spectral file structure:\n"
            f"{file_path}\n\n"
            "The file must contain at least two columns:\n"
            "  1. Wavelength\n"
            "  2. Fluorescence intensity"
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

    if not np.any(valid):

        raise ValueError(
            f"No valid numerical spectral data were found in:\n"
            f"{file_path}\n\n"
            "Possible causes:\n"
            "  - The number of header rows is incorrect.\n"
            "  - The wavelength column contains text.\n"
            "  - The intensity column contains text or missing values.\n"
            "  - The selected file does not have the expected format."
        )

    invalid_wavelengths = np.sum(
        ~np.isfinite(wavelengths)
    )

    invalid_intensities = np.sum(
        ~np.isfinite(spectrum)
    )

    if invalid_wavelengths > 0 or invalid_intensities > 0:

        raise ValueError(
            f"Invalid numerical values were found in:\n"
            f"{file_path}\n\n"
            f"Invalid wavelength values: {invalid_wavelengths}\n"
            f"Invalid intensity values: {invalid_intensities}\n\n"
            "Please check the number of header rows and the "
            "contents of the wavelength and intensity columns."
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

    The parent folder name is used as the experimental
    variable name (condition).

    Examples:
        A1.txt       -> A1
        A2.txt       -> A2
        Control.txt -> Control
    """

    if not file_paths:

        raise ValueError(
            "No spectral files were selected.\n\n"
            "Please select at least one file."
        )

    wavelengths = None
    spectra = []
    labels = []

    # Alphabetical sorting.
    # The filename is treated only as a label.
    file_paths = sorted(
        file_paths,
        key=lambda path: path.stem.lower()
    )

    # For one-spectrum-per-file data, the folder name
    # is the experimental variable name.
    condition = file_paths[0].parent.name

    for file_path in file_paths:

        current_wavelengths, spectrum = (
            load_single_spectrum_file(
                file_path,
                skiprows
            )
        )

        if wavelengths is None:

            wavelengths = current_wavelengths

        else:

            if len(wavelengths) != len(current_wavelengths):

                raise ValueError(
                    f"The wavelength axis in '{file_path.name}' "
                    "has a different number of points from the "
                    "other selected spectra.\n\n"
                    f"Expected: {len(wavelengths)} points\n"
                    f"Found: {len(current_wavelengths)} points"
                )

            if not np.allclose(
                wavelengths,
                current_wavelengths
            ):

                raise ValueError(
                    f"The wavelength values in "
                    f"'{file_path.name}' do not match "
                    "the other selected spectra.\n\n"
                    "All spectra used together must have "
                    "the same wavelength axis."
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
            f"{data_directory}\n\n"
            "Please make sure that the repository contains "
            "a 'data' folder next to this script."
        )

    if not data_directory.is_dir():

        raise NotADirectoryError(
            f"The expected data path is not a directory:\n"
            f"{data_directory}\n\n"
            "Please make sure that 'data' is a folder."
        )

    # Select data folder
    folder = select_data_folder(
        data_directory
    )

    if not folder.exists():

        raise FileNotFoundError(
            f"The selected data folder does not exist:\n"
            f"{folder}"
        )

    if not folder.is_dir():

        raise NotADirectoryError(
            f"The selected data path is not a directory:\n"
            f"{folder}"
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
            "for spectral phasor analysis.\n\n"
            f"Number of spectra loaded: {len(spectra)}\n\n"
            "Please select at least two spectra."
        )

    if len(wavelengths) != spectra.shape[1]:

        raise ValueError(
            "The wavelength array and spectra have "
            "incompatible dimensions.\n\n"
            f"Number of wavelengths: {len(wavelengths)}\n"
            f"Number of intensity points per spectrum: "
            f"{spectra.shape[1]}"
        )

    if not np.all(np.isfinite(wavelengths)):

        raise ValueError(
            "The wavelength axis contains invalid values.\n\n"
            "All wavelength values must be finite numerical values."
        )

    if not np.all(np.isfinite(spectra)):

        raise ValueError(
            "One or more spectra contain invalid intensity values.\n\n"
            "All fluorescence intensity values must be finite "
            "numerical values."
        )

    if len(wavelengths) < 2:

        raise ValueError(
            "The spectral data contain fewer than two valid "
            "wavelength points.\n\n"
            "At least two wavelength points are required "
            "for spectral phasor analysis."
        )

    if not np.all(np.diff(wavelengths) != 0):

        raise ValueError(
            "The wavelength axis contains repeated wavelength values.\n\n"
            "Each wavelength value must be unique."
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

    The experimental variable name is determined by the data
    organization:

        One spectrum per file:
            selected folder name.

        Multiple spectra in one file:
            first cell of the first row of the file.

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
        condition,
        _
    ) = load_spectral_data()

    x_values = extract_x_values(
        labels
    )

    # IMPORTANT:
    # Keep the condition returned by load_spectral_data().
    #
    # One spectrum per file:
    #     condition = selected folder name
    #
    # Multiple spectra in one file:
    #     condition = first row, first column
    #
    # For example:
    #     [Denaturant],0.000,0.150,0.300,...
    #
    # gives:
    #     condition = "[Denaturant]"

    return (
        wavelengths,
        spectra,
        labels,
        x_values,
        condition
    )