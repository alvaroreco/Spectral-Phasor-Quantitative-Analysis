from pathlib import Path
import numpy as np
import pandas as pd


# ============================================================
# BASIC FILE READERS
# ============================================================

def load_txt(file, skiprows=0, delimiter=None):
    """
    Read a text file containing numerical data.

    Parameters
    ----------
    file : str or Path
        Path to the file.
    skiprows : int, optional
        Number of rows to skip at the beginning.
    delimiter : str, optional
        Column delimiter. If None, NumPy attempts to detect
        whitespace-separated columns.

    Returns
    -------
    data : numpy.ndarray
        Numerical data contained in the file.
    """

    file = Path(file)

    return np.loadtxt(
        file,
        skiprows=skiprows,
        delimiter=delimiter
    )


def load_csv(file, skiprows=0):
    """
    Read a CSV file.

    Returns
    -------
    data : pandas.DataFrame
    """

    file = Path(file)

    return pd.read_csv(
        file,
        skiprows=skiprows
    )


def load_excel(file, sheet_name=0):
    """
    Read an Excel file.

    Returns
    -------
    data : pandas.DataFrame
    """

    file = Path(file)

    return pd.read_excel(
        file,
        sheet_name=sheet_name
    )


# ============================================================
# SPECTRUM READERS
# ============================================================

def load_spectrum(file, skiprows=0, delimiter=None):
    """
    Read a file containing one spectrum.

    The first column is assumed to contain wavelength and
    the second column intensity.

    Returns
    -------
    wavelength : numpy.ndarray
    intensity : numpy.ndarray
    """

    data = load_txt(
        file,
        skiprows=skiprows,
        delimiter=delimiter
    )

    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(
            f"File '{file}' does not contain at least two columns."
        )

    wavelength = data[:, 0]
    intensity = data[:, 1]

    return wavelength, intensity


def load_jasco_spectrum(file, skiprows=19):
    """
    Read a JASCO text file containing one fluorescence spectrum.

    By default, the first 19 rows are skipped.
    The first column is wavelength and the second column intensity.

    Returns
    -------
    wavelength : numpy.ndarray
    intensity : numpy.ndarray
    """

    return load_spectrum(
        file,
        skiprows=skiprows
    )


# ============================================================
# MULTIPLE SPECTRA IN ONE FILE
# ============================================================

def load_multispectrum_file(
    file,
    skiprows=0,
    delimiter=None,
    labels=None
):
    """
    Read a file containing several spectra.

    Expected structure:

        Wavelength   Spectrum1   Spectrum2   Spectrum3
        340          1234        1456        1623
        339          1240        1461        1630
        ...

    The first column is wavelength.
    Every subsequent column is interpreted as one spectrum.

    Returns
    -------
    wavelength : numpy.ndarray
    spectra : numpy.ndarray
        Shape = (number_of_spectra, number_of_wavelengths)
    labels : list
        Names assigned to the spectra.
    """

    data = load_txt(
        file,
        skiprows=skiprows,
        delimiter=delimiter
    )

    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(
            f"File '{file}' does not contain multiple spectra."
        )

    wavelength = data[:, 0]
    spectra = data[:, 1:].T

    if labels is None:
        labels = [
            f"Spectrum {i + 1}"
            for i in range(spectra.shape[0])
        ]

    if len(labels) != spectra.shape[0]:
        raise ValueError(
            "The number of labels does not match the number "
            "of spectra."
        )

    return wavelength, spectra, labels


# ============================================================
# FOLDER WITH ONE SPECTRUM PER FILE
# ============================================================

def load_folder(
    folder,
    extension=".txt",
    skiprows=0,
    delimiter=None,
    sort=True
):
    """
    Read all spectrum files contained in a folder.

    Each file is assumed to contain one spectrum:
        column 1 = wavelength
        column 2 = intensity

    Returns
    -------
    wavelength : numpy.ndarray
    spectra : numpy.ndarray
        Shape = (number_of_spectra, number_of_wavelengths)
    labels : list
        File names without extensions.
    """

    folder = Path(folder)

    files = list(folder.glob(f"*{extension}"))

    if len(files) == 0:
        raise FileNotFoundError(
            f"No '{extension}' files were found in '{folder}'."
        )

    if sort:
        files = sorted(files)

    spectra = []
    labels = []
    wavelength = None

    for file in files:

        current_wavelength, intensity = load_spectrum(
            file,
            skiprows=skiprows,
            delimiter=delimiter
        )

        if wavelength is None:
            wavelength = current_wavelength

        else:
            if not np.array_equal(
                wavelength,
                current_wavelength
            ):
                raise ValueError(
                    f"Wavelength axis in '{file}' differs "
                    "from the other files."
                )

        spectra.append(intensity)

        labels.append(file.stem)

    spectra = np.asarray(spectra)

    return wavelength, spectra, labels


def load_jasco_folder(
    folder,
    extension=".txt",
    skiprows=19,
    sort=True
):
    """
    Read a folder containing JASCO spectrum files.

    Each file is assumed to contain one spectrum and the
    first 19 rows are skipped.

    Returns
    -------
    wavelength : numpy.ndarray
    spectra : numpy.ndarray
    labels : list
    """

    return load_folder(
        folder=folder,
        extension=extension,
        skiprows=skiprows,
        sort=sort
    )


# ============================================================
# AUTOMATIC FOLDER READER
# ============================================================

def load_data(
    folder="data",
    extension=".txt",
    skiprows=0,
    delimiter=None
):
    """
    General-purpose reader for a folder containing spectrum
    files.

    Currently, the folder is interpreted as containing one
    spectrum per file.

    Returns
    -------
    wavelength : numpy.ndarray
    spectra : numpy.ndarray
    labels : list
    """

    return load_folder(
        folder=folder,
        extension=extension,
        skiprows=skiprows,
        delimiter=delimiter
    )