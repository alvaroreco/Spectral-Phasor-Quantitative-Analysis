from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from spectral_phasor import (
    G,
    S,
    theta,
    M,
    create_phasor_plot
)


# ============================================================
# DATA
# ============================================================

DATA_DIR = Path(__file__).parent / "data"

files = list(DATA_DIR.glob("*"))

print("Files found in data/:")
for file in files:
    print(f"  {file.name}")

def load_spectrum(file):

    extension = file.suffix.lower()

    if extension == ".csv":

        data = pd.read_csv(file)

        wavelength = data.iloc[:, 0].to_numpy()
        intensity = data.iloc[:, 1].to_numpy()

    elif extension == ".txt":

        data = np.loadtxt(file)

        wavelength = data[:, 0]
        intensity = data[:, 1]

    elif extension == ".xlsx":

        data = pd.read_excel(file)

        wavelength = data.iloc[:, 0].to_numpy()
        intensity = data.iloc[:, 1].to_numpy()

    else:

        raise ValueError(
            f"Unsupported file format: {extension}"
        )

    return wavelength, intensity

for file in files:

    wavelength, intensity = load_spectrum(file)

    G_value = G(wavelength, intensity)
    S_value = S(wavelength, intensity)

    print(file.name)
    print(f"Gλ = {G_value:.4f}")
    print(f"Sλ = {S_value:.4f}")