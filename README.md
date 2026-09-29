# Spectral-Phasor-Quantitative-Analysis

Python code for transforming fluorescence emission spectra into spectral phasor coordinates and performing quantitative analysis of reversible protein denaturation.

This repository contains the code used for the analyses presented in:

**Recoulat Angelini AA, Malacrida L, Gonzalez Flecha FL.**
*Linear combination properties of spectral phasors allow thermodynamic analysis of reversible denaturation of multidomain and membrane proteins.*

## Overview

Spectral phasor analysis provides a model-independent representation of fluorescence emission spectra. In this approach, each fluorescence spectrum is transformed into a pair of spectral phasor coordinates, \(G_{\lambda}\) and \(S_{\lambda}\).

The code in this repository allows users to:

* read fluorescence spectral data from CSV files;
* calculate spectral phasor coordinates;
* generate fluorescence spectra and spectral phasor plots;
* visualize phasor trajectories during protein denaturation;
* identify thermodynamic states along a denaturation trajectory;
* calculate relative distances between experimental phasors;
* fit the experimental data to a thermodynamic model of reversible denaturation;
* obtain estimates of \(\Delta G\) and the denaturant dependence parameter \(m\).

The implementation of the spectral phasor calculation follows the mathematical formulation used in the associated publication.

## Requirements

The code was developed and tested using **Python 3.9.7** with the following package versions:

* NumPy 1.20.3
* SciPy 1.7.1
* Matplotlib 3.4.3
* pandas 1.3.4

Two installation options are provided.

### Option 1: Conda environment (recommended)

An `environment.yml` file is provided to create a complete environment with the required Python version, packages, and Tk support for graphical output.

If Anaconda or Miniconda is installed, create the environment from the repository directory:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate spectral-phasor
```

The scripts can then be run directly:

```bash
python phasor_plot.py
```

or:

```bash
python phasor_thermodynamic_analysis.py
```

### Option 2: pip

Alternatively, the required Python packages can be installed using:

```bash
pip install -r requirements.txt
```

The `requirements.txt` file specifies the package versions used for the analyses.

## Input data

Fluorescence spectra should be provided as CSV files.

For data containing multiple spectra, the first column should contain the emission wavelengths and the remaining columns should contain the fluorescence intensity values for each experimental condition.

For example:

```text
Wavelength,x=0.000,x=0.150,x=0.300,x=0.375
300.0,0.0902,0.0903,0.0925,0.0941
301.0,0.1222,0.1224,0.1251,0.1267
302.0,0.1351,0.1360,0.1382,0.1395
...
```

The experimental conditions can be specified in the column headers using the format:

```text
x=0.150
```

The `data_reader.py` module extracts these values automatically.

The input data should be placed in the `data/` directory.

## Spectral phasor analysis

The script

```text
phasor_plot.py
```

reads the spectral data, calculates the spectral phasor coordinates \(G_{\lambda}\) and \(S_{\lambda}\), and generates:

1. fluorescence emission spectra;
2. the corresponding spectral phasor trajectory.

Run the script from the repository directory:

```bash
python phasor_plot.py
```

The script displays the generated figures and provides the option to save them in the `Figures/` directory.

## Thermodynamic analysis

The script

```text
phasor_thermodynamic_analysis.py
```

extends the spectral phasor analysis to the quantitative analysis of reversible protein denaturation.

The workflow includes:

1. loading the fluorescence spectra and experimental denaturant concentrations;
2. calculating spectral phasor coordinates;
3. selecting the thermodynamic states along the phasor trajectory;
4. calculating the relative distance of each experimental phasor between the selected initial and final states;
5. fitting the relative-distance data to a thermodynamic model of reversible denaturation;
6. obtaining estimates of \(\Delta G\) and \(m\) for each transition.

Run the analysis with:

```bash
python phasor_thermodynamic_analysis.py
```

The program interactively requests the temperature, energy units, and thermodynamic states used for the analysis.

## Output

The analysis generates:

* fluorescence emission spectra plots;
* spectral phasor plots;
* thermodynamic fitting plots for each transition;
* a table containing the fitted \(\Delta G\) and \(m\) parameters.

Figures can be saved to the `Figures/` directory from within the scripts.

## Repository structure

```text
Spectral-Phasor-Quantitative-Analysis/
│
├── README.md
├── requirements.txt
├── environment.yml
├── data_reader.py
├── phasor_functions.py
├── phasor_plot.py
├── phasor_thermodynamic_analysis.py
├── data/
└── Figures/
```

### Main modules

#### `data_reader.py`

Functions for reading fluorescence spectral data and extracting experimental conditions from CSV files.

#### `phasor_functions.py`

Functions for calculating spectral phasor coordinates and generating fluorescence and spectral phasor plots.

#### `phasor_plot.py`

Main script for calculating and plotting spectral phasor trajectories from fluorescence emission spectra.

#### `phasor_thermodynamic_analysis.py`

Main script for performing quantitative thermodynamic analysis of reversible denaturation using spectral phasor trajectories.

## Reproducibility

The repository contains the Python implementation used for the analyses reported in the associated publication.

The mathematical implementation of the spectral phasor coordinates and thermodynamic analysis is kept consistent with the analysis described in the publication.

The Python version and package versions used for the analyses are specified in both `requirements.txt` and `environment.yml` to facilitate reproducibility across different systems.

## Citation

If you use this code, please cite the associated publication:

**Recoulat Angelini AA, Malacrida L, Gonzalez Flecha FL.**
*Linear combination properties of spectral phasors allow thermodynamic analysis of reversible denaturation of multidomain and membrane proteins.*
