# Spectral-Phasor-Quantitative-Analysis

Python tools for spectral phasor analysis of fluorescence spectra, including visualization, phasor trajectories, and thermodynamic analysis of reversible transitions.

This repository contains the analysis code associated with:

> Recoulat Angelini AA, Malacrida L, Gonzalez Flecha FL. *Linear combination properties of spectral phasors allow thermodynamic analysis of reversible denaturation of multidomain and membrane proteins.*

---

## Overview

Spectral phasors provide a model-independent representation of fluorescence spectra.

For a fluorescence spectrum $F(\lambda)$, the spectral phasor coordinates are calculated as:

$$
G_{\lambda} =
\frac{\sum_{\lambda} F(\lambda)
\cos\left[\frac{2\pi(\lambda-\lambda_0)}{L}\right]}
{\sum_{\lambda} F(\lambda)}
$$

$$
S_{\lambda} =
\frac{\sum_{\lambda} F(\lambda)
\sin\left[\frac{2\pi(\lambda-\lambda_0)}{L}\right]}
{\sum_{\lambda} F(\lambda)}
$$

where:

* $F(\lambda)$ is the fluorescence intensity or weight associated with each wavelength $\lambda$;
* $\lambda_0$ is the reference wavelength;
* $L$ is the scale parameter that determines the periodicity of the cosine and sine functions;
* the numerator calculates the intensity-weighted sum of the corresponding trigonometric function;
* the denominator normalizes the result by the total fluorescence intensity.

The repository provides two main analysis workflows:

* **Spectral phasor analysis and visualization**
* **Thermodynamic analysis of reversible transitions**

---

## Requirements

The repository was developed and tested with:

* Python 3.9.x
* NumPy 1.20.3
* SciPy 1.7.1
* Matplotlib 3.4.3
* pandas 1.3.4

The recommended installation method is the **Conda environment provided with the repository**.

If Conda is not already installed, **Miniconda is sufficient and is recommended**. The full Anaconda distribution is not required.

Using a dedicated Conda environment avoids conflicts with the Python installation already present on the computer.

The repository does **not** require replacing or modifying the system Python installation.

---

## Installing Conda

The recommended way to run this repository is through a dedicated Conda environment.

If Conda is not already installed on the computer, **Miniconda** is recommended. Miniconda provides Conda and the basic packages required to create the environment, without installing the full Anaconda distribution.

### Installing Miniconda

Download and install Miniconda for your operating system from the official Miniconda documentation:

https://docs.conda.io/projects/miniconda/en/latest/

For Windows, download the **Windows 64-bit** installer unless a different architecture is specifically required.

During installation, the default options are generally sufficient. There is no need to replace or modify an existing Python installation.

After installation, open **Anaconda Prompt** or a new **PowerShell** window and verify that Conda is available:

```powershell
conda --version
```

A version number should be displayed, for example:

```text
conda 25.x.x
```

If `conda` is not recognized in PowerShell immediately after installation, close and reopen the terminal and try again.

Once Conda is available, proceed to [Getting the repository](#getting-the-repository) and then create the dedicated environment using the `environment.yml` file provided with the repository.


## Getting the repository

The programs must be available on the local computer before they can be run.

### Option 1: Clone the repository with Git

If Git is installed, clone the repository using:

```powershell
git clone <repository-url>
````

Then enter the repository directory:

```powershell
cd "Spectral-Phasor-Quantitative-Analysis"
```

Replace `<repository-url>` with the URL of this repository.

### Option 2: Download the repository

Alternatively, download the repository as a ZIP file from GitHub and extract it to a local directory.

After extracting the files, open **Anaconda Prompt** or **PowerShell** with Conda available and navigate to the extracted repository directory:

```powershell
cd "path\to\Spectral-Phasor-Quantitative-Analysis"
```

The repository directory should contain files such as:

```text
README.md
environment.yml
requirements.txt
data_reader.py
phasor_functions.py
phasor_plot.py
phasor_thermodynamic_analysis.py
```

---

## Running the analysis

The recommended way to run the programs is to use the Conda environment provided with the repository.

### Using Conda

Open **Anaconda Prompt** or **PowerShell** with Conda available and navigate to the repository folder:

```powershell
cd "path\to\Spectral-Phasor-Quantitative-Analysis"
```

Create the environment:

```powershell
conda env create -f environment.yml
```

Activate the environment:

```powershell
conda activate spectral-phasor
```

Verify the Python version:

```powershell
python --version
```

The repository was developed and tested with Python 3.9.x.

The required packages can be checked with:

```powershell
python -c "import numpy, scipy, matplotlib, pandas; print('Environment OK')"
```

Then run the desired analysis.

For spectral phasor visualization:

```powershell
python phasor_plot.py
```

For thermodynamic analysis:

```powershell
python phasor_thermodynamic_analysis.py
```

### Using an existing Python 3.9 installation

If Python 3.9.x is already installed and a dedicated environment is not required, the packages can alternatively be installed with:

```powershell
python -m pip install -r requirements.txt
```

Then run:

```powershell
python phasor_plot.py
```

or:

```powershell
python phasor_thermodynamic_analysis.py
```

Using Python versions substantially newer than the version specified by the repository is not recommended because the pinned package versions may not be compatible with newer Python releases.

---

## Data organization

The programs support two explicitly selected data organizations.

The user chooses the organization when the program starts.

### 1. One spectrum per file

Each file contains one fluorescence spectrum.

For example:

```text
data/
└── Example 2 states/
    ├── Native.txt
    ├── Intermediate.txt
    └── Unfolded.txt
```

Each file must contain at least two columns:

```text
wavelength    intensity
300           0.123
301           0.145
302           0.178
...
```

The filename is used as the spectrum label:

```text
Native.txt       → Native
Intermediate.txt → Intermediate
Unfolded.txt     → Unfolded
```

The labels are used as identifiers in the spectral phasor workflow and therefore do not need to be numerical.

For this organization, the **name of the selected data folder** is used as the experimental variable name (`condition`).

For example:

```text
Example 2 states
```

becomes the experimental variable name used for the corresponding output and plot axis.

---

### 2. Multiple spectra in one file

Several fluorescence spectra can be stored in a single file.

The expected structure is:

```text
[Denaturant],0.000,0.150,0.300,0.375,0.450,0.525,0.550,...
300,0.090186038,0.090306682,0.092502228,0.097594762,...
301,0.122232232,0.122379331,0.125055946,0.131259734,...
302,0.162669039,0.162844769,0.166041836,0.173445700,...
...
```

The structure is:

```text
first row:
    first cell      → experimental variable name
    remaining cells → spectrum labels

first column:
    wavelengths

remaining columns:
    fluorescence spectra
```

For the example above:

```text
condition = "[Denaturant]"
```

and the spectrum labels are:

```text
0.000
0.150
0.300
0.375
0.450
...
1.000
```

These numerical labels are used as the experimental $x$ values in the thermodynamic analysis.

Thus:

```text
experimental variable:
    [Denaturant]

x-values:
    0.000
    0.150
    0.300
    0.375
    ...
    1.000
```

The first row is therefore not treated as ordinary spectral data.

For thermodynamic analysis, all spectrum labels in this organization must be numerical.

Examples:

```text
0.000
0.026
0.051
0.100
1e-3
```

Labels such as:

```text
Native
Intermediate
Unfolded
```

are not valid for thermodynamic analysis because they cannot be converted into numerical experimental values.

---

## Data input workflow

When a program is started, the user is asked to select:

1. the data folder;
2. the data organization;
3. the file format;
4. the file or files;
5. the number of header rows to skip.

### Data folder

The program displays the available folders inside:

```text
data/
```

The selected folder contains the spectral files to be analyzed.

### Data organization

The user chooses:

```text
1. One spectrum per file
2. Multiple spectra in one file
```

The two formats are handled separately and are not automatically detected.

### File format

The supported formats are:

```text
CSV
TXT
DAT
Excel (.xlsx)
Excel (.xls)
```

### File selection

For **One spectrum per file**, multiple files can be selected.

Pressing Enter loads all available files.

Specific files or ranges can also be selected, for example:

```text
1-4
```

or:

```text
1,3,5
```

or:

```text
1-4,7,10-12
```

For **Multiple spectra in one file**, one file is selected.

### Header rows

The number of rows to skip before reading the data can be specified.

Pressing Enter selects:

```text
0
```

header rows.

---

## Spectral phasor analysis

The script:

```text
phasor_plot.py
```

performs the following steps:

1. Loads the selected fluorescence spectra.
2. Calculates $G_{\lambda}$ and $S_{\lambda}$.
3. Plots the fluorescence spectra.
4. Generates the spectral phasor plot.
5. Displays the phasor trajectory.
6. Optionally exports the phasor coordinates.
7. Optionally saves the generated figures.

The spectrum labels are used only as identifiers in this workflow and therefore do not need to be numerical.

For example:

```text
A1
A2
Control
Treatment
```

are valid labels for spectral phasor visualization.

---

## Thermodynamic analysis

The script:

```text
phasor_thermodynamic_analysis.py
```

extends the spectral phasor analysis to reversible thermodynamic transitions.

For this analysis, the experimental values associated with the spectra must be numerical.

The spectra are first sorted according to the numerical value of the experimental variable.

The user then selects the thermodynamic states defining the transitions.

For each transition, the program calculates the relative spectral phasor distance:

$$
r(x) =
\frac{
d_{\mathrm{initial},x}
}{
d_{\mathrm{initial},\mathrm{final}}
}
$$

where:

* $d_{\mathrm{initial},x}$ is the Euclidean distance between the phasor coordinates of the initial state and the spectrum measured at experimental condition $x$;
* $d_{\mathrm{initial},\mathrm{final}}$ is the Euclidean distance between the initial and final states.

The relative phasor distance is then fitted using the thermodynamic model:

$$
r(x) =
\frac{
I_{\mathrm{final}}
}{
I_{\mathrm{final}}
+
I_{\mathrm{initial}}
\exp\left[
\frac{\Delta G-mx}{RT}
\right]
}
$$

where:

* $r(x)$ is the relative spectral phasor distance;
* $I_{\mathrm{initial}}$ is the integrated fluorescence of the initial state;
* $I_{\mathrm{final}}$ is the integrated fluorescence of the final state;
* $\Delta G$ is the standard free-energy difference;
* $m$ describes the dependence of the free energy on the experimental variable $x$;
* $R$ is the gas constant;
* $T$ is the absolute temperature;
* $x$ is the experimental variable.

The fit provides estimates of:

* $\Delta G$;
* the standard error of $\Delta G$;
* $m$;
* the standard error of $m$.

---

## Thermodynamic state selection

The initial and final states are selected automatically as the first and last spectra after numerical sorting.

Intermediate states can be selected by entering their spectrum numbers.

The selected states must:

* correspond to different spectra;
* be entered in increasing experimental order;
* define valid transitions containing sufficient experimental points.

For a transition between two selected states, the spectra between the initial and final states are used as experimental points for the thermodynamic fit.

The final state itself is not included as an experimental fitting point for that transition.

---

## Temperature and energy units

The thermodynamic analysis allows the user to specify:

* temperature in Kelvin;
* energy units.

The available energy units are:

* J/mol
* kJ/mol
* cal/mol
* kcal/mol

The corresponding gas constants are handled internally by the program.

---

## Output

All generated files are stored in:

```text
output/
```

The program does not require the user to save data or figures.

Saving is optional.

---

### Spectral phasor output

When phasor data are exported, the output contains:

```text
condition    G    S
```

followed by the spectrum labels and their corresponding phasor coordinates.

For example:

```text
[Denaturant]    G    S
0.000           0.1234    0.4567
0.150           0.1251    0.4582
0.300           0.1280    0.4610
```

No redundant numerical index is added to the output.

The spectrum labels remain identifiers in the spectral phasor workflow.

---

### Thermodynamic analysis output

The thermodynamic analysis can export:

* a human-readable summary;
* the complete phasor data;
* relative phasor distance data for each transition.

The fitted curves are used internally to generate the thermodynamic plots but are **not exported as separate text files**.

Typical output files are:

```text
<name>_summary.txt
<name>_phasor.txt
<name>_transition_1_rd_vs_x.txt
```

Additional transition files are generated when more than one transition is analyzed.

There is no:

```text
<name>_transition_1_fit.txt
```

file because the smooth fitted curves are not exported.

The summary contains information including:

* selected data folder;
* experimental variable;
* temperature;
* energy units;
* loaded spectra;
* selected thermodynamic states;
* spectral phasor coordinates;
* fitted thermodynamic parameters;
* final transition results.

---

## Figures

The programs can optionally save the generated figures in the:

```text
output/
```

directory.

Depending on the analysis, these may include:

* fluorescence spectra;
* spectral phasor plots;
* phasor trajectories;
* thermodynamic transition plots.

Figures are saved at:

```text
300 dpi
```

The user can select one of two formats when saving figures:

```text
1 - PNG
2 - TIFF (LZW compressed)
```

TIFF figures are saved using LZW compression.

### Figure titles

Thermodynamic transition plots are displayed on screen with their corresponding state-transition title.

When the figure is saved to disk, the title is temporarily removed so that the saved figure contains **no title**.

This allows the displayed figure to be identified during analysis while keeping the saved figure suitable for incorporation into manuscripts or other layouts where titles may be added separately.

---

## Error handling

The thermodynamic analysis includes user-oriented error messages for common problems, including:

* incorrect data folder or files;
* incompatible file format;
* incorrect number of header rows;
* non-numerical thermodynamic labels;
* insufficient spectra;
* invalid thermodynamic state selection;
* insufficient experimental points;
* invalid temperature or energy-unit settings;
* unsuccessful thermodynamic fitting;
* output-directory errors.

Technical error information is also displayed to facilitate troubleshooting.

---

## Repository structure

```text
Spectral-Phasor-Quantitative-Analysis/
│
├── README.md
├── requirements.txt
├── environment.yml
│
├── data/
│   ├── Example 2 states/
│   └── Example 3 states/
│
├── output/
│
├── data_reader.py
├── phasor_functions.py
├── phasor_plot.py
└── phasor_thermodynamic_analysis.py
```

### Main modules

#### `data_reader.py`

Handles:

* data-folder selection;
* data-organization selection;
* file-format selection;
* file selection;
* header-row selection;
* CSV, TXT, DAT, XLSX, and XLS input;
* one-spectrum-per-file datasets;
* multiple-spectra-per-file datasets;
* extraction and validation of numerical experimental values;
* determination of the experimental variable name (`condition`).

For one spectrum per file, the experimental variable name is the selected folder name.

For multiple spectra in one file, the experimental variable name is the first cell of the first row of the file.

---

#### `phasor_functions.py`

Contains the core spectral phasor calculations and plotting functions, including:

* $G_{\lambda}$;
* $S_{\lambda}$;
* phasor magnitude;
* phasor angle;
* spectra plotting;
* phasor plotting;
* phasor trajectories;
* state connections.

---

#### `phasor_plot.py`

Runs the spectral phasor visualization workflow.

This includes loading the data, calculating spectral phasors, displaying the spectra and phasor plots, and optionally exporting data and figures.

---

#### `phasor_thermodynamic_analysis.py`

Runs the thermodynamic analysis workflow, including:

* numerical ordering of spectra;
* thermodynamic state selection;
* relative phasor distance calculation;
* thermodynamic fitting;
* parameter estimation;
* data export;
* figure generation;
* optional PNG or LZW-compressed TIFF figure saving.

The thermodynamic fitting curves are retained in memory for visualization but are not exported as separate data files.

---

## Citation

If you use this repository or the associated analysis in your work, please cite:

> Recoulat Angelini AA, Malacrida L, Gonzalez Flecha FL. *Linear combination properties of spectral phasors allow thermodynamic analysis of reversible denaturation of multidomain and membrane proteins.*

---

## Contact

**Alvaro A. Recoulat Angelini**

Universidad de Buenos Aires
Facultad de Farmacia y Bioquímica

Email: **[arecoulat@docente.ffyb.uba.ar](mailto:arecoulat@docente.ffyb.uba.ar)**
