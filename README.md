# Spectral-Phasor-Quantitative-Analysis

Python tools for spectral phasor analysis of fluorescence spectra, including visualization, phasor trajectories, and thermodynamic analysis of reversible transitions.

This repository contains the analysis code associated with:

> Recoulat Angelini AA, Malacrida L, Gonzalez Flecha FL. *Linear combination properties of spectral phasors allow thermodynamic analysis of reversible denaturation of multidomain and membrane proteins.*

---

## Overview

Spectral phasors provide a model-independent representation of fluorescence spectra.

For a fluorescence spectrum \(I(\lambda)\), the spectral phasor coordinates are calculated as:

$$
G_{\lambda} =
\sum
\left(
\frac{I(\lambda)}{\sum I(\lambda)}
\right)
\cos
\left[
2\pi
\frac{\lambda-\lambda_{\min}}
{\lambda_{\max}-\lambda_{\min}}
\right]
$$

$$
S_{\lambda} =
\sum
\left(
\frac{I(\lambda)}{\sum I(\lambda)}
\right)
\sin
\left[
2\pi
\frac{\lambda-\lambda_{\min}}
{\lambda_{\max}-\lambda_{\min}}
\right]
$$

The repository provides two main analysis workflows:

* **Spectral phasor analysis and visualization**
* **Thermodynamic analysis of reversible transitions**

---

## Requirements

The repository requires **Python 3.9.x**.

The code was developed and tested with:

* Python 3.9.7
* NumPy 1.20.3
* SciPy 1.7.1
* Matplotlib 3.4.3
* pandas 1.3.4

The exact dependency versions are provided in:

* `requirements.txt`
* `environment.yml`

### Recommended installation

A dedicated Conda environment is recommended.

You do **not** need to replace or downgrade the Python installation on your computer. The Conda environment contains its own Python 3.9 installation and does not modify the system Python.

From the repository directory, run:

```powershell
conda env create -f environment.yml
conda activate spectral-phasor
```

Verify the Python version:

```powershell
python --version
```

It should report Python 3.9.x.

You can also verify that the required packages are available:

```powershell
python -c "import numpy, scipy, matplotlib, pandas; print('Environment OK')"
```

The expected output is:

```text
Environment OK
```

### Alternative installation

If Python 3.9.x is already installed and you are not using Conda, the dependencies can be installed with:

```powershell
python -m pip install -r requirements.txt
```

The Conda environment is recommended because the repository uses specific versions of the scientific Python packages.

---

## Input data

The analysis scripts support two different ways of organizing spectral data:

1. **One spectrum per file**
2. **Multiple spectra in a single file**

The following file formats are supported:

* CSV (`.csv`)
* Text (`.txt`)
* DAT (`.dat`)
* Excel (`.xlsx`)
* Excel (`.xls`)

Input data should be organized inside the `data/` directory. Multiple experimental datasets can be stored in different subfolders.

For example:

```text
data/
├── Example 2 states/
├── Example 3 states/
└── Another dataset/
```

The user selects the appropriate data folder when the program starts.

The selected folder name is also used as the name of the experimental variable in the thermodynamic analysis.

---

## Data organization

### One spectrum per file

When one spectrum is stored per file, the **file name is used as the spectrum label**.

For example, a dataset may be organized as:

```text
data/
└── Dataset/
    ├── spectrum_1.csv
    ├── spectrum_2.csv
    ├── spectrum_3.csv
    └── spectrum_4.csv
```

Each file contains one fluorescence spectrum.

The program uses the file name (without the file extension) as the spectrum identifier.

For example:

```text
spectrum_1
spectrum_2
spectrum_3
spectrum_4
```

For **spectral phasor plotting**, these labels are identifiers only and can be numerical or non-numerical.

For **thermodynamic analysis**, however, the labels must represent numerical values of the experimental variable. The program converts these labels to numerical values and sorts the spectra numerically before performing the analysis.

This prevents alphabetical file ordering from affecting the experimental sequence.

For example, numerical labels such as:

```text
0.000
0.026
0.051
0.100
```

can be used for thermodynamic analysis.

Labels such as:

```text
A1
A2
Sample1
```

cannot be used as thermodynamic experimental values.

---

### Multiple spectra in a single file

When multiple spectra are stored in a single file, the file should contain the spectrum labels in the first row after the selected header rows.

The first column contains the wavelength values.

A generic structure is:

```text
Condition    Label_1    Label_2    Label_3
300.0        0.102      0.115      0.128
301.0        0.104      0.117      0.131
302.0        0.107      0.119      0.134
...
```

Here:

* the first row identifies the individual spectra;
* the first column contains the wavelength values;
* the remaining columns contain the corresponding fluorescence intensities.

For thermodynamic analysis, the spectrum labels must be numerical experimental values.

For example:

```text
Condition    0.000    0.026    0.051    0.100
```

is valid for thermodynamic analysis.

The program uses the selected data folder as the experimental variable name rather than the first cell of the data table.

---

## Interactive data selection

The programs guide the user through the data selection process.

The workflow is:

1. **Select the data folder**
2. **Select the input file format**
3. **Select the data organization**

   * One spectrum per file
   * Multiple spectra in a single file
4. **Select the files to analyze**
5. **Select the number of header rows to skip**

When multiple files are available, the program allows the user to select:

* all files by pressing `Enter`;
* individual files by their index;
* ranges of files;
* combinations of individual files and ranges.

For example:

```text
1-4
```

selects files 1 through 4, while:

```text
1,3,5
```

selects files 1, 3, and 5.

---

## Spectral phasor analysis

The script:

```text
phasor_plot.py
```

performs the following steps:

1. Loads the selected fluorescence spectra.
2. Calculates \(G_{\lambda}\) and \(S_{\lambda}\).
3. Plots the fluorescence spectra.
4. Generates the spectral phasor plot.
5. Displays the phasor trajectory.
6. Optionally exports the phasor coordinates.
7. Optionally saves the generated figures.

The spectrum labels are used only as identifiers in this workflow and therefore do not need to be numerical.

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
r(x)
=
\frac{
d_{\mathrm{initial},x}
}{
d_{\mathrm{initial},\mathrm{final}}
}
$$

where \(d\) represents the Euclidean distance between spectral phasor coordinates.

The relative distance is fitted using:

$$
r(x)
=
\frac{
I_{\mathrm{final}}
}{
I_{\mathrm{final}}
+
I_{\mathrm{initial}}
\exp
\left[
\frac{\Delta G-mx}{RT}
\right]
}
$$

where:

* \(G_{\lambda}\) and \(S_{\lambda}\) are the spectral phasor coordinates;
* \(r(x)\) is the relative spectral phasor distance;
* \(I_{\mathrm{initial}}\) is the integrated fluorescence of the initial state;
* \(I_{\mathrm{final}}\) is the integrated fluorescence of the final state;
* \(\Delta G\) is the standard free-energy difference;
* \(m\) describes the dependence of the free energy on the experimental variable;
* \(R\) is the gas constant;
* \(T\) is the absolute temperature;
* \(x\) is the experimental variable.

The fit provides estimates of:

* \(\Delta G\)
* the standard error of \(\Delta G\)
* \(m\)
* the standard error of \(m\)

---

## Thermodynamic state selection

The initial and final states are selected automatically as the first and last spectra.

Intermediate states can be selected by entering their spectrum numbers.

The selected states must:

* correspond to different spectra;
* be entered in increasing experimental order;
* define valid transitions containing sufficient experimental points.

For a transition between two selected states, the spectra between the initial and final states are used as experimental points for the thermodynamic fit.

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

The program does not require the user to save data or figures. Saving is optional.

### Spectral phasor output

When phasor data are exported, the output contains:

```text
condition    G    S
```

followed by the spectrum labels and their corresponding phasor coordinates.

For example:

```text
0.000    0.1234    0.4567
0.026    0.1251    0.4582
0.051    0.1280    0.4610
```

No redundant numerical index is added to the output.

### Thermodynamic analysis output

The thermodynamic analysis can export:

* a human-readable summary;
* the complete phasor data;
* relative phasor distance data for each transition;
* smooth fitted curves for each transition.

Typical output files are:

```text
<name>_summary.txt
<name>_phasor.txt
<name>_transition_1_rd_vs_x.txt
<name>_transition_1_fit.txt
```

Additional transition files are generated when more than one transition is analyzed.

The summary contains:

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

The programs can optionally save the generated figures in the `output/` directory.

Depending on the analysis, these may include:

* fluorescence spectra;
* spectral phasor plots;
* phasor trajectories;
* thermodynamic transition plots.

Figures are saved as PNG files at 300 dpi.

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
* file-format selection;
* file selection;
* header-row selection;
* CSV, TXT, DAT, XLSX, and XLS input;
* one-spectrum-per-file datasets;
* multiple-spectra-per-file datasets;
* extraction and validation of numerical experimental values.

#### `phasor_functions.py`

Contains the core spectral phasor calculations and plotting functions, including:

* \(G_{\lambda}\)
* \(S_{\lambda}\)
* phasor magnitude;
* phasor angle;
* spectra plotting;
* phasor plotting;
* phasor trajectories;
* state connections;
* figure saving.

#### `phasor_plot.py`

Runs the spectral phasor visualization workflow.

#### `phasor_thermodynamic_analysis.py`

Runs the thermodynamic analysis workflow, including:

* numerical ordering of spectra;
* thermodynamic state selection;
* relative phasor distance calculation;
* thermodynamic fitting;
* parameter estimation;
* data export;
* figure generation.

---

## Reproducibility

For reproducible analysis, use the Conda environment provided with the repository:

```powershell
conda env create -f environment.yml
conda activate spectral-phasor
```

The package versions are explicitly specified in the environment configuration.

---

## Citation

If you use this repository or the associated analysis in your work, please cite:

> Recoulat Angelini AA, Malacrida L, Gonzalez Flecha FL.
> *Linear combination properties of spectral phasors allow thermodynamic analysis of reversible denaturation of multidomain and membrane proteins.*

---

## Contact

**Alvaro A. Recoulat Angelini**
Universidad de Buenos Aires
Facultad de Farmacia y Bioquímica

Email: **[arecoulat@docente.ffyb.uba.ar](mailto:arecoulat@docente.ffyb.uba.ar)**