# Spectral-Phasor-Quantitative-Analysis

Python code for transforming fluorescence emission spectra into spectral phasor coordinates and performing quantitative analysis of reversible protein denaturation.

This repository contains the code used for the analyses presented in:

**Recoulat Angelini AA, Malacrida L, Gonzalez Flecha FL.**
*Linear combination properties of spectral phasors allow thermodynamic analysis of reversible denaturation of multidomain and membrane proteins.*

---

## Overview

Spectral phasor analysis provides a model-independent representation of fluorescence emission spectra. In this approach, each fluorescence spectrum is transformed into a pair of spectral phasor coordinates, \(G_{\lambda}\) and \(S_{\lambda}\).

The code in this repository allows users to:

* read fluorescence spectral data from several common file formats;
* calculate spectral phasor coordinates;
* generate fluorescence emission spectra and spectral phasor plots;
* visualize phasor trajectories during protein denaturation;
* select thermodynamic states along a denaturation trajectory;
* calculate relative phasor distances between selected states;
* perform thermodynamic fitting of reversible denaturation transitions;
* obtain estimates of \(\Delta G\) and the denaturant dependence parameter \(m\);
* export numerical results and fitted curves for further analysis.

The implementation of the spectral phasor calculation follows the mathematical formulation used in the associated publication.

---

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

---

## Input data

The repository supports fluorescence spectra organized in two ways:

1. **One spectrum per file**
2. **Multiple spectra in a single file**

The user selects the data organization interactively when running the scripts.

### Supported file formats

The following input formats are supported:

* CSV (`.csv`)
* Text (`.txt`)
* DAT (`.dat`)
* Excel (`.xlsx`)
* Excel (`.xls`)

The input data should be placed inside the `data/` directory.

Multiple subfolders can be used to organize different datasets or experimental conditions.

For example:

```text
data/
├── Example 2 states/
├── Example 3 states/
└── GuCl/
```

The name of the selected data folder is used as the experimental variable name in the thermodynamic analysis.

---

## Multiple spectra in a single file

When several spectra are stored in one file, the first row after the specified header rows should contain the spectrum labels.

The first column should contain the emission wavelengths, and the remaining columns should contain fluorescence intensity values.

For example:

```text
Condition    0.000    0.150    0.300    0.375
300.0        0.0902   0.0903   0.0925   0.0941
301.0        0.1222   0.1224   0.1251   0.1267
302.0        0.1351   0.1360   0.1382   0.1395
...
```

The first cell of the first data row identifies the type of experimental variable, while the remaining cells identify the individual spectra.

For thermodynamic analysis, these spectrum labels **must be numerical**, because they are used as the experimental \(x\)-values.

Examples of valid labels:

```text
0.000
0.026
0.051
0.100
0.250
```

Labels such as:

```text
A1
A2
Sample1
x=0.150
```

are not valid for thermodynamic analysis.

The program will display an informative error message if non-numerical labels are encountered.

---

## One spectrum per file

When one spectrum is stored per file, the file name is used as the spectrum label.

For example:

```text
data/
└── GuCl/
    ├── 0.000.csv
    ├── 0.026.csv
    ├── 0.051.csv
    ├── 0.100.csv
    └── 0.250.csv
```

For thermodynamic analysis, the file names must therefore represent numerical experimental values.

The program sorts the spectra numerically according to these values before performing the thermodynamic analysis.

This prevents alphabetical file ordering from affecting the experimental sequence.

---

## Interactive data selection

The scripts guide the user through the data-loading process.

The workflow includes:

1. selecting the data folder;
2. selecting the input file format;
3. selecting the data organization;
4. selecting individual files or all files;
5. specifying the number of header rows to skip;
6. loading the spectra.

The user can select all files or specific files using file numbers or ranges.

For example:

```text
1-4
```

or:

```text
1,3,5
```

---

# Spectral phasor analysis

The script

```text
phasor_plot.py
```

performs the basic spectral phasor analysis.

Run the script from the repository directory:

```bash
python phasor_plot.py
```

The script:

1. loads the fluorescence spectra;
2. calculates the spectral phasor coordinates \(G_{\lambda}\) and \(S_{\lambda}\);
3. generates the fluorescence emission spectra;
4. generates the spectral phasor plot;
5. displays the phasor trajectory;
6. optionally saves the numerical phasor data;
7. optionally saves the generated figures.

### Spectrum labels

For `phasor_plot.py`, spectrum labels are used only as identifiers.

Therefore, they do **not** need to be numerical.

For example, the following labels are valid:

```text
A1
A2
A3
Control
Sample 1
Sample 2
```

---

# Thermodynamic analysis

The script

```text
phasor_thermodynamic_analysis.py
```

extends the spectral phasor analysis to the quantitative analysis of reversible protein denaturation.

Run the analysis with:

```bash
python phasor_thermodynamic_analysis.py
```

The workflow includes:

1. loading the fluorescence spectra and experimental \(x\)-values;
2. sorting the spectra according to increasing experimental \(x\);
3. calculating spectral phasor coordinates;
4. selecting the thermodynamic states along the phasor trajectory;
5. calculating the relative distance of each experimental phasor between the selected initial and final states;
6. fitting the relative-distance data to a thermodynamic model of reversible denaturation;
7. obtaining estimates of \(\Delta G\) and \(m\);
8. generating thermodynamic fitting plots;
9. optionally exporting all numerical results.

---

## Thermodynamic states

The first and last spectra are automatically selected as thermodynamic states.

The user can optionally select intermediate states by entering their spectrum numbers.

For example:

```text
Intermediate states: 4,7
```

The states must be entered in increasing order.

The program then analyzes each transition between consecutive selected states:

```text
State 1 ↔ State 2
State 2 ↔ State 3
State 3 ↔ State 4
```

This allows multi-state denaturation trajectories to be analyzed as a series of thermodynamic transitions.

---

## Thermodynamic model

The relative phasor distance is fitted using the thermodynamic model:

$$
r(x)=
\frac{I_{\mathrm{final}}}
{I_{\mathrm{final}}+
I_{\mathrm{initial}}
\exp\left[
\frac{\Delta G-mx}{RT}
\right]}
$$

where:

* \(r(x)\) is the relative spectral phasor distance;
* \(I_{\mathrm{initial}}\) is the integrated fluorescence of the initial state;
* \(I_{\mathrm{final}}\) is the integrated fluorescence of the final state;
* \(\Delta G\) is the standard free-energy difference;
* \(m\) describes the dependence of the free energy on the experimental variable;
* \(R\) is the gas constant;
* \(T\) is the absolute temperature;
* \(x\) is the experimental variable.

The initial and final fluorescence values are obtained by numerical integration of the corresponding spectra.

---

## Temperature and energy units

The program interactively requests the temperature in Kelvin.

The available energy units are:

```text
1 - J/mol
2 - kJ/mol
3 - cal/mol
4 - kcal/mol
```

The corresponding gas constants are automatically assigned.

---

## Error handling

The thermodynamic analysis includes checks for common input and fitting problems.

Informative error messages are provided for situations including:

* invalid or missing experimental \(x\)-values;
* non-numerical thermodynamic labels;
* fewer than two spectra;
* inconsistent numbers of spectra, labels, or \(x\)-values;
* invalid temperature;
* invalid energy-unit selection;
* repeated experimental \(x\)-values within a transition;
* insufficient experimental points for a thermodynamic fit;
* identical initial and final phasor coordinates;
* non-finite phasor coordinates;
* invalid integrated fluorescence values;
* failure of the thermodynamic fit to converge;
* problems writing output files.

The program reports the relevant error without changing the underlying mathematical analysis.

---

# Output

All generated files are stored in the:

```text
output/
```

directory.

The spectral phasor analysis can generate:

* fluorescence spectra figures;
* spectral phasor figures;
* numerical \(G\) and \(S\) data.

The thermodynamic analysis can generate:

* fluorescence spectra figures;
* spectral phasor figures;
* one thermodynamic fitting figure for each transition;
* a human-readable summary of the complete analysis;
* a table of spectral phasor coordinates;
* experimental relative-distance data;
* smooth fitted curves;
* fitted \(\Delta G\) and \(m\) values with their estimated standard errors.

---

## Thermodynamic analysis export files

When the user selects the option to export all analysis data, the following files are generated.

### Summary

```text
<name>_summary.txt
```

Contains:

* selected data folder;
* experimental variable;
* temperature;
* gas constant;
* energy units;
* loaded spectra;
* spectrum labels;
* experimental \(x\)-values;
* selected thermodynamic states;
* spectral phasor coordinates;
* integrated initial and final fluorescence;
* initial fitting parameters;
* fitted \(\Delta G\);
* fitted \(m\);
* estimated standard errors;
* final compact results table.

### Phasor data

```text
<name>_phasor.txt
```

Contains:

```text
condition    G    S
x1           ...  ...
x2           ...  ...
x3           ...  ...
```

### Experimental relative-distance data

For each transition:

```text
<name>_transition_1_rd_vs_x.txt
<name>_transition_2_rd_vs_x.txt
...
```

These files contain the experimental \(x\)-values and relative phasor distances.

### Fitted curves

For each transition:

```text
<name>_transition_1_fit.txt
<name>_transition_2_fit.txt
...
```

These files contain the \(x\)-values used to generate the smooth thermodynamic fit and the corresponding fitted relative distances.

---

# Repository structure

```text
Spectral-Phasor-Quantitative-Analysis/
│
├── README.md
├── requirements.txt
├── environment.yml
│
├── data_reader.py
├── phasor_functions.py
├── phasor_plot.py
├── phasor_thermodynamic_analysis.py
│
├── data/
│   ├── Example 2 states/
│   ├── Example 3 states/
│   └── ...
│
└── output/
```

---

## Main modules

### `data_reader.py`

Functions for:

* selecting the data folder;
* selecting the input file format;
* selecting files;
* reading CSV, TXT, DAT, XLSX, and XLS files;
* handling header rows;
* loading one or multiple spectra;
* extracting spectrum labels;
* converting thermodynamic spectrum labels into numerical experimental \(x\)-values.

For thermodynamic analysis, labels are required to be numerical.

---

### `phasor_functions.py`

Functions for:

* calculating spectral phasor coordinates;
* generating fluorescence spectra plots;
* generating spectral phasor plots;
* plotting phasor trajectories;
* plotting connections between selected thermodynamic states;
* saving figures.

---

### `phasor_plot.py`

Main script for basic spectral phasor analysis.

It calculates and displays:

* fluorescence emission spectra;
* spectral phasor coordinates;
* spectral phasor trajectories.

Spectrum labels are treated only as identifiers in this analysis.

---

### `phasor_thermodynamic_analysis.py`

Main script for quantitative thermodynamic analysis using spectral phasor trajectories.

It performs:

* numerical sorting of spectra according to experimental \(x\);
* thermodynamic state selection;
* relative phasor-distance calculations;
* thermodynamic fitting;
* estimation of \(\Delta G\) and \(m\);
* generation of thermodynamic fitting plots;
* export of numerical analysis results.

---

# Reproducibility

The repository contains the Python implementation used for the analyses reported in the associated publication.

The mathematical implementation of the spectral phasor coordinates and thermodynamic analysis is kept consistent with the analysis described in the publication.

The Python version and package versions used for the analyses are specified in both:

```text
requirements.txt
```

and:

```text
environment.yml
```

to facilitate reproducibility across different systems.

---

# Citation

If you use this code in your research, please cite the associated publication:

**Recoulat Angelini AA, Malacrida L, Gonzalez Flecha FL.**

*Linear combination properties of spectral phasors allow thermodynamic analysis of reversible denaturation of multidomain and membrane proteins.*

---

# Contact

For questions regarding the code, implementation, or reproducibility of the analyses, please contact:

**Alvaro A. Recoulat Angelini**
Universidad de Buenos Aires
Facultad de Farmacia y Bioquímica

Email: **[arecoulat@docente.ffyb.uba.ar](mailto:arecoulat@docente.ffyb.uba.ar)**