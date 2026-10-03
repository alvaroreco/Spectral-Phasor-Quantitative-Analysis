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

* $G_{\lambda}$;
* $S_{\lambda}$;
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

## Running the analysis

The recommended way to run the programs is to use the Conda environment provided with the repository. This keeps the required Python version and package versions isolated from the system Python installation.

### Using Conda

Open **Anaconda Prompt** or **PowerShell** with Conda available and navigate to the repository folder:

```powershell
cd path\to\Spectral-Phasor-Quantitative-Analysis
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

Then run the desired analysis:

```powershell
python phasor_plot.py
```

or:

```powershell
python phasor_thermodynamic_analysis.py
```

### Using an existing Python 3.9 installation

If Python 3.9.x is already installed, the required packages can alternatively be installed with:

```powershell
python -m pip install -r requirements.txt
```

Then run the scripts with:

```powershell
python phasor_plot.py
```

or:

```powershell
python phasor_thermodynamic_analysis.py
```

Python 3.12 or newer is not recommended for this repository because the pinned package versions were selected for the Python 3.9 environment used during development.

### PowerShell execution

If PowerShell does not recognize `conda`, initialize Conda for PowerShell from an Anaconda Prompt:

```powershell
conda init powershell
```

Then restart PowerShell and activate the environment:

```powershell
conda activate spectral-phasor
```

The system Python installation does not need to be modified or replaced.


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
