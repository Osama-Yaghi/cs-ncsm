# CS-NCSM: Complex Scaling Extension for No-Core Shell Model (NCSM)

This repository contains an extension of the **No-Core Shell Model (NCSM)**, designed to calculate the energy levels of light nuclei. The core code (`manyeffv3b.f`) is developed by a collaboration of nuclear theory teams, primarily led by Petr Navratil at TRIUMF, and cannot be shared directly. To obtain access to this base code, please contact **Guillaume Hupin** at [guillaume.hupin@gmail.com](mailto:guillaume.hupin@gmail.com). I share here my contributions. Some of the tools can be run on Example data files that I add without needing the core code.

## 📂 Components

### 1. **Fortran Modules: Complex Scaling Extension**
The core **Fortran modules** in this repository are designed to extend the NCSM framework with **complex scaling**. These modules are integrated into the original **Fortran-based NCSM code** (`manyeffv3b.f`), but **they cannot run as standalone**; they require the base code to function. 

- **n3lo_qp_c.f**: This module modifies the **I-N3LO interaction** with **complex scaled momentum** and softened generators. The parameters of the module need to be **refitted** using neutron-neutron (NN) phase shifts to ensure accurate calculations.

### 2. **Python Tools**

#### **fit.py**: NCSM Input/Output Management & Base Code Execution

This Python script provides an interface for organizing inputs and outputs and running the base code, which can be very complicated and computationally demanding.

- **Functionality**: Simplifies the management of input and output files and helps run the core **Fortran** code, integrating various functions to streamline the workflow.
- **Usage**: The script helps with organizing large-scale calculations and automating tasks in the nuclear structure calculations.
#### **results_analysis/**: Data Analysis & Visualization

This folder contains Python tools for analyzing and visualizing the raw output of the **manyeffv3b** program. It provides functionalities for analyzing eigenvalues and resonance data, helping users to interpret and visualize the results of their NCSM calculations.

**Key Scripts**:
- **spectrum_reader.py**: A tool to read, organize, and visualize eigenvalue data from many-body nuclear structure calculations.
    - **Functions**:
        - `readfiles(nucleus, name1, replace=False)`: Reads eigenvalue data and stores it in a pandas DataFrame.
        - `plot_complex_plane(num_p, first_param, second_param)`: Visualizes eigenvalues in the complex plane.
        - `track_state_rotation(states=[])`: Tracks the rotation of states and determines resonances.
        - And many more to help with data visualization and analysis.
    - For more detailed information, run `python spectrum_reader.py` or `help()` within the script.
    
- **resonance_search.py**: A script used to read and analyze resonance data, comparing it with benchmark results (e.g., R-matrix analysis).
    - **Functions**:
        - `energy_levels(file_names, nucleus, bar_style, plot_thresholds, comparison)`: Compares resonance energy levels with benchmark data.
        - `quick(nuclei1)`: Calculates the critical angle of resonances based on benchmark data.
    - For more detailed information, run `python resonance_search.py` or `help()` within the script.

### 3. **2b_fitting/**: NN Interaction Fitting

This folder contains a Python tool for fitting nuclear-nucleon (NN) interactions in momentum space over a fixed mesh. This tool is required for the **CS-NCSM** calculations.

**Key Features**:
- **Functionality**: Fits the NN interaction for different angular momentum (J) and isospin projections.
- **Input**: Reads data from `fit_input/` (generated via the `generate_interaction_mesh` subroutine in the `fitting_2b` module).
- **Output**: The fitting parameters are written to the `lmfit_parameters/` directory and are used by the parent code in the CS-NCSM model.
- **Parallelization**: The fitting procedure supports parallel execution to handle large datasets efficiently.
- **Example Usage**:
    ```bash
    python -c 'import fitting_lmfit_cost as flc; flc.fit_parallel(0, 4, 40, "np")'
    ```
    This runs a parallel fitting over J values between 0 and 4 with 40 mesh points for each axis and an isospin projection of `'np'`.

For more information, see the `README.md` inside the `2b_fitting/` folder.

---

## 📑 Notes

- **External Code Requirements**: The Fortran modules require the original `manyeffv3b.f` code, which is proprietary to the TRIUMF collaboration. To access the code, contact **Guillaume Hupin** at [guillaume.hupin@gmail.com](mailto:guillaume.hupin@gmail.com).
- **Complex Scaling**: The primary modification in this repository extends the NCSM code with **complex scaling** techniques.

---

## 🛠️ Requirements

- **Fortran 90** for running the core NCSM code.
- **Python 3.x** for the post-processing and analysis scripts.
- **Necessary Python libraries**: `numpy`, `scipy`, `pandas`, `matplotlib` (for plotting), `multiprocessing` (for parallel execution).

---

## 📄 License

 see the LICENSE file for details.

---

## 📧 Contact

For any questions or further information about the repository or collaborations, please contact **Osama Yaghi**.  If you're looking to access the original **Fortran** code, please contact **Guillaume Hupin** at [guillaume.hupin@gmail.com](mailto:guillaume.hupin@gmail.com).

