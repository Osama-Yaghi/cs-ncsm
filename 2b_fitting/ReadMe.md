# NN Interaction Fitting Script

This Python script is designed to fit nuclear-nucleon (NN) interactions in momentum space using a custom fitting procedure. It operates on a mesh grid and employs various functions to fit the interaction over different values of total angular momentum (J) and isospin projections.

The script is intended to be used from the terminal with simple Python commands, making it easy to integrate into larger workflows. The fitting procedure is highly customizable and can be run in parallel for better performance.

## Overview

- **Functionality**: Fits the NN interaction for various values of total angular momentum (J) and isospin projections.
- **Input**: The script reads input files from a designated directory (`fit_input/`), which contain data for the NN interaction. Generate this data using the generate_interaction_mesh subroutine in the fitting_2b module in the parent code.
- **Output**: It writes the fitting parameters to a directory (`lmfit_parameters/`), creating files that represent the fitted interaction. This file is read by the inter_2b_fit_quadrature and inter_2bA_fit_quadrature subroutines in the parent code.

## Requirements

- Python 3.x
- Necessary Python libraries (e.g., `numpy`, `scipy`, `multiprocessing`)


## Example Usage

You can call the fitting functions directly from the command line using Python's `-c` flag.

For example, to run a parallel fitting over a range of J values:

python -c 'import fitting_lmfit_cost as flc; flc.fit_parallel(0, 4, 40, "np")'

This command fits the interaction for J values between 0 and 4, using 40 mesh points on each axis and an isospin projection of `'np'`.

## Example results

example input and output are stored in the example folder. Copy the content of Example/fit_input into the fit_input/ directory and run the code.
You should get results similar to those in Example/lmfit_parameters and Example/lmfit_graphs
## Output

After running any of the fitting functions, the output files containing the fitted parameters will be stored in the `lmfit_parameters/` directory.

## Notes

- The fitting procedure is designed to be flexible and customizable.
- You can modify the `run_fit_2` function if you want to experiment with different fitting methods.
- The script can be run in parallel to reduce computation time, especially when working with large datasets or multiple isospin projections.

## License

This project is licensed under the MIT License — see the LICENSE file for details.
