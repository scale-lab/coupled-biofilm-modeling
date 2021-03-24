# Coupled Biofilm Modelling

This repository contains python code to perform mathematical modelling on coupled biofilms. We start with a baseline two biofilm coupled model that is a copy of the one proposed in the supplemental papers for the [*Coupling between distant biofilms and emergence of nutrient time-sharing* paper](https://science.sciencemag.org/content/356/6338/638). While we start from this baseline, this repository tries to build on this two biofilm model and try and model N number of biofilms and how they may be coupled. These models will ultimately help us build experiments for creating logic gates and other forms of computation using the coupled nature of biofilms and their oscillating behavior.

### Repository Contents

- **parameters.py:** this file contains the hyperparameters defined in the linked paper above and at the moment, is exactly the same as Table 1 in the supplemental notes.
- **functions.py**: this file contains some basic functions that will be necessary in all modelling of biofilms and is again, derived from the supplemental notes of the paper above.
- **two_biofilm_coupling.py**: this file takes the previous two files and models the differential equations described in the supplemental notes of the paper above for two coupled biofilms represented as modified Kuramoto oscillators.
- **three_biofilm_coupling.py**: this files is WIP, but is a light extrapolation of the two biofilm model.

### Running the Repository 

1. From Github, clone the repository and install the dependencies listed below
2. Inside of the repository, navigate to **two_biofilm_coupling.py** which will output a variety of graphs that describes the coupling behavior between the two biofilms specified in the file.
3. Under the comment, Initial conditions, we can change out the initial conditions that we pass into our two biofilm system. We can also change the values in the **parameters.py** file which contain some hyperparameters defined by the above paper. This is how we can emulate different situations that we can influence.

### Dependencies

- An installation of Python
- Uses packages: **matplotlib, seaborn, numpy**



