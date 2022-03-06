# Coupled Biofilm Modelling

This repository contains python code to perform mathematical modelling on coupled biofilms. We start with a baseline two biofilm coupled model that is a direct copy of the proposed in the supplemental papers for the [*Coupling between distant biofilms and emergence of nutrient time-sharing* paper](https://science.sciencemag.org/content/356/6338/638). While we start from this baseline, this repository tries to build on this two biofilm model and model 'N' number of biofilms and how they may be coupled. These models will ultimately help us build experiments for creating logic gates and other forms of computation using the coupled nature of biofilms and their oscillating behavior.

At the moment, the repository contains a two-biofilm model that is a recreation of the *Coupled Biofilm Paper* and a slightly extrapolated model that uses three biofilms instead of two biofilms in a linear orientation. We can also accomodate weighted adjacency matrices, in which it is possible
to construct a differential equation model around them given the matrix and hyperparameters to 
describe the state of the physical system.

### Dependencies

- An installation of Python
- Packages: **matplotlib, seaborn, numpy, pandas, scipy.integrate**

### Repository Contents

- **$ROOT/analysis/**: This folder contains python scripts that analyze CSV files produced in the /output files. These are generally pandas-heavy files that manipulate .CSV file data for data science applications and graphing sake.
- **$ROOT/img/:** This folder contains all of the media inside of this repo and is labelled accordingly. At the moment, it contains gifs of different animations, gate level descriptions, and other visualizations that have been used in order to better understand what we are working on.
- **$ROOT/output/:** This folder contains the output of any of the python scripts that produce anything. That is, it contains the .csv files generated from the sweeps, the cleaned .csv files of those .csv files and more. 
- **$ROOT/sweep/:** This folder primarily contains scripts that are used in order to sweep a certain set of parameters in an attempt to get a better understanding of what is going on and these scripts produce .csv files for analysis later on.
- **parameters.py:** this file contains the hyperparameters defined in the linked paper above and at the moment, is exactly the same as Table 1 in the supplemental notes.
- **functions.py**: this file contains some basic functions that will be necessary in all modelling of biofilms and is again, derived from the supplemental notes of the paper above.
- **two_biofilm_coupling.py**: this file takes the previous two files and models the differential equations described in the supplemental notes of the paper above for two coupled biofilms represented as modified Kuramoto oscillators.
- **three_biofilm_coupling.py**: Initial three biofilm coupling model using nearest neighbors Kuramoto coupling. This file functions off of the functions and parameters defined in parameters.py and functions.py. It functions as a manual way of running one particular configuration of our model.
- **n_biofilms_coupling.py**: This is an extrapolation of the three biofilm coupling system into 
n biofilms! We use a weighted adjacency matrix to describe the graph, in which it is parsed and turned
into a differential equation model by our functions. 
- **animate.py**: This file contains the script that creates .gif files of the three biofilm animations that we can use in order to better visualize what is going on. This is done using the animation functions inside of matplotlib and are outputted into the /output folder.

### Running the Repository 

1. From Github, clone the repository and install the dependencies listed below
2. Inside of the repository, navigate to **two_biofilm_coupling.py** which will output a variety of graphs that describes the coupling behavior between the two biofilms specified in the file.
3. Under the comment, Initial conditions, we can change out the initial conditions that we pass into our two biofilm system. We can also change the values in the **parameters.py** file which contain some hyperparameters defined by the above paper. This is how we can emulate different situations that we can influence. 
4. Feel free to go around and change some of the hyperparameters in the **three_biofilm_coupling.py** 
or any of the biofilm coupling models.
5. For the csv_analysis and the graphing, all we really need to do to be able to run it on new .csvs created from the three_biofilm_testing is to change the filename global variables in the files.

