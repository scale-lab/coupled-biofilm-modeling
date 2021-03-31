# Coupled Biofilm Modelling

This repository contains python code to perform mathematical modelling on coupled biofilms. We start with a baseline two biofilm coupled model that is a copy of the one proposed in the supplemental papers for the [*Coupling between distant biofilms and emergence of nutrient time-sharing* paper](https://science.sciencemag.org/content/356/6338/638). While we start from this baseline, this repository tries to build on this two biofilm model and try and model N number of biofilms and how they may be coupled. These models will ultimately help us build experiments for creating logic gates and other forms of computation using the coupled nature of biofilms and their oscillating behavior.

![gate_level_description](https://github.com/scale-lab/coupled-biofilm-modeling/blob/main/img/gate_level_description.PNG)

**Figure 1**: Initial gate level description of three biofilm model with three hyperparameters (g, k, delta) given the assumption that we can, at runtime, set the phases of biofilm 1 and biofilm 3.

### Dependencies

- An installation of Python
- Uses packages: **matplotlib, seaborn, numpy**
- 

### Repository Contents

- **parameters.py:** this file contains the hyperparameters defined in the linked paper above and at the moment, is exactly the same as Table 1 in the supplemental notes.
- **functions.py**: this file contains some basic functions that will be necessary in all modelling of biofilms and is again, derived from the supplemental notes of the paper above.
- **two_biofilm_coupling.py**: this file takes the previous two files and models the differential equations described in the supplemental notes of the paper above for two coupled biofilms represented as modified Kuramoto oscillators.
- **three_biofilm_coupling.py**: Initial three biofilm coupling model using nearest neighbors Kuramoto coupling. This file functions off of the functions and parameters defined in parameters.py and functions.py. It functions as a manual way of running one particular configuration of our model.
- **three_biofilm_testing.py**: This is a file that takes the three_biofilm_coupling model from the previous file and then runs it iteratively over hyperparameters defined in the file and outputs a .CSV file containing key information characterizing a particular communication strength, competition strength and glutamate concentration and initial phases between biofilms 1 and biofilms 3. Note that in this case, biofilm 2 is held constant at phase of 0 at the start of simulations
- **three_biofilm_csv_analysis.py**: This file takes in the previously created .csv file and runs it through a pandas dataframe to graph what data we have gathered. Right now, it takes the .csv file and separates the different initial phases and graphs the rest of the information separately on four different subplots. The phase differences for biofilms 1 and 2 at steady state then represent the binary color mapping of the graphs. In addition, at the end of the script, we also create a truncated, more informative CSV file from the raw data that uses a threshold to determine if a particular combination has biofilms that are in phase or out of phase. It then uses this binary encoded information on all four possible initial phase combinations to produce a 2 input truth table where we characterize logic gates off of this information which is written into the *_gates.csv file.
- **graph_three_biofilm_gates.py**: This is the file that takes the previous gates.csv file and then graphs based off of the three hyperparameters (g, delta, k) like before, but now combines all four of the different phase plots with a different colormap. This colormap dictates which points correspond to what logical gates. This graph is provided above.

### Running the Repository 

1. From Github, clone the repository and install the dependencies listed below
2. Inside of the repository, navigate to **two_biofilm_coupling.py** which will output a variety of graphs that describes the coupling behavior between the two biofilms specified in the file.
3. Under the comment, Initial conditions, we can change out the initial conditions that we pass into our two biofilm system. We can also change the values in the **parameters.py** file which contain some hyperparameters defined by the above paper. This is how we can emulate different situations that we can influence. 
4. Feel free to go around and change some of the hyperparameters in the **three_biofilm_testing.py** file to create a new .csv file containing information that may be different from what .csvs I already have produced
5. For the csv_analysis and the graphing, all we really need to do to be able to run it on new .csvs created from the three_biofilm_testing is to change the filename global variables in the files.

### Code Issues

- In the testing file, to determine phase difference, I only use the last 10 points in the solved phase vectors but in order to approach correct mean values, I should be using a larger amount of points, but this is dependent on the time vector.
- The testing file runs quite slowly but it is a highly parallelizable program, meaning that with a couple of threads, I should be able to cut down quite a bit on the time to run programs, meaning we can run much more complex programs
- Right now, the graphs and images that I have produced only look at the phase difference between biofilms 1 and 2. We need to determine in our computing system, what does it mean to be in phase and out of phase for the second biofilm. Maybe we the middle biofilm to be out of phase or in phase with both of the other biofilms? 

### Higher Level Conceptual Questions

- Right now, all of the testing and graphing I have done so far depends on the ability to set initial phases of biofilms which I am not sure whether that is possible, but something to think about
- At the moment, we determine some threshold for assigning some value to in phase or out of phase but I have just chosen a pretty arbitrary value (1 radian) at the moment without much thought. In the future, it may be a good idea to determine the best value for thresholding this value
  - This brings into the question that this type of computing is nonregenerative. That is, after one gate, in-phase might be something like 2.8 radians out of phase, where true out of phase is pi radians. We might have to think about this more in order to get any type of conventional computing that makes sense here.
- At the moment, this biofilm computation is kind of different from conventional computing in that inputs are destroyed when something is computed. This means that if multiple biofilms need to use the same computation, we run into huge issues. We also run into issues in that this input being related to the output and the output being related to the input makes the computation more of a graph rather than a tree, which is not typical of conventional computing
  - I wonder whether there are papers on graph type conventional computing? 
  - Maybe this is where we need to realize that this idea of gates and conventional computing applied to biofilm computing may not work as intended and that maybe we should instead be trying to analyze coupling of many biofilms instead to do something a little more unconventional.





