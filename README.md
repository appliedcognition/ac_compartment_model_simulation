This repository contains the [mrgSolve](https://mrgsolve.org/) simulation script used for the MedRxiv manuscript "The glymphatic system clears amyloid beta and tau from brain to plasma in humans".

To run this script, install the missing R dependencies and run the script. We recommend using [RStudio](https://posit.co/download/rstudio-desktop/) to simplify the process. The output of the simulation is saved in the `output` folder.

The output consists in the following files:

* `model/biomarker_comp_model.cpp`: The mrgSolve model parameters.
* Plots `INCREASE_Ab40.pdf/tiff` and `DECREASE_Ab40.pdf/tiff`: depict the modelled effects of increasing and decreasing the Ab40 concentration.
* `df_compart_BIO_P2_SIMULATION.csv`: Shows all the simulation steps
* `df_compart_BIO_P2_SIMULATION_RATIOS.csv`: Shows all the simulation ratios
 
Copyright 2025 [Applied Cognition Inc](https://appliedcognition.com/).