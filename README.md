# access_idealised
Create input files for idealised ACCESS-ESM1.5 runs

# Generate inputs:
./generate_inputs.sh

First, specify an experiment name. Default is solo150.
This script generates the input files for the ocean, atmosphere, ice and coupler. 
It takes inputs from previously configured idealised inputs generated by Dietmar Dommenget.

**Need to reconfigure the CODEDIR and ANCILDIR if setting up for a new user.**

# Create run directory
./make_rundirs.sh

Again, specify experiment name. Default is solo150.

**Need to reconfigure the ANCILDIR and RUNDIR if setting up for a new user.**

This sets up the run directories for a new simulation.
Will be configured in:

~/ACCESS/$EXP 

(where $EXP=solo150 by default)

