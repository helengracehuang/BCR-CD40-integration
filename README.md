# BCR-CD40-integration Model
"A molecular network model reveals non-monotonic integration of BCR and CD40 signals for controlling B-cell proliferation"

A combination of Receptor, NFkB, Proliferation, and Apoptosis modules.
We recommend running the model on a server with at least 32 threads. To run the model, one can run the main2.jl file. There are a few options one can set:
* `-v` for the type of multi-thread parallelization, where the option are: `"nonthread"`, `"thread"`, `"spawn"`. `"nonthread"` will run the model linearly and does not parallelize it, while `"thread"` will parallelize the model on static schedule, and `"spawn"` on dynamic schedule.
* `-o` for the destination of output .txt cell lineage file
* `-c` for the destination of output .jld signaling dynamics file (include output of nuclear RelA and cRel)
* `-i` for the destination of output steady state file. In the case it is combined with `-r`, the path will be used to reload from previous steady states, if the parameter distribution & pre-stimulation was already done. This can be used when you would like to rerun a simulation from a .jld file generated previously.
* `-r` for reloading from previous steady states

Example:
```
export JULIA_NUM_THREADS=64 # set number of threads to be used
home_dir="/path/to/dir/BCELL_PROJECT/"
modifier="lineages_125_CD40A_H62" # set the file name for outputs
julia $home_dir'scripts/main2.jl' -v "spawn" -o $home_dir'results/'$modifier'.txt' -i $home_dir'data/steady_'$modifier'.jld' -c $home_dir'data/cells_'$modifier'.jld' >> $home_dir'job-logs/'$modifier'.out'

```
## Scripts:
- `main3.jl`: function for running the simulations and saving results
- `ConstantParams2.jl`: constant parameters, including stimulus doses, stimulus delay, simulation time, scaling factors, etc.
- `ReactionRates3.jl`: reaction rate parameters for all module
- `ODE_Receptor5.jl`: ODE equations for BCR and CD40 receptor modules
- `ODE_NFkB3.jl`: ODE equations for NFkB module
- `ODE_Apoptosis2.jl`: ODE equations for Apoptosis module
- `ODE_Differentiation.jl`: ODE equations for Differentiation module
- `ODE_Proliferation.jl`: ODE equations for Cell Cycle module
- `SimulateFunctions4+.jl`: pre-simulation and simulation functions
- `HelperFunctions.jl`: helper functions for Michaelis-Menten and Hill functions, as well as parameter distributions

## Plotting scripts for each figure:

### Figure 1:

**Panel C,D**: `NFkB trajectories (dose response + composition).ipynb`

**Panel E,G**: `NFkB trajectories (ignore cell fates).ipynb` to plot NFkB trajectories from intermediate .jld files

### Figure 2:

**Panel B**: `calcModelFit.R` to plot the population dynamics by generation

**Panel C**: Excel sheets

**Panel F**: `calcModelFit.R` to calculate RMSD between model vs. experiment, then `RMSDheatmap.R` to plot from the tabulated results

### Figure 3:

**Panel A,C**: `calcModelFit.R` to plot the population dynamics by generation

**Panel B,D**: Excel sheets

**Panel E**: `calcModelFit.R` to calculate RMSD between model vs. experiment, then used Excel to tabulate the results and plot

**Panel F,G,I,J**: Excel sheets

**Panel H,K**: `calcReproducibility.R` to calculate RMSD between experimental replicates, then used Excel to tabulate the results and plot

### Figure 4:

**Panel B,D**: `calcModelFit.R` to plot the population dynamics by generation

**Panel C,E**: Excel sheets

**Panel F**: `calcModelFit.R` to calculate RMSD between model vs. experiment, then used Excel to tabulate the results and plot

**Panel G-J**: `NFkB trajectories (cell fates low vs high).ipynb` to plot NFkB trajectories from intermediate .jld files

### Figure 5:

**Panel A**: `calcModelFit.R` to plot the population dynamics by generation

**Panel B,C,D**: `plotTd2.R` for both Kaplan-Meier curve for each CD40 and BCR dose (B,C) and the bar graph of # survived cells at 24hrs (D)

**Panel E-J**: `plotFateLandscape.R` for the fate map. Need to adjust the code to specify with or without AICD

**Panel K-M**: `plotFateLandscapeDiff.R` to plot the difference between 2 fate maps from E-J

### Figure 6:

**Panel B**: `plotPopulationSize2.R` to plot the relative population size over time

**Panel C**: Excel sheets

### Figure 7:

**Panel B**: `BclXL trajectories (colored line).ipynb` to plot single-cell Bcl-xL trajectories, colored with caspase 8 level

**Panel C**: `BclXL trajectories (colored line).ipynb` to plot single-cell Bcl-xL trajectories, colored with NFkB level

**Panel D**: `BclXL trajectories (heatmap).ipynb` (last 8 blocks) to plot the violin plot of RelA, cRel, and Bcl-xL peak activity between dead and live cells

**Panel E-I**: `plotFateLandscape.R` for the fate map. Need to adjust the code to specify for with AICD

**Panel G-J**: `plotFateLandscapeDiff.R` to plot the difference between 2 fate maps from E-I
