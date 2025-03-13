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
home_dir="/path/to/dir/BCR-CD40-integration/"
modifier="lineages_125_CD40A_H62" # set the file name for outputs
julia $home_dir'Simulation_scripts/main3.jl' -v "spawn" -o $home_dir'results/'$modifier'.txt' -i $home_dir'data/steady_'$modifier'.jld' -c $home_dir'data/cells_'$modifier'.jld' >> $home_dir'job-logs/'$modifier'.out'

```
## Simulation Scripts:
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

## Image Analysis Scripts (ImageJ Macro):
- `Select_AOI_segmentation.ijm`: semi-automate the cell segmentation process based on a selected area of interests (AOI). After automatic thresholding, the code prompts the user to manually inspect the segmentation and adjust if needed. Finally, fluorescence intensity values are measured for the RelA, cRel, and H2B channels for whole cell, nucleus, and cytoplasm.
    - **Here’s a pseudo-code description of all the processing steps:**
        1. Duplicate selected AOI region & build a dead cell mask from DRAQ7-APC channel (optional, as only some conditions stained for DRAQ7)
        2. Build a "WholeCell_Raw" image from 2 or 3 fluorescence channels
        3. Build a mask of the edges of the imaging grid
        4. Whole cell segmentation
            1. apply global auto threshold method (mostly **Otsu**, occasionally **Moments**)
            2. morphological processing: fill holes + watershed + erosion
            3. filter based on cell sizes and circularity
            4. check if ROI falls in black edge mask or dead cell mask
        5. Nuclear segmentation
            1. apply global auto threshold method (mostly **Moments**)
            2. morphological processing: fill holes + watershed + erosion
            3. check if ROI falls in black edge mask
        6. Match nucleus with best available whole cell
            1. has to overlap at least 30%
            2. delete unmatched nuclei and whole cells
        7. Clean up nuclei and whole cell, and compute cytoplasm
        8. Fluorescence intensity measurement
    - In the output **Results** table from `Select_AOI_segmentation.ijm`:
        - `Ch` = fluorescence channel
            - when 3 channels are open:
                - 1 = cRel-TFP
                - 2 = H2B-mCherry
                - 3 = RelA-mVenus
            - when 4 channels are open (DRAQ7 as far red), then DRAQ7 becomes 1 and the rest are moved fown
        - `Group` = cell morphology features
            - 0 = AOI selection (unimportant for quantification)
            - 1 = whole cell
            - 2 = nucleus
            - 3 = cytoplasm
        - `Label` contains both morphology information + cell ID, but can be modified to include only cell ID
- `Select_AOI_CellDeath.ijm`: nuclear segmentation based on H2B channel, and dead cell segmentation based on DRAQ7 channel. Final output is the number of total cells, the number of dead cells, and the percentage of dead cells out of all.
- `Export_RepPNG_wScale.ijm`: export representative images (281 x 281 pixel) in PNG format with scale bar (Figure S3). The representative image is brightness-thresholded based on median background intensity (minimum) and a batch correction factor (maximum). The script includes the min and max values for each channel at each timepoint.

## Plotting scripts for each figure:

Source data for all panels noted “Excel sheets” are available on BioStudies public database ([www.ebi.ac.uk/biostudies](http://www.ebi.ac.uk/biostudies)) . See  “Data Availability” section of the corresponding paper for accession numbers.

- **Figure 1:**
    - **Panel C,D**: `NFkB trajectories (dose response + composition).ipynb`
    - **Panel E,G**: `NFkB trajectories (ignore cell fates).ipynb` to plot NFkB trajectories from intermediate .jld files
    - **Panel F,H**: Excel sheets
- **Figure 2:**
    - **Panel B**: `calcModelFit.R` to plot the population dynamics by generation
    - **Panel C**: Excel sheets
    - **Panel D,E**: `calcModelFit.R` to calculate RMSD between model vs. experiment, then `RMSDheatmap.R` to plot the heatmap from Excel-tabulated results
- **Figure 3:**
    - **Panel A,C**: `calcModelFit.R` to plot the population dynamics by generation
    - **Panel B,D**: Excel sheets
    - **Panel E**: `calcModelFit.R` to calculate RMSD between model vs. experiment, then used Excel to tabulate the results and plot
    - **Panel F,G,I,J**: Excel sheets
    - **Panel H,K**: `calcReproducibility.R` to calculate RMSD between experimental replicates, then used Excel to tabulate the results and plot
    - **Panel L**: Excel sheets
    - **Panel M**: `plotImageAnalysis5.R` to plot cell area (from Microscopy) in line plots with error bars
- **Figure 4:**
    - **Panel B**: `calcModelFit.R` to plot the population dynamics by generation
    - **Panel C**: Excel sheets
    - **Panel D**: `calcModelFit.R` to calculate RMSD between model vs. experiment, then used Excel to tabulate the results and plot
    - **Panel E,G**: `NFkB trajectories (cell fates low vs high).ipynb` to plot NFkB trajectories from intermediate .jld files
    - **Panel F,H**: `plotImageAnalysis8.R` to plot RelA and cRel fluorescence (from Microscopy) in violin and line plots overlayed with Western blot quantification
- **Figure 5:**
    - **Panel A**: `calcModelFit.R` to plot the population dynamics by generation
    - **Panel B**: Excel sheets
    - **Panel C,D,E**: `plotTd2.R` for both Kaplan-Meier curve for each CD40 and BCR dose (C,D) and the bar graph of # survived cells at 24hrs (E)
    - **Panel F-K**: `plotFateLandscape.R` for the fate map. Need to adjust the code to specify with or without AICD
    - **Panel L-N**: `plotFateLandscapeDiff.R` to plot the difference between 2 fate maps from F-K
- **Figure 6:**
    - **Panel B**: `plotPopulationSize2.R` to plot the relative population size over time
    - **Panel C**: Excel sheets
- **Figure 7:**
    - **Panel B**: `BclXL trajectories (colored line).ipynb` to plot single-cell Bcl-xL trajectories, colored with caspase 8 level
    - **Panel C**: `BclXL trajectories (colored line).ipynb` to plot single-cell Bcl-xL trajectories, colored with NFkB level
    - **Panel D**: `BclXL trajectories (heatmap).ipynb` (last 8 blocks) to plot the violin plot of RelA, cRel, and Bcl-xL peak activity between dead and live cells
    - **Panel E-I**: `plotFateLandscape.R` for the fate map. Need to adjust the code to specify for with AICD
    - **Panel G-J**: `plotFateLandscapeDiff.R` to plot the difference between 2 fate maps from E-I
    - **Panel D**: `BclXL trajectories (heatmap).ipynb` (last 8 blocks) to plot the violin plot of RelA, cRel, and Bcl-xL peak activity between dead and live cells
    - **Panel E-I**: `plotFateLandscape.R` for the fate map. Need to adjust the code to specify for with AICD
    - **Panel G-J**: `plotFateLandscapeDiff.R` to plot the difference between 2 fate maps from E-I
