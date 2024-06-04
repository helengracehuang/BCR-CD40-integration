# For sequential stimulation, make BCR-induced IKK stop at 1,3,5,8hrs but stop AICD effects later due to delay
# Include libraries / packages required
#--------------------------------------------
using DifferentialEquations;

using ArgParse; # Argument parsing

# Benchmarking & profiling packages
using BenchmarkTools;
using StatProfilerHTML;

# Visualization packages
using DataFrames;
using Plots;
using Gadfly;
using Cairo;
using JLD;

# Include source files
#--------------------------------------------
include("ReactionRates3.jl");
include("HelperFunctions.jl");

include("ODE_Receptor5.jl"); # CHANGED from ReactionRates3.jl
include("ODE_NFkB3.jl");
include("ODE_Apoptosis2.jl");
include("ODE_Differentiation.jl");
include("ODE_Proliferation.jl");

include("SimulateFunctions4+.jl");

# Argument parsing w/ command line input
#--------------------------------------------
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--version", "-v"
            help = "nonthread, thread (static schedule), or spawn (dynamic schedule) version of lineage simulation"
            arg_type = String
            default = "nonthread"
            required = true
        "--initial", "-i"
            help = "reaction rates and steady states for each starting cell after initialization (.jld format)"
            arg_type = String
            default = "initial.jld"
            required = true
        "--cells", "-c"
            help = "simulated outputs from each cell (.jld format)"
            arg_type = String
            default = "cells.jld"
            required = true
        "--output", "-o"
            help = "output file name for simulated cell lineages"
            arg_type = String
            default = "output.txt"
            required = true
        "--param", "-p"
            help = "whether to load a parameter set from a .jl file"
            arg_type = String
            default = ""
            required = false
        "--reload", "-r"
            help = "whether to reload from previous steady states, provided by the --initial command"
            arg_type = Bool
            action = :store_true
            default = false
            required = false
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
const version = get(parsed_args, "version", "nonthread");
const steady_fn = get(parsed_args, "initial", "initial.jld");
const cells_fn = get(parsed_args, "cells", "cells.jld");
const output_fn = get(parsed_args, "output", "output.txt");
const param_fn = get(parsed_args, "param", "");
const reload = get(parsed_args, "reload", false);

if (param_fn != "")
    include(param_fn);
    print("Parameter set loaded!");
end

# Define standard reaction rates
#--------------------------------------------
rates = Matrix{Float64}(undef, TOTAL_SPECIES, TOTAL_PARAMS);
setAllRates!(rates);
const Srates = rates;

# Define ODEs for network
#--------------------------------------------
# Attempt to speed up the solver by splitting into 2 parts
# time-independent ODE (pre-simulation: phase = 1, only receptor & NFkB needs steady state simulation)
function computeNetworkNettFluxes!(nettFlux, concentration, (Srates, reactionFlux), time)
    computeReceptorNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    computeNFkBNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    # computeProlifNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1); #!!!!!!! remove if burn-in is taking too long
    nothing
end

function computeNetworkNettFluxes_AP1!(nettFlux, concentration, (Srates, reactionFlux), time)
    computeApoptosisNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    computeProlifNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    nothing
end

# time-dependent ODE (simulation, w/ delay: phase = 2)
function computeNetworkNettFluxes!(nettFlux, concentration, delay, (birthday, Srates, reactionFlux, historicFlux), time)
    computeReceptorNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2, time); #*** changed to add delay to receptor module *** add "delay, historicFlux, " if using Receptor4.jl
    computeNFkBNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2, delay, historicFlux, time);
    # computeDiffNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2);
    nothing
end

function computeNetworkNettFluxes_AP2!(nettFlux, concentration, (birthday, Srates, reactionFlux, inputCurves), time)
    computeApoptosisNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2, time, birthday, inputCurves);
    computeProlifNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2);
    nothing
end


# Callback function to detect cell death, mitotic and differentiation events
#----------------------------------------------------------------------
function condition(out, u, t, integrator)
    out[1] = u[CPARP] - 2500;
    out[2] = u[CDH1] - 0.2;
    # out[3] = u[IRF4] - 0.65*u[BCL6] - 1.2;
    out[3] = t - CD40L_DELAY;
    out[4] = t - BCR_END; # CHANGED
end

function affect!(integrator, index) # change back to terminate! after NFkB trajectory tracking
    if (index == 1)
        integrator.u[TOTAL_SPECIES] = 1;
        terminate!(integrator);
    elseif (index == 2)
        if (integrator.u[CYCB] > 2.0 && integrator.u[GEN] < MAX_GEN)
            integrator.u[MASS] /= 2.0;
            integrator.u[GEN] += 1.0;
            integrator.u[TOTAL_SPECIES] = 2;
            terminate!(integrator);
        end
    # elseif (index == 3) # remove effects of differentiation by commenting out termination
    #     integrator.u[TOTAL_SPECIES] = 3;
    #     # terminate!(integrator);
    elseif (index == 3)
        integrator.u[CD40L] = CD40L_DOSE;
        # integrator.u[ANTIGEN] = 0; # CHANGED
    elseif (index == 4) # CHANGED
        integrator.u[ANTIGEN] = 0; # CHANGED
    end
    nothing
end

cellFate = VectorContinuousCallback(condition, affect!, nothing, 4, save_positions=(true, true)); # CHANGED

# Set up a structure to hold cells
#-------------------------------------------------------------------------
allCells = Vector{Cell}(undef, FOUNDER_CELL_NUM);
if reload
    allCells = load(steady_fn, "allCells");
    print("Old steady states loaded!")
else
    initializeFounderCells!(Srates, allCells);
    JLD.save(steady_fn, "allCells", allCells);
end

# Define default delay functions (t < tau)
#--------------------------------------------
# const delay(historicFlux, p, t) = (historicFlux .= 0.0);
const delay(p, t; idxs=nothing) = typeof(idxs) <: Number ? 0.0 : zeros(TOTAL_SPECIES);

# Output parameters
#--------------------------------------------
print("ANTIGEN_DOSE:", '\t', ANTIGEN_DOSE, '\n');
print("CD40L_DOSE:", '\t', CD40L_DOSE, '\n');
print("CD40L_DELAY:", '\t', CD40L_DELAY, '\n');
print("BCR_DEATH:", '\t', BCR_DEATH, '\n');
print("IKK_MOD:", '\t', IKK_MOD, '\n');
print("MYCTHR:", '\t', MYCTHR, '\n');
print("BCL2THR:", '\t', BCL2THR, '\n');
# print("CYCDTHR:", '\t', CYCDTHR, '\n');
print("GROWTHR:", '\t', GROWTHR, '\n');
print("MAX_GEN:", '\t', MAX_GEN, '\n');

# Simulate cell lineages
#--------------------------------------------
if version == "nonthread"
    @time Simulate_nonthreaded!(allCells, Srates, delay);
elseif version == "thread"
    @time Simulate_threaded!(allCells, Srates, delay);
elseif version == "spawn"
    @time Simulate_spawned!(allCells, Srates, delay);
else
    print("Please input the correct version -v: nonthread, thread, or spawn.");
end

# @time Simulate_nonthreaded!(allCells, Srates, delay);
# @time Simulate_threaded!(allCells, Srates, delay);
# @time Simulate_spawned!(allCells, Srates, delay);
JLD.save(cells_fn, "allCells", allCells);

# Output information about all cells (for visualization)
#--------------------------------------------
out = open(output_fn, "w");
write(out, "birthday", '\t', "current_idx", '\t', "parent_idx", '\t', "generation", '\t', "fate", '\t', "fate_t", '\t', "abs_fate_t", '\t', "daughter_1_idx", '\t', "daughter_2_idx", '\n');
for i in 1:lastindex(allCells)
    write(out, string(allCells[i].birthday), '\t', string(allCells[i].current_idx), '\t', string(allCells[i].parent_idx), '\t', string(allCells[i].generation), '\t', string(allCells[i].fate), '\t', string(allCells[i].fate_t), '\t', string(allCells[i].abs_fate_t), '\t', string(allCells[i].daughter_1_idx), '\t', string(allCells[i].daughter_2_idx), '\n');
end
close(out);
