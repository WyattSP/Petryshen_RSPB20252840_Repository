# Run a grid search through all parameter combinations in grid
# They will all get save to individual folders

# == Library Import == 
# Set directrory to src
# Run a simulation and save the output matrix.
src_paths = "../src/"


## == Internal Functions ==
save_path_sims = "../Simulation_Outputs_RW/RW_Dep_Cut"
save_path_ideal = "../Simulation_Outputs_RW/"

# Set paths 
import os
import sys
os.chdir(src_paths)
sys.path.append(os.getcwd()) 

# From src
## Diversity Metric Functions
from diversityMetrics import *
import divDynFunctions
from pareto import *

## Haar Fluctuation Analysis Functions
from haarFluctuationAnalysis import Haar_hebert
from crossHaarCorrelation import CHFC

# Other Imports
import pandas as pd
import numpy as np
from itertools import product

# Depositional Model
def depositional_model_time(cutoffs, max_time, deposition_integration):

    # Externally Set Parameters
    cutoff = cutoffs # This is what we will be setting for the pareto_cutoff

    # Internally Set Parameters
    gamma = 0.4 # Shape parameter used in Schumer and Jerolmack 2009
    xmin = 1 # minimum hiatus duration (year); using integer values so must be int
    deposition_rate = 1000  # cm/yr
    deposition_duration = deposition_integration  # 1000 years
    deposition_thickness = deposition_rate * deposition_duration  # deposition = 10 cm
    max_time = max_time  # max total simulated time in kyr (stop condition); Base on selection from UNTB simulation

    # Set Save Directories
    times = [0]
    elevation = [0]
    obs_record =[]
    t = 0
    s = 0
    step = 0

    while t < max_time:
        # Setting the hiatuses
        while True:
            hiatus = rpareto(n=1, lam=gamma, a=xmin).astype(int)[0]
            if hiatus <= cutoff:
                break

        # Check if next step would exceed max_time
        if t + hiatus + deposition_duration > max_time:
            break 

        # Record dynamics within step
        t += hiatus
        obs_record.extend([False] * hiatus)  # Append erosion / nondeposition
        t += deposition_duration
        obs_record.extend([True] * deposition_duration)  # Append deposition
        s += deposition_thickness
        times.append(t)
        elevation.append(s)
        step += 1

    # Will return up to the while statement
    return np.asarray(obs_record), np.asarray(times), np.asarray(elevation)

# Density-Independent Speciation Model
# Density Independent Neutral Model
def density_independent_model(K, steps, speciation_pr, neutral_pr, extinction_pr):
    # Check probability sum
    # P(0) + P(1) + P(2) = 1
    p_sum = speciation_pr + neutral_pr + extinction_pr
    if not np.isclose(p_sum, 1.0):
        raise ValueError("Probabilities incorrectly set")

    # Itialize community of K species 
    i_comm = np.arange(1, K + 1)  # Species labeled from 1 to K
    community = [i_comm.copy()]  # Start with initial generation

    # Iterate through i steps 
    for _ in range(steps - 1):
        current_comm = community[-1]
        m_species = np.unique(current_comm) # You may be able to drop this

        if len(m_species) > 0:
            # Draw single event per species
            events = np.random.choice([0, 1, 2], size=K, p=[speciation_pr, neutral_pr, extinction_pr])
            bins = np.bincount(events, minlength=3) # bins[0] = speciation; bins[1] = nothing; bins[2] = extinction;
            n_speciation = bins[0]
            n_extinction = bins[2]

            # Extinctions events
            if n_extinction >= len(m_species):
                surviving_species = np.array([], dtype=int)
            else:
                surviving_species = np.random.choice(m_species, size=len(m_species) - n_extinction, replace=False)
            
            # Speciation events
            max_sp = max(m_species) if len(m_species) > 0 else 0
            new_species = np.arange(max_sp + 1, max_sp + 1 + n_speciation)

            #Next generation
            next_comm = np.concatenate((surviving_species, new_species))
            community.append(next_comm)
        else:
            # All species went extinct, keep empty generation
            community.append(np.array([], dtype=int))

    return community

## == Define Parameters ==
num_runs = 25

# Define Parameter Grid
pareto_cutoff = [100, 1000, 10000]
deposition_integration = [10, 100, 1000]
param_grid = list(product(pareto_cutoff, deposition_integration))

## == Begin Simulations ==

## == Set Speciation and Extinction Rates ==
speciation_rate = np.round(0.005,3)
extinction_rate = np.round(0.005,3)

# Make sure we have a simulation that doesn't go extinct
while True:
    RW_Output = density_independent_model(250, 25000, speciation_rate, 0.99, extinction_rate)
    richness = [len(i) for i in RW_Output]

    if 0 in richness:
        print("Zero richness found — rerunning simulation...")
        continue
    else:
        break

print(f'Speciation Rate: {speciation_rate}; Extinction Rate: {extinction_rate}')

## == Metrics for the perfectly sampled community using second-for-third == 
# Ideal Haar fluctuations for the community 
comm_true = RW_Output.copy()
# Diversity 
dynMat = divDynFunctions.binnedToExpandedFad(comm_true)
metricMatrix = divDynFunctions.counts(dynMat)
divRT = divDynFunctions.getdivRT(metricMatrix)

# Second for third extinction and origination
ext2f3 = divDynFunctions.getExt2f3(metricMatrix)
org2f3 = divDynFunctions.getOrg2f3(metricMatrix)

# Replace non-integer values with zero
ext2f3[np.isnan(ext2f3) | np.isinf(ext2f3)] = 0
org2f3[np.isnan(org2f3) | np.isinf(org2f3)] = 0

# New times are the bin transformed times
times = np.arange(0,len(richness))
org2f3_ideal = org2f3[1:-1]
ext2f3_ideal = ext2f3[1:-1]
divRT_ideal = divRT[1:-1]

# Get midpoints between stratigraphic intervals
mid_point_times = (times[:-1] + np.diff(times)/2)
time_an_ideal = mid_point_times[1:-1]

# Haar Fluctuation Analysis
Thd, full_Div, _ = Haar_hebert(np.array(divRT_ideal), time_an_ideal, Scales = np.array([0.1, 4]), q = 1, overlap = 0, returnFull = True, allq = False)
Tex, full_extRT, _ = Haar_hebert(np.array(ext2f3_ideal) * 10, time_an_ideal, Scales = np.array([0.1, 4]), q = 1, overlap = 0, returnFull = True, allq = False)
Tor, full_oriRT, _ = Haar_hebert(np.array(org2f3_ideal) * 10, time_an_ideal, Scales = np.array([0.1, 4]), q = 1, overlap = 0, returnFull = True, allq = False)
time_axis = Thd[:,0]
# Correlation Origination/Extinction
true_corr_org_ext = CHFC(np.array(full_oriRT),np.array(full_extRT))
true_corr_org_div = CHFC(np.array(full_oriRT),np.array(full_Div))
true_corr_ext_div = CHFC(np.array(full_extRT),np.array(full_Div))

# Save this output
ideal_data = {
    "org2f3": org2f3_ideal,
    "ext2f3": ext2f3_ideal,
    "divRT": divRT_ideal,
    "mid_point_times": mid_point_times,
    "time_an": time_an_ideal,
    "haar_Div": Thd,
    "haar_oriRT": Tex,
    "haar_extRT": Tor,
    "corr_org_ext": true_corr_org_ext,
    "corr_org_div": true_corr_org_div,
    "corr_ext_div": true_corr_ext_div
}


filename = f"Ideal_Case_RW_Dep_Cut.npz"
np.savez_compressed(os.path.join(save_path_ideal, filename), ideal_data)


## == Run Model ==

# Run intial simulation for data before aggregation
time_v = np.arange(0,len(richness))

# Set desktop path
results_summary = []

for cutoff, depositional_steps in param_grid:
    save_comb = []

    for run_id in range(1, num_runs + 1):
        try:
            print(f"Running simulation with cutoff={cutoff}, dep={depositional_steps}")


            #### -------------------------------------
            # Run depositional model
            obs_record, times, elevation = depositional_model_time(cutoffs = cutoff, max_time = 25000, deposition_integration = depositional_steps)

            # Plug in Random Walk model here
            comm_model = RW_Output.copy()
            o_div = [len(i) for i in comm_model]
            # Set all arrays to empty if erosion occurs
            for i in range(len(obs_record)):
                state = obs_record[i]
                if state == False:
                    comm_model[i] = np.array([])

            # Filter based on hiatuses and plot
            o_div_filt = [len(i) for i in comm_model]
            non_zero_values = [val for val in o_div_filt if val > 0]
            #min_non_zero = min(non_zero_values) if non_zero_values else 0

            # Combine elements of array based on boundaries
            new_occ_list = []
            for i in range(len(times) - 1):
                agg_list = comm_model[int(times[i]):int(times[i+1])]
                if len(agg_list) <= 1:
                    new_occ_list.append(np.array(list(agg_list)))
                else:
                    aggregated_beds = np.unique(np.concatenate(agg_list)) 
                    new_occ_list.append(aggregated_beds)

            # Diversity Metrics
            # Diversity 
            dynMat = divDynFunctions.binnedToExpandedFad(new_occ_list)
            # Convert in matrix
            metricMatrix = divDynFunctions.counts(dynMat)
            # Diversity metrics
            divRT = divDynFunctions.getdivRT(metricMatrix)
            
            # Second for third extinction and origination
            ext2f3 = divDynFunctions.getExt2f3(metricMatrix)
            org2f3 = divDynFunctions.getOrg2f3(metricMatrix)

            # Replace non-integer values with zero
            ext2f3[np.isnan(ext2f3) | np.isinf(ext2f3)] = 0
            org2f3[np.isnan(org2f3) | np.isinf(org2f3)] = 0

            # New times are the bin transformed times
            # Drop first and last values so no NAN present
            org2f3 = org2f3[1:-1]
            ext2f3 = ext2f3[1:-1]
            divRT = divRT[1:-1]

            # Get midpoints between stratigraphic intervals
            mid_point_times = (times[:-1] + np.diff(times)/2)
            time_an = mid_point_times[1:-1]

            # Haar Fluctuation Analysis
            # Haar Calculations
            haar_Div, full_Div, _ = Haar_hebert(np.array(divRT), time_an, Scales = np.array([0.1, 4]), q = 1, overlap = 0, returnFull = True, allq = False)
            haar_extRT, full_extRT, _ = Haar_hebert(np.array(ext2f3) * 10, time_an, Scales = np.array([0.1, 4]), q = 1, overlap = 0, returnFull = True, allq = False)
            haar_oriRT, full_oriRT, _ = Haar_hebert(np.array(org2f3) * 10, time_an, Scales = np.array([0.1, 4]), q = 1, overlap = 0, returnFull = True, allq = False)
            time_axis = haar_Div[:,0]
            # Correlation Origination/Extinction
            corr_org_ext = CHFC(np.array(full_oriRT),np.array(full_extRT))
            corr_org_div = CHFC(np.array(full_oriRT),np.array(full_Div))
            corr_ext_div = CHFC(np.array(full_extRT),np.array(full_Div))

            #### -------------------------------------
            # --- Save full outputs ---
            sim_data = {
                "org2f3": org2f3,
                "ext2f3": ext2f3,
                "divRT": divRT,
                "mid_point_times": mid_point_times,
                "time_an": time_an,
                "o_div_filt": o_div_filt,
                "o_div": o_div,
                "haar_Div": haar_Div,
                "haar_oriRT": haar_oriRT,
                "haar_extRT": haar_extRT,
                "corr_org_ext": corr_org_ext,
                "corr_org_div": corr_org_div,
                "corr_ext_div": corr_ext_div
            }

            # Append simulation iteration to list
            save_comb.append(sim_data)

        except Exception as e:
            print(f"Error for Pareto cutoff={cutoff}, depositional steps = {depositional_steps}, run={run_id}: {e}")

    # After 50 iterations Save full save_comb list
    filename = f"RW_cutoff{cutoff}_dep{depositional_steps}.npz"
    np.savez_compressed(os.path.join(save_path_sims, filename), save_comb)

    # --- Save summary ---
    results_summary.append({
        "cutoff": cutoff,
        "dep_rate": depositional_steps,
        "filename": filename,
        "type:": "Stochastic SP Model"
    })
        

# Save summary DataFrame
summary_df = pd.DataFrame(results_summary)
summary_file = os.path.join(save_path_sims, "grid_search_summary.pkl")
summary_df.to_pickle(summary_file, compression="gzip")

print(f"Simulations complete. Summary saved to:\n{summary_file}")
