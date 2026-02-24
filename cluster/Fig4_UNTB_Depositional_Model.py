# Run a grid search through all parameter combinations in grid
# They will all get save to individual folders

# == Library Import == 
# Set directrory to src
# Run a simulation and save the output matrix.
src_paths = "../src/"

## == Internal Functions ==
save_path_sims = "../Simulation_Outputs_UNTB/UNTB_Dep_Cut"
save_path_ideal = "../Simulation_Outputs_UNTB/"

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
import untbPython
import samplingFunctions
import divDynFunctions
## Haar Fluctuation Analysis Functions
from haarFluctuationAnalysis import Haar_hebert
from crossHaarCorrelation import CHFC

# Other Imports
import pandas as pd
import numpy as np
from itertools import product
from numba import njit, set_num_threads

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

# Neutral Model
# Redefintion of UNTB function to include actual number of mutations per individual time step
@njit
def untb_Hankin_Phylo_V2(initial_population, mutation_probability, death_rate, generations, keep = False):
    # Initalize
    pop_size = len(initial_population) # Size of Jm
    
    # Allocate output matrix
    out_matrix = np.full(shape = (generations, pop_size), fill_value = np.nan)    
    a = initial_population.copy()
    out_matrix[0,:] = a # Time 1 is initial population

    # Speciation Tracker
    ansector = np.full(shape = (generations, death_rate), fill_value = np.nan)    
    descendent = np.full(shape = (generations, death_rate), fill_value = np.nan)   

    # Tracker number of new mutations
    total_mutations = [0]

    # Max Species Tracker
    max_species_id = np.max(a)

    # Simulate
    for i in range(1, generations):
    
        died = np.random.choice(pop_size, death_rate, replace = False) # index position of deaths
        mutated = np.random.uniform(0.0, 1, size=len(died)) < mutation_probability # probability new mutation arises
        n_mutations = np.sum(mutated) # new species
        n_deaths = np.sum(~mutated) # replaced species
        total_mutations.append(n_mutations)

        if n_deaths > 0:
            a[died[~mutated]] = np.random.choice(a, n_deaths, replace = True)
        if n_mutations > 0:
            new_species = np.arange(1, n_mutations + 1) + max_species_id 
            max_species_id += n_mutations

            ansector[i,0:n_mutations] = a[died[mutated]]
            a[died[mutated]] = new_species 
            descendent[i,0:n_mutations] = new_species
            
        out_matrix[i,:] = a
        
    return out_matrix, ansector, descendent, np.asarray(total_mutations)

## == Define Parameters ==
num_runs = 25

# Define Parameter Grid
pareto_cutoff = [100, 1000, 10000]
deposition_integration = [10, 100, 1000]
param_grid = list(product(pareto_cutoff, deposition_integration))

## == Begin Simulations ==
sim_length = 30000 # Running the simulation for 50,000 time steps
J = 1000 # Community size 
v = 0.05 # Speciation rate
d = 100
cutoffs = np.array([5000])
model_steps = np.max(cutoffs) + sim_length # Extending the simulation time
time_v = np.arange(cutoffs[0], cutoffs[0] + sim_length) - cutoffs[0]

# Simulate Neutral Ecology
expanded_community = np.repeat(1,J)
untb_results, ancs, desc, _ = untb_Hankin_Phylo_V2(expanded_community, mutation_probability = v, death_rate = d, generations = model_steps, keep = False)

# Make sure we have a simulation that doesn't go extinct
species_list_all = np.arange(1, np.max(untb_results)+1, dtype = np.int64)
set_num_threads(8) # Set thread number
spec_dyn = samplingFunctions.retrieve_species_dynamics_numba(untb_results, species_id = species_list_all)

## == Metrics for the perfectly sampled community using second-for-third == 
# Ideal Haar fluctuations for the community 
cut = cutoffs[0]
# Bias the output
test_spec = np.array(spec_dyn.copy(), dtype=np.int64)
presence = test_spec[cut:] # Discard burn-in period
# Get metrics
# Filter out time steps you will not observe
occ_list = convert_to_dictionary(presence)

# Diversity 
dynMat = divDynFunctions.binnedToExpandedFad(occ_list)
metricMatrix = divDynFunctions.counts(dynMat)
divRT = divDynFunctions.getdivRT(metricMatrix)

# Second for third extinction and origination
ext2f3 = divDynFunctions.getExt2f3(metricMatrix)
org2f3 = divDynFunctions.getOrg2f3(metricMatrix)

# Replace non-integer values with zero
ext2f3[np.isnan(ext2f3) | np.isinf(ext2f3)] = 0
org2f3[np.isnan(org2f3) | np.isinf(org2f3)] = 0

# New times are the bin transformed times
times = np.arange(0,len(divRT), dtype=np.float64)
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


filename = f"Ideal_Case_UNTB_Dep_Cut.npz"
np.savez_compressed(os.path.join(save_path_ideal, filename), ideal_data)


## == Run Model ==

# Run intial simulation for data before aggregation
time_v = np.arange(0,len(divRT_ideal))

# Set desktop path
results_summary = []

for cutoff, depositional_steps in param_grid:
    save_comb = []

    for run_id in range(1, num_runs + 1):
        try:
            print(f"Running simulation with cutoff={cutoff}, dep={depositional_steps}")

            ###
            # Retrive UNTB outputs
            test_spec = np.array(spec_dyn.copy(), dtype=np.int64)
            presence = test_spec[cutoffs[0]:]
            occ_list = convert_to_dictionary(presence[time_v,:])
            o_div = [len(i) for i in occ_list]

            #### -------------------------------------
            # Run depositional model
            obs_record, times, elevation = depositional_model_time(cutoffs = cutoff, max_time = 25000, deposition_integration = depositional_steps)

            # Bias the output
            # Set all arrays to empty if erosion occurs
            for i in range(len(obs_record)):
                state = obs_record[i]
                if state == False:
                    occ_list[i] = np.array([])

            # Filter based on hiatuses and plot
            o_div_filt = [len(i) for i in occ_list]
            non_zero_values = [val for val in o_div_filt if val > 0]
            #min_non_zero = min(non_zero_values) if non_zero_values else 0

            # Combine elements of array based on boundaries
            new_occ_list = []
            for i in range(len(times) - 1):
                agg_list = occ_list[int(times[i]):int(times[i+1])]
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
    filename = f"UNTB_cutoff{cutoff}_dep{depositional_steps}.npz"
    np.savez_compressed(os.path.join(save_path_sims, filename), save_comb)

    # --- Save summary ---
    results_summary.append({
        "cutoff": cutoff,
        "dep_rate": depositional_steps,
        "filename": filename,
        "type:": "Neutral Model"
    })
        
# Save summary DataFrame
summary_df = pd.DataFrame(results_summary)
summary_file = os.path.join(save_path_sims, "grid_search_summary.pkl")
summary_df.to_pickle(summary_file, compression="gzip")

print(f"Simulations complete. Summary saved to:\n{summary_file}")
