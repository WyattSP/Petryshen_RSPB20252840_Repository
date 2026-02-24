# Functions adapted from the R package described in Hankin, RKS. 2007. Introducing untb: an R package for simulating ecological drift under the unified neutral theory of biodiversity.
# Run in R to keep all analysis unified
# https://github.com/RobinHankin/untb/blob/master/R/untb.R
# Implementation is about x10 the R function with no optmization
import numpy as np
from numba import njit

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

@njit
def return_richness(untb_results):
    r = untb_results.shape
    # Save vector
    tmp = np.full(shape = (r[0], 2), fill_value=np.nan)
    for i in range(r[0]):
        tmp[i, 0] = i
        tmp[i, 1] = len(np.unique(untb_results[i,:]))
    return tmp

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