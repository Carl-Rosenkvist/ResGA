0#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 00:04:15 2024

@author: carl
"""
import glob
import os
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm

from experiment_data_loader import ExperimentalDataLoader
from emulator import Emulator_fast

cross_dir = "../experimental_data/cross_sections"

reactions = ['K̅⁰+K⁰+p+p', 'K⁺+Σ⁺','K⁺+Σ⁻', 'K⁺+p+Λ', 'K⁺+n+Σ⁺', 'K⁰+p+Σ⁺']
#reactions = ['K⁺+Σ⁻','K⁰+Λ']

loader = ExperimentalDataLoader(cross_dir)

gens_path = "../generations"
df = pd.read_csv(gens_path + "/scores.csv")




# Group by 'generation' and find the maximum score for each generation
best = df.loc[df.groupby('generation')['score'].idxmax()]

plt.figure(figsize=(10, 6))
generations = best.generation

plt.plot(generations,best.score,label="score")
#reaction =reactions[3]
reaction = "score"
plt.plot(generations,best[reaction],label=reaction)
plt.xlabel('Generation')
plt.ylabel('Score')
#plt.yscale("log")
plt.title('Scores Across Generations')
plt.legend()
plt.show()


# Show the plot




best = df[df.score ==max(df.score)]
name = best.name.values[0]
gen = best.generation.values[0]


worst = df[df.score ==min(df.score)]
worst_name = worst.name.values[0]
worst_gen = worst.generation.values[0]


# Group the DataFrame by 'generation' and find the best performer within each group
best_performers = df.loc[df.groupby('generation')['score'].idxmax()].iloc[-10:]

# Extract the 'name' and 'generation' columns from the best performers
best_names_and_gens = best_performers[['name', 'generation']]

# Convert the result to a list of tuples
best_names_and_gens_list = list(best_names_and_gens.itertuples(index=False, name=None))
paths = [gens_path+"/generation-"+str(gen)+"/"+name for name,gen in best_names_and_gens_list]
best_path = gens_path+"/generation-"+str(gen)+"/"+name
# Get all CSV files in the directory



def getDATA(path):
    data_files = glob.glob(path + "/*.csv")
    # Initialize an empty dictionary to store numpy arrays
    data_dict = {}
    reactions = []
    # Read each CSV file into a numpy array and store in dictionary
    for file in data_files:
        filename = file[len(path)+1:-4]  # Extract filename without .csv extension
        reactions.append(filename)
        #data = np.loadtxt(file, delimiter=',', dtype=float)  # Read data into numpy array
        data = pd.read_csv(file)
        data_dict[filename] =data
        
    return data_dict


datas = [getDATA(path) for path in paths]
default = getDATA("../../genfit/default")
best = getDATA(best_path)

# Create a colormap
cmap = cm.get_cmap('viridis', len(datas[-5:]))

        





def ploting(label1,label2):
    for reaction in best.keys():
        plt.figure()
        plt.title(reaction)
        energies = loader.data[reaction].ssqrt.values
        sigs = loader.data[reaction].sigma.values 
        yerror = loader.data[reaction].yerror
        
        plt.scatter(default[reaction].ssqrt,default[reaction].sigma,label=label1,color="blue",marker="^")
        plt.scatter(best[reaction].ssqrt,best[reaction].sigma,label=label2,color="green",marker="x")

        plt.grid(True)

        plt.errorbar(energies, sigs, yerr=yerror, fmt='o', label="Data", capsize=2,color="red")

        plt.legend()
        plt.show()

    

    return 

ploting("SMASH-3.2","SMASH-3.2-reversed")
