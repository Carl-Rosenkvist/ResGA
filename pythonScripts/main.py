#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:44:36 2023

@author: Carl Rosenkvist
"""
from multiprocessing import Pool
from genetic import *
from emulator import Emulator_fast
import os
import argparse
import numpy as np
import math

def parse_args():
    parser = argparse.ArgumentParser(description='Genetic algorithm parameters')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--max_generations', type=int, default=2, help='Maximum number of generations')
    parser.add_argument('--population_size', type=int, default=4, help='Population size')
    parser.add_argument('--shave', type=float, default=0.5, help='Shave value')
    parser.add_argument('--wd', type=str, default='/Users/carl/Phd/genfit', help='Working directory')
    parser.add_argument('--score_path', type=str,default="../generations/scores.csv", help='you must pick a place to save the scores')
    parser.add_argument('--case', type=str,default="default", help='Chose between using all genes or just mass or width')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    max_generations = args.max_generations
    population_size = args.population_size
    shave = args.shave
    threads = args.threads
    cross_dir = "../experimental_data/cross_sections"
    emulator_path ="../Emulator-fast"
    wd = args.wd
    score_path = args.score_path
    case = args.case
    
    initial_rate = 0.9 # Starting at 90%
    final_rate = 0.1    # Ending at 10%
    decay_constant = 0.02  # Adjust this to control the rate of decay

    
    # Start the genetic algorithm with generation 0
    start = 0
    if not os.listdir("../generations"):
        population, names = initPopulation(population_size,case)
    else:
        population, names = chromosFromNewest()
        start = int(newestGeneration().split("-")[-1])
        
        
    eml = Emulator_fast(emulator_path, cross_dir)
    reactions = eml.reactions_pp + eml.reactions_piminus_proton+eml.reactions_piplus_proton
    if(start == 0):
        with open(score_path, 'w') as f:
            header = 'generation,name,score,' + ','.join(reactions) + '\n'
            f.write(header)
            
            
            

    children = population
    parents = []

    for generation in range(start, max_generations):
        if case == "default" or case == "branch+mass":
            generation_path = f'{wd}/generations/generation-{generation}'
        else:
            generation_path = f'{wd}/stability/{case}'

        if not os.path.exists(generation_path):
            os.makedirs(generation_path)
        makeOrder(children, generation_path,case)

        
        updated_population = eml.process_chromos(
            children,
            generation_path,
            threads
            )

        # Replace the original population with the updated one
        population = updated_population

    
        score_board = [(chromo.name, chromo.score,chromo.scores) for chromo in population]
        with open(score_path, 'a') as f:
            for i in score_board:
                line = f'{generation},{i[0]},{i[1]}'
                for reaction in reactions:
                    line+= ","+str(i[2][reaction]) 
                f.write(line+"\n")

        population += parents
        if(max_generations > 1):
            #mutation_rate = decaying_mutation_rate(initial_rate, final_rate, decay_constant, generation)
            #print(mutation_rate)
            mutation_rate = 0.1
            parents, children = nextPopulation(generation_path, population_size, shave, population,mutation_rate)
