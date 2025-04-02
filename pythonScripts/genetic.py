#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 10:27:40 2023

@author: Carl Rosenkvist
"""

import numpy as np
import pandas as pd
from decayFILE import *
from score import chromoScore
import glob
import shutil
import random
import os
import subprocess
import math
ref_path = "./ref.csv"
orig_path="../decaymodes.txt"
abc='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
modes = ['Λ K','Σ K','N φ', 'N a₀(980)', "N f₀(980)"]

def decaying_mutation_rate(initial_rate, final_rate, decay_constant, generation):
    """
    Calculate the decaying mutation rate for a given generation.

    Parameters:
    initial_rate (float): The initial mutation rate (e.g., 1.0 for 100%).
    final_rate (float): The final mutation rate (e.g., 0.1 for 10%).
    decay_constant (float): The decay constant.
    generation (int): The current generation number.

    Returns:
    float: The mutation rate for the given generation.
    """
    return (initial_rate - final_rate) * math.exp(-decay_constant * generation) + final_rate


def randomName(usedNames=[]):
    name =''
    for i in range(0,9):
        choice = np.random.randint(0,len(abc))
        name+=abc[choice]
    if name not in usedNames:
        return name
    else: 
        return randomName(usedNames)
    
class chromosome:
    def __init__(self,name,genes:list = []):
        self.name = name
        self.score = 0
        self.scores = {}
        if genes:
            self.genes = self.make_genes(genes)
    def make_genes(self,genes):
        geneDict={}
        for gene in genes:
            geneDict[gene[:-1]] = gene[-1]
        
        return geneDict
        

def crossover(parent1: chromosome, parent2: chromosome,usedNames):
    genes1 = parent1.genes
    genes2 = parent2.genes
    child_genes = {}
    for key in genes1.keys():
        if np.random.rand() > 0.5:
            child_genes[key] = genes1[key]
        else:
            child_genes[key] = genes2[key]
            
    # Create and return the child chromosome
    name = randomName(usedNames)
    child_chromosome = chromosome(name)
    child_chromosome.genes=child_genes
    return child_chromosome

def mutation(child: chromosome,mutation_rate):
    for key in child.genes.keys():
        if np.random.rand()<mutation_rate:
            child.genes[key] = np.random.rand()
    return 




def defaultGENE(path,modes):
    genes = []
    df = pd.read_csv(path)
    for i in df.index:
        if df.at[i,"mode"] in modes:
            genes.append((df.at[i,"resonance"],df.at[i,"mode"],df.at[i,"branching ratio"]))
    
    
    df = pd.read_csv("../particle_ref.csv")
    
    for i in df.index:
        genes.append((df.at[i,"name"],"mass",df.at[i,"mass"]))
        genes.append((df.at[i,"name"],"width",df.at[i,"width"]))
    return genes


def genesFromFile(path):
    genes = []
    ref_part = pd.read_csv("../particle_ref.csv")
    ref_decay = pd.read_csv("../ref.csv")
    decaymodes = makeDF(path + "/decaymodes.txt")
    particles = makeDF_particles(path + "/particles.txt")
    for row in ref_decay.itertuples(index=True):
        resonance = row.resonance
        mode = row.mode
        upper = row.upper
        lower = row.lower
        value = decaymodes[(decaymodes.resonance == resonance) & (decaymodes["mode"] == mode)]["branching ratio"]
        
        if not value.empty:
            # Get the first "branching ratio" value
            value = value.iloc[0]
    
        value = (value - lower) / (upper - lower)  
        genes.append((resonance, mode, value))
    
    for row in ref_part.itertuples(index=True):
        name = row.name
        mass_low = row.mass_low
        mass_high = row.mass_high

        width_low = row.width_low
        width_high = row.width_high
        
        mass = float(particles[(particles.name == name)]["mass"].iloc[0])
        width = float(particles[(particles.name == name)]["width"].iloc[0])
        
        mass = (mass-mass_low)/(mass_high-mass_low)
        width = (width-width_low)/(width_high-width_low)
        genes.append((name,"mass",mass))
        genes.append((name,"width",width))
    return genes

def chromosFromNewest():
    gen = newestGeneration()
    paths = glob.glob(gen+"/*")
    chromos = [chromosome(path.split("/")[-1],genesFromFile(path)) for path in paths]
    scores = pd.read_csv("../generations/scores.csv")
    newest = scores[scores["generation"] == float(gen.split("-")[-1])]
    names = []
    for chromo in chromos:
        names.append(chromo.name)
        
    return chromos,names


def branchmassGENE(path,modes):
    genes = []
    df = pd.read_csv(path)
    for i in df.index:
        if df.at[i,"mode"] in modes:
            genes.append((df.at[i,"resonance"],df.at[i,"mode"],df.at[i,"branching ratio"]))
    
    
    df = pd.read_csv("../particle_ref.csv")
    
    for i in df.index:
        genes.append((df.at[i,"name"],"mass",df.at[i,"mass"]))
        #genes.append((df.at[i,"name"],"width",df.at[i,"width"]))
    return genes



def branchGENE(path,modes):
    genes = []
    df = pd.read_csv(path)
    for i in df.index:
        if df.at[i,"mode"] in modes:
            genes.append((df.at[i,"resonance"],df.at[i,"mode"],df.at[i,"branching ratio"]))
    return genes






def widthGENE(path):
    genes = []
    
    df = pd.read_csv("../particle_ref.csv")
    
    for i in df.index:
        genes.append((df.at[i,"name"],"width",df.at[i,"width"]))
    return genes

def massGENE(path):
    genes = []
    
    df = pd.read_csv("../particle_ref.csv")
    
    for i in df.index:
        genes.append((df.at[i,"name"],"mass",df.at[i,"mass"]))
    return genes





def renorm(df,info):
    df = df.copy()
    resonances = np.unique(info["resonance"].values)    
    for i,resonance in enumerate(resonances):
        temp = info[ info["resonance"]==resonance]
        modes = temp["mode"].values

        one= 1
        weight = (one - temp["new branching ratio"].sum()) / (one - temp["branching ratio"].sum()) 
        for j in df.index:
            if df.at[j,"resonance"] == resonance:
                if not df.at[j,"mode"] in modes:
                    df.at[j,"branching ratio"]*=weight
        
                        
                if df.at[j,"mode"] in modes:
    
                    
                    newBR = temp[temp["mode"]==df.at[j,"mode"]]["new branching ratio"].values[0]
               
                    df.at[j,"branching ratio"]=newBR
                
    return df



def makeINFO(ref,chromo):
    """
    Parameters
    ----------
    ref : TYPE
        DESCRIPTION.
    chromo : chromosome
        DESCRIPTION.

    Returns
    -------
    info : TYPE
        DESCRIPTION.

    """
    new_br = []
    info = ref.copy()
    for i in ref.index:
        res = ref.at[i,"resonance"]
        mode = ref.at[i,"mode"]
        key = (res,mode)
        new_br.append(chromo.genes[key])
    High = ref["upper"].fillna(0.03)
    Low = ref["lower"].fillna(0.001)
    new_br = np.array(new_br) *(High - Low) + Low
    info["new branching ratio"] = new_br
    return info


def makeINFO_particle(ref,chromo,case="default"):
    info = ref.copy()
    new_width = []
    new_mass = []
    for i in ref.index:    
        name = ref.at[i,"name"]
        mass_key = (name,"mass")
        width_key = (name,"width")
        if case!="mass":
            width = chromo.genes[width_key]
            new_width.append(width) 
        if case!="width":
            mass = chromo.genes[mass_key]
            new_mass.append(mass)
    
    if case!="mass":
        new_width = np.array(new_width)
        width_low,width_high = ref["width_low"],ref["width_high"]
        new_width = new_width*(width_high - width_low) + width_low
        info["new width"] = new_width
    
    if case!="width":
        new_mass = np.array(new_mass)
        mass_low,mass_high = ref["mass_low"],ref["mass_high"]
        new_mass = new_mass*(mass_high - mass_low) + mass_low
        info["new mass"] = new_mass  
         
    
    
    return info



def randomChromo(refChromo,name):
    rando = chromosome(name)
    dic = {}
    for key in refChromo.genes.keys():
        
        dic[key] = np.random.rand()
        #dic[key] = 0.0
    rando.genes=dic
    return rando


def makeDECAYFILE(ref,chromo,folder):
    new_folder = "{}/{}".format(folder,chromo.name)
    
    if os.path.exists(new_folder):
        pass
    else:
        os.makedirs(new_folder)

    new_file="{}/{}/{}.txt".format(folder,chromo.name,"decaymodes")
   

    

    FILE = decayFILE(orig_path)
    
    # This part is confusing put necesarry. 
    # Essentially, we have ensure that the value given from the chromosome
    # are not normalized by smash. Hence, we normalize them ourself.
    
    # Transform the branching ratios as : 
    # new_br = np.array(new_br) * (High - Low) + Low.
    info = makeINFO(ref,chromo)
    
    original = FILE.df  
    
    # Use the information from Info to renormalize the branching ratios 
    # before smash.
    newdf = renorm(original, info)
    # Make the new decaymodes file at the specifed place.
    FILE.df = newdf
    FILE.makeDecayBlocks(new_file)
    
    
    return newdf

def makePARTICLEFILE(ref,chromo,folder,case="default"):

    new_folder = "{}/{}".format(folder,chromo.name)
    #print(new_folder)
    if os.path.exists(new_folder):
        pass
    else:
        os.makedirs(new_folder)

    new_file="{}/{}/{}.txt".format(folder,chromo.name,"particles")

    info = makeINFO_particle(ref, chromo,case)
    original = makeDF_particles("../particles.txt")
    for i in info.index:
        name = info.at[i,"name"]
        if case!="mass":
            width = info.at[i,"new width"]
            original.loc[original["name"]==name, "width"] = width
        if case!="width":
            mass = info.at[i,"new mass"]
            original.loc[original["name"]==name, "mass"] = mass      
    makeParticleFILE(new_file, original)
    return info

def initPopulation(numberOfChromosomes,case="default"):
    chromos = []
    names = []
    if case == "default":
        ADAM = chromosome("adam",defaultGENE("../ref.csv",modes))
    elif case == "mass":
        ADAM = chromosome("adam",massGENE("../ref.csv"))
    elif case == "width":
        ADAM = chromosome("adan",widthGENE("../ref.csv"))
    elif case == "branch":
        ADAM = chromosome("adam",branchGENE("../ref.csv",modes))
    elif case == "branch+mass":
        ADAM = chromosome("adam",branchmassGENE("../ref.csv",modes))
    else:
        raise ValueError("Invalid case specified") 
                          
    for i in range(0,numberOfChromosomes):
        name = randomName(names)
        names.append(name)
        rando = randomChromo(ADAM,name)
        print(name)
        chromos.append(rando)
    
    return chromos,names




## Reads the reference csv and creates chromosomes whith "genes" as speficied in the reference csv.
def makeOrder(population,path,case="default"):
    ref = pd.read_csv("../ref.csv")
    pref = pd.read_csv("../particle_ref.csv")
    if case == "default":
        for chromo in population:
            makeDECAYFILE(ref, chromo,path)
            makePARTICLEFILE(pref, chromo, path)
    elif case=="mass" or case == "width":
        for chromo in population:
            makePARTICLEFILE(pref, chromo, path,case)
            folder="{}/{}/".format(path,chromo.name)
            shutil.copy(orig_path, folder)
            
    elif case=="branch":
        for chromo in population:
            makeDECAYFILE(ref, chromo,path)
            folder="{}/{}/".format(path,chromo.name)
            shutil.copy("../particles.txt", folder)
    elif case == "branch+mass":
        for chromo in population:
            makeDECAYFILE(ref, chromo,path)
            makePARTICLEFILE(pref, chromo, path,"mass")
            folder="{}/{}/".format(path,chromo.name)
            #shutil.copy(orig_path, folder)
    else:
        raise ValueError("Invalid case specified") 
        
    
    return
def saveParents(parents,path):
    
    return 

## Indeeded to be used incase the script ened prematurely. Also used for debugging. 
## Reads the decaymode files in the specifed path and returns a list of chromosomes.
def orderToChromosomes(order_path):
    ref = pd.read_csv("../ref.csv")
    
    chromosomes_paths = glob.glob(order_path+"/*.txt")

    ADAM = chromosome("adam",defaultGENE("../ref.csv",modes))
    info = makeINFO(ref,ADAM)
    population = []
    for i,chromo_path in enumerate(chromosomes_paths):
       
        name = chromo_path.split("/")[-1].split(".")[0]
        FILE = decayFILE(chromo_path)
        df = FILE.df
        rando = chromosome(name)
        dic = {}
        for key in ADAM.genes.keys():
            temp = ref[ref["resonance"]==key[0]]
            High = temp[temp["mode"] == key[1]]["upper"].values[0]
            Low = temp[temp["mode"] == key[1]]["lower"].values[0]
            if High != High or Low != Low:
                raise Exception("Bounds in referance file are missing")

  
            
            index = df[(df["resonance"] == key[0]) & (df["mode"] == key[1])].index[0]
            BR = df.at[index,"branching ratio"]
            q = (BR-Low)/(High-Low)
            dic[key] = q
        rando.genes=dic
        population.append(rando)
        
            
    return population


def selection(shavePercent,population):
    sorted_pop = sorted(population, key=lambda x: x.score)
    selec = len(sorted_pop) - int(shavePercent*len(sorted_pop))
    #print(selec)
    return sorted_pop[selec:]
    #return sorted_pop





def newestGeneration():
  
    paths=glob.glob("../generations/generation-*")
    sorted_paths=[]
    dic={}
    for path in paths:
        dic[path]=int(path.split("-")[-1])

    dic=sorted(dic.items(), key=lambda x: x[1])
    
    for key in dic:
        sorted_paths.append(key[0])


    return sorted_paths[-1]

def nextGeneration():
    numb = 0
    try:
        newest = newestGeneration()
        numb = int(newest.split("-")[-1])
    except:
        numb = -1

    nextGen = "../generations/generation-{}".format(numb+1)
    
    return nextGen


def nextPopulation(gen_path,populationSize,shave,population,mutation_rate=0.1):
    ### Selection and creation of new chromosomes ###
    new_pop = selection(shave,population)
    names = [chromo.name for chromo in new_pop]
    children = []

    while len(new_pop) + len(children) < populationSize:
        sample=random.sample(new_pop,2)
        child = crossover(sample[0], sample[1], names)
        mutation(child,mutation_rate)
        children.append(child)
        names.append(child.name)
    return new_pop,children,


def runOrder(gen_path,threads,population,rot_s):
    os.makedirs("{}/order".format(gen_path))
    makeOrder(population,gen_path)
    command = "../shellScripts/order.sh {} {}".format(threads,gen_path)

    for energy in rot_s:
        command+=" {}".format(energy)
    print(subprocess.run([command], shell=True))

    return

def lastScore():
    population= orderToChromosomes(newestGeneration()+"/order")
    for chromo in population:
        chromo.score = chromoScore(gen_path+"/{}".format(chromo.name))
        
    score_board = [(chromo.name,chromo.score) for chromo in population]

    with open(gen_path+"/scores.csv","w") as f:
        f.write("name,score\n")
        for i in score_board:
            f.write("{},{}\n".format(i[0],i[1]))
    return





def rot_s_to_plab(rot_s, m):
    s = rot_s
    plab = (((s - 4 * m**2)/(2 * m))**0.5)
    return plab
