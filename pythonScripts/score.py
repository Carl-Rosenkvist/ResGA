#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 16:46:38 2023

@author: carl
"""

import numpy as np
import pandas as pd 
import yaml
import matplotlib.pyplot as plt
import glob
import os
import subprocess


cross_dir = "/home/carl/PhD/smash-analysis-extras/experimental_data/cross_sections"
pp_cs_path=glob.glob(cross_dir+"/proton_proton/*.exp")

file_path= cross_dir+'/proton_proton/lambda_proton_kplus.exp'


def get_super(x): 
    normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
    super_s = "ᴬᴮᶜᴰᴱᶠᴳᴴᴵᴶᴷᴸᴹᴺᴼᴾQᴿˢᵀᵁⱽᵂˣʸᶻᵃᵇᶜᵈᵉᶠᵍʰᶦʲᵏˡᵐⁿᵒᵖ۹ʳˢᵗᵘᵛʷˣʸᶻ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾"
    res = x.maketrans(''.join(normal), ''.join(super_s)) 
    return x.translate(res) 

def load_lines(file_path):
    
    with open(file_path) as f:
        lines = f.readlines()
    lines = [i.split() for i in lines]
    return lines


def load_exp(exp_path):
    
    with open(exp_path) as f:
        lines = f.readlines()
    lines = np.array([i.split() for i in lines if i[0]!="#"])
    p_lab = lines[1:,0].astype(float)
    sigma = lines[1:,1].astype(float)
    yerror = lines[1:,2].astype(float)
    
    
    mt=mp=0.938270
    ssqrt=(((p_lab**2 + mp**2)**0.5 + mt)**2 - p_lab**2)**0.5
    dic = {}
    dic["ssqrt"]=ssqrt
    dic["sigma"]=sigma
    dic["yerror"]=yerror
    
    return pd.DataFrame(dic)


def find_modes(modes,columns):
    
    dic = {}
    for mode in modes:
        listed_mode = []
        for count,char in enumerate(mode):
            if char == '⁺' or char == '⁻':
                listed_mode[count-1] += get_super(char)
            else:
                listed_mode.append(char)
                
        for count,column in enumerate(columns):
            is_in = True
            for particle in listed_mode:
                is_in *= particle in column
            if is_in:
                dic[column] = count
            
                
                
    return dic
def load_energy(config_path):
    with open(config_path, "r") as stream:
        try:
            return yaml.safe_load(stream)["Modi"]["Collider"]["Sqrtsnn"]
        except yaml.YAMLError as exc:
            print(exc)

    
    return


def load_analysis(analysis_path):
    dat="/xs_final_individual.dat"
    dat_res="/xs_final_individual_res.dat"
    lines = load_lines(analysis_path+dat)
    #print(lines)
    columns = np.array(lines[3])[2:]
    vals = np.array(lines[4]).astype(float)[1:]
    lines_res = load_lines(analysis_path+dat_res)
    columns_res = np.array(lines_res[3])[2:]
    vals_res = np.array(lines_res[4]).astype(float)[1:]
    dic = {}
    dic["modes"] = np.append(columns,columns_res)
    dic["vals"] = np.append(vals,vals_res)
    
    return pd.DataFrame(dic)


def findCross(path):
  
    mode = "K⁺+p+Λ"
    analysis = load_analysis(path+"/analysisFiles")
    try:
        sigma = analysis[analysis["modes"]==mode]["vals"].values[0]
    except:
        sigma = 0.0
    energy = load_energy(path+"/data/config.yaml")
    return sigma,energy

def CalcXsections(chromo_path,threads):
    #chromo_path = "../generations/generation-0/NFsUZE4G0/"
    output_path = chromo_path+"/analysis"
    
    paths = glob.glob(chromo_path+'/sqrt-*/data/collisions_binary.bin')
    command = "./../CalcXsection/build/CalcXsection {} {}".format(threads,output_path)
    
    for path in paths:
        command+= " "+path
    #print(command)
    subprocess.run([command], shell=True)
    
    return 
    
    
    
    
    


def score(sig,energy,exp):
    
  
    
    energy_xp = exp["ssqrt"].values
    sig_xp = exp["sigma"].values
    sig_err = exp["yerror"].values
    
    indx = np.argmin(np.abs(energy_xp - energy))
    
    sig_exp = sig_xp[indx]
    sig_err = sig_err[indx]
    chi_h = ((sig_exp+sig_err - sig)**2)/(sig_exp+sig_err)
    chi_m = ((sig_exp - sig)**2)/sig_exp
    chi_l = ((sig_exp-sig_err - sig)**2)/(sig_exp-sig_err)
    return np.min([chi_h,chi_m,chi_l])
    
def chromoScore(chromo_path):
    exp=load_exp(file_path)
    scores = []
    output_path = chromo_path+"/analysis/output.csv"
    result = pd.read_csv(output_path)
    sigmas = result["Sigma"].values
    energies = result["Energy"].values
    for i in range(0,len(energies)):
        
        scores.append(score(sigmas[i],energies[i],exp))
        
        
    if sum(scores) != 0:
        return (sum(scores)/len(scores))**(-1)
    else: 
        return float("inf")

    


