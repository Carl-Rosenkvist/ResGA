#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 10:27:40 2023

@author: Carl Rosenkvist
"""

import numpy as np
import pandas as pd
#import pdg
#import ast
#api = pdg.connect()

def getLines(path):

    lines = ""
    with open(path) as f:
        start = False
        for line in f:
            if line[0]!="#":
                lines+=line
            
    return lines
      
def blocks(path):
    lines = getLines(path).split("\n")
    blocks = []
    block = []
    for line in lines:
        if line != '':
            block.append(line)
        if line == '':
            if block:
                blocks.append(block)
            block = []
    return blocks
            
def makeDF(path):
    b = blocks(path)
    res = []
    l = []
    modes = []
    br = []
    for block in b:
        
       
        for line in block[1:]:
            res.append(block[0].split()[0])
           
            #res.append(block[0])
            br.append(float(line.split()[0]))
            l.append(line.split()[1])
            
            mode = line.split()[2] +' '+line.split()[3]
            try:
                add_mode = line.split()[4]
            except:
                add_mode = "#"
            if add_mode !="#":
                mode+=' '+add_mode
            modes.append(mode)
    dic = {}
    dic["resonance"] = np.array(res)
    dic["branching ratio"] = np.array(br)
    dic["l"] = np.array(l)
    dic["mode"] = np.array(modes)
    return pd.DataFrame(dic)
        


class decayFILE:
    def __init__(self,path=None,df=None):
        if df is None:
            if path is None:
                raise Exception("No path nor Dataframe given")
            else:
                self.df = makeDF(path)
        else:
            self.df = df
    
    def makeDecayBlocks(self,name):
        resonances = self.df["resonance"].values
        modes = self.df["mode"].values
        l = self.df["l"].values
        br =self.df["branching ratio"].values
        
       

        res_uniq = []
        for res in resonances:
            if res not in res_uniq:
                res_uniq.append(res)
        
        #res_uniq =resonances
        with open("{}".format(name),"w") as f:
            for res in res_uniq:
                mask = np.where(resonances == res)[0]
                #print(res)
                f.write(res+"\n")
                for i in mask:

                    br_str = '{:f}'.format(br[i])
                    #br_str = br[i]
               
                    new_line = "{0: <6} {1: <2} {2: <10}".format(br_str,l[i],modes[i])
                    f.write(new_line)
                    f.write("\n")
                f.write("\n")
                    
        return
    def normalizeDF(self):
        resonances = np.unique(self.df["resonance"].values)
        for res in resonances:
            temp = self.df[self.df["resonance"]==res]
            brsum = temp["branching ratio"].sum()
            for i in self.df.index:
                if self.df.at[i,"resonance"] == res:
                    self.df.at[i,"branching ratio"]/=brsum
        
        
        return




def makeDF_particles(path):
    lines = getLines(path).split("\n")
    names = []
    masses = []
    widths = []

    parities = []
    pdgs = []
    for line in lines:
        if line=='':
            continue
        
        split = line.split()
        names.append(split[0])
        masses.append(split[1])
        widths.append(split[2])
        parities.append(split[3])
        temp = []
        for pdg in split[4:]:
            if pdg[0]=="#":
                break
            else:
                temp.append(pdg)
        pdgs.append(temp)
    dic = {}
    dic["name"] = np.array(names)
    dic["mass"] = np.array(masses)
    dic["width"] = np.array(widths)
    dic["parity"] = np.array(parities)
    dic["pdg(s)"] = pdgs
    return pd.DataFrame(dic)
        
def getProperties(name):
    mass_bounds = None
    width_bounds = None
    if(name[0] == "Δ"):
        name="Delta"+name[1:]
    
    if(name!="N(1900)"):
        particle=api.get_particle_by_name(name)
    else:
        particle = api.get('B144')
    for p in particle.properties():
        if p.description.split()[1:] == ['BREIT-WIGNER', 'MASS']:
            mass_lower=float(p.display_value_text.split()[0])/1000
            mass_upper = float(p.display_value_text.split()[-1])/1000
            mass_bounds=(mass_lower,mass_upper)
        if p.description.split()[1:] == ['BREIT-WIGNER', 'WIDTH']:
            width_lower=float(p.display_value_text.split()[0])/1000
            width__upper=float(p.display_value_text.split()[-1])/1000
            width_bounds = (width_lower,width__upper)
    return mass_bounds,width_bounds

def getResonances():
    df = makeDF_particles('../particles.txt')
    res_df = df[df['name'].str.contains('[NΔ]\(\d+\)')]
    res_df[['mass_bounds', 'width_bounds']] = res_df['name'].apply(lambda x: pd.Series(getProperties(x)))
    

    # Split 'mass_bounds' into 'mass_low' and 'mass_high'
    res_df[['mass_low', 'mass_high']] = pd.DataFrame(res_df['mass_bounds'].tolist(), index=res_df.index)
    
    # Split 'width_bounds' into 'width_low' and 'width_high'
    res_df[['width_low', 'width_high']] = pd.DataFrame(res_df['width_bounds'].tolist(), index=res_df.index)
    
    # Drop the original columns if needed
    res_df = res_df.drop(['mass_bounds', 'width_bounds'], axis=1)
    return res_df

def createParticleRef():
    df = getResonances()
    
    df.to_csv("../particle_ref.csv",index=False)
    
    return
def getParticleRef():
    df = pd.read_csv("../particle_ref.csv")
    #df["mass_bounds"] = df["mass_bounds"].apply(ast.literal_eval)
    #df["width_bounds"] = df["width_bounds"].apply(ast.literal_eval)
    # Assuming df is your DataFrame

    
    return df

def makeParticleFILE(path, df):
    with open(path, 'w') as output_file:
        output_file.write("# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n\n")
        for i in df.index:
            name = df.at[i, "name"]
        
            # Convert the mass to a float (if it's a string)
            mass = float(df.at[i, "mass"])
            
            # Convert the width to a float (if it's a string)
            width = float(df.at[i, "width"])
            
            parity = df.at[i, "parity"]
            
            # Convert each element in the pdg(s) list to a string
            pdgs = df.at[i, "pdg(s)"]
            
            pdgs = " ".join(["{:<7}".format(str(pdg)) for pdg in pdgs])
            # Write data to the file with fixed-width columns
            new_line = "{:<15} {:<6.3f} {:<12.3e} {:<6} {}".format(name, mass, width, parity,pdgs)
    
            
            output_file.write(new_line)
            
            output_file.write("\n")
        return 
        
        

class ParticlesFILE:
    def __init__(self,path=None,df=None):
        if df==None:
            if path == None:
                raise Exception("No path nor Dataframe given")
            else:
                self.df = makeDF_particles(path)
        else:
            self.df = df    

