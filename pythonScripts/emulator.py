import os
import pandas as pd
import numpy as np
from experiment_data_loader import ExperimentalDataLoader
import concurrent.futures
import os
import glob
import sys
import matplotlib.pyplot as plt
import subprocess
import math
import numpy as np
from scipy.interpolate import UnivariateSpline


class Emulator_fast:
    def __init__(self, emulator_path,cross_dir):
        self.reactions_pp = ["K⁺+p+Λ", "K⁰+p+Σ⁺", "K⁺+n+Σ⁺", "K̅⁰+K⁰+p+p","K⁺+p+Σ⁰"]
        self.reactions_piminus_proton = ["K⁺+Σ⁻","K⁰+Λ"]
        self.reactions_piplus_proton = ["K⁺+Σ⁺"]
        self.reactions_pions = self.reactions_piminus_proton + self.reactions_piplus_proton
        self.emulator_path = emulator_path
        self.cross_dir = cross_dir
        self.loader = ExperimentalDataLoader(cross_dir)
        self.exp = self.loader.data 
        self.string_xsections = {}
        self.string_energies = {}
        self.extract_string_xsections()

        
    def process_chromo(self,chromo, gen_path):
       
        path = f'{gen_path}/{chromo.name}'
        self.emulate_chromo(chromo,path)
       
        return 
    def process_chromos(self, chromosomes, gen_path,max_workers=1):
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_chromo = {executor.submit(self.process_chromo, chromo, gen_path): chromo for chromo in chromosomes}
            for future in concurrent.futures.as_completed(future_to_chromo):
                chromo = future_to_chromo[future]
                try:
                    result_chromo = future.result()
                # Handle the result if needed
                except Exception as exc:
                    print(f'Chromosome processing failed: {exc}')
        return chromosomes
    
    def emulate_chromo(self,chromo,chromo_path):
        decaymodes_path = chromo_path+"/decaymodes.txt"
        particles_path = chromo_path+"/particles.txt"
        arguments = [chromo_path,particles_path,decaymodes_path]
        try:
            result = subprocess.run(
                [self.emulator_path + "/build/run"] + arguments,
                check=True,
                capture_output=True,
                text=True
                )
        except subprocess.CalledProcessError as e:
            print("Error:", e)
            print("Standard Output:", e.stdout)
            print("Standard Error:", e.stderr)  # This will give you the detailed error message
        
        # Get all CSV files in the directory
        data_files = glob.glob(chromo_path + "/*.csv")

        
        # Initialize an empty dictionary to store numpy arrays
        data_dict = {}
        reactions = []
        # Read each CSV file into a numpy array and store in dictionary
        for file in data_files:
            filename = file[len(chromo_path)+1:-4]  # Extract filename without .csv extension
            reactions.append(filename)
            #data = np.loadtxt(file, delimiter=',', dtype=float)  # Read data into numpy array
            data = pd.read_csv(file)
            #energies = self.exp[filename]["ssqrt"].values
            energies = data["ssqrt"]

            strings = self.getStringXsection(filename, energies)
   
            data["sigma"]+=strings
            data_dict[filename] = data 
        
        scores = {reaction: self.score(data_dict[reaction], reaction) for reaction in reactions}
        score = sum(scores.values())/len(scores)
        
        #print(scores)
        
        if score == 0:
            score = np.inf
        else:
            score = 1/score
        chromo.score = score
        chromo.scores = scores
        return 
    def calculate_weights(self,reaction, sqrt_s):
        if reaction in self.reactions_piminus_proton+self.reactions_piplus_proton:
            weights = (1 - self.probability_transit_high(sqrt_s, 1.9, 2.1))
        else:  
            weights = (1 - self.probability_transit_high(sqrt_s, 3.5, 4.5))
    
        return weights
    def score(self, df, reaction):
        scores = []
        sig_xps = self.exp[reaction]["sigma"].values
        sqrt_s = self.exp[reaction]["ssqrt"].values
        sig_errs = self.exp[reaction]["yerror"].values
        weights = 0
        for index, row in df.iterrows():
            sig = row['sigma']
            energy = row['ssqrt']

            # Find the closest experimental energy value
            i = np.argmin(np.abs(sqrt_s - energy))
            sig_exp = sig_xps[i]
            sig_err = sig_errs[i]
            distance = np.abs(sig-sig_exp)
            weight = self.calculate_weights(reaction, energy)
            weights+=weight
            score = distance/sig_exp
            score = score*weight
            if(reaction in self.reactions_pions or reaction == "K⁺+p+Λ"):
                score*=20
                if(reaction == "K⁺+Σ⁻"):
                    score*=3
                if(reaction == "K⁺+p+Λ"):
                    score*=40
                    
            
            scores.append(score)

        # Compute the final score as the average of the individual scores
        return sum(scores)/weights

    
    def extract_string_xsections(self):

       
       piminus_proton = self.getData("default", "piminus_proton")
       piplus_proton = self.getData("default", "piplus_proton")
       proton_proton = self.getData("default", "proton_proton")
       
       # Initialize the dictionary
       self.string_xsections = {}

       # Spline fitting for piplus_proton
       sqrt_s_piplus = piplus_proton["sqrt_s"].values
       for reaction in self.reactions_piplus_proton:
           cross_section = self.probability_transit_high(sqrt_s_piplus, 1.9, 2.1) * piplus_proton[reaction]
           spline = UnivariateSpline(sqrt_s_piplus, cross_section, s=0.01)
           self.string_xsections[reaction] = spline
       
       # Spline fitting for piminus_proton
       sqrt_s_piminus = piminus_proton["sqrt_s"].values
       for reaction in self.reactions_piminus_proton:
           cross_section = self.probability_transit_high(sqrt_s_piminus, 1.9, 2.1) * piminus_proton[reaction]
           spline = UnivariateSpline(sqrt_s_piminus, cross_section, s=0)
           self.string_xsections[reaction] = spline
       
       # Spline fitting for proton_proton
       sqrt_s_proton = proton_proton["sqrt_s"].values
       for reaction in self.reactions_pp:
           cross_section = self.probability_transit_high(sqrt_s_proton, 3.5, 4.5) * proton_proton[reaction]
           spline = UnivariateSpline(sqrt_s_proton, cross_section, s=0)
           self.string_xsections[reaction] = spline

       return

    def getStringXsection(self, reaction, sqrt_s_values):
       if reaction not in self.string_xsections:
           raise ValueError(f"No spline fit found for reaction: {reaction}")
   
       # Get the spline
       spline = self.string_xsections[reaction]
   
       # Evaluate the spline at the given sqrt_s values
       res = spline(sqrt_s_values)
       res[res < 0] = 0
       return res
   
    
   
    def getData(self,folder,coll):

        
        path_template = "../saves/{}/{}/*/xs_final_individual.dat".format(folder,coll)
        
        if coll == "proton_proton":
            reactions = self.reactions_pp
        elif coll == "piminus_proton":
            reactions = self.reactions_piminus_proton
        elif coll == "piplus_proton":
            reactions = self.reactions_piplus_proton
            
        paths = glob.glob(path_template)
        
        
        data = {}
        data["sqrt_s"] = []
        
        for reaction in reactions:
            data[reaction] = []
            data[reaction +"_error"] = []
            
        for path in paths : 
            with open(path,"r") as file:
                lines = file.readlines()
                
                columns = np.array(lines[3].split())[1:]
                if(len(lines) !=5):
                    #print(path)
                    continue
    
                values = lines[4].split()
                data["sqrt_s"].append(values[0])
    
    
                for reaction in reactions:
                    position = np.where(columns==reaction)[0]
                    if position.size == 1:
                        if(columns[position[0]] !=reaction):
                            print("wtf")
                        data[reaction].append(values[position[0]])
                        data[reaction +"_error"].append(values[position[0]+1])
    
                    else:
                          data[reaction].append(0)
                          data[reaction +"_error"].append(0)
                          
        for key in data.keys():
            data[key] = np.array(data[key],dtype=float)
        df = pd.DataFrame(data)
        df = df.sort_values(by="sqrt_s")
        return df
    def probability_transit_high(self,sqrt_s_, region_lower, region_upper):
        # Ensure sqrt_s_ is a numpy array
        sqrt_s_ = np.asarray(sqrt_s_)

        # Initialize the probability array with zeros
        prob = np.zeros_like(sqrt_s_)

        # Condition 1: sqrt_s_ < region_lower
        mask1 = sqrt_s_ <= region_lower
        prob[mask1] = 0.0

        # Condition 2: sqrt_s_ > region_upper
        mask2 = sqrt_s_ >= region_upper
        prob[mask2] = 1.0

        # Condition 3: region_lower <= sqrt_s_ <= region_upper
        mask3 = ~mask1 & ~mask2
        x = (sqrt_s_[mask3] - 0.5 * (region_lower + region_upper)) / (region_upper - region_lower)
        x = np.clip(x, -0.5, 0.5)  # Ensure x is within [-0.5, 0.5]
        prob[mask3] = 0.5 * (np.sin(np.pi * x) + 1.0)

        # Ensure probabilities are within [0, 1]
        prob = np.clip(prob, 0.0, 1.0)

        return prob

