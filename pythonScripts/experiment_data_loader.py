import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
def remove_comment(data):
  

    if "#" in data:
        i = data.index("#")
        return data[:i]
    else:
        return data
    
   


class ExperimentalDataLoader:
    def __init__(self, cross_dir):
        self.cross_dir = cross_dir
        self.pp_exps = {
            "K⁺+p+Λ": "proton_proton/lambda_proton_kplus.exp",
            "K⁰+p+Σ⁺": "proton_proton/sigmaplus_proton_kzero.exp",
            "K⁺+n+Σ⁺": "proton_proton/sigmaplus_neutron_kplus.exp",
            "K̅⁰+K⁰+p+p" : "proton_proton/proton_proton_kzero_kzero.exp",
            "K⁺+Σ⁺" : "piplus_proton/sigmaplus_kplus.exp",
            "K⁺+Σ⁻" : "piminus_proton/sigmaminus_kplus.exp",
            "K⁰+Λ" : "piminus_proton/lambda_kzero.exp",
            "K⁺+p+Σ⁰" :"proton_proton/sigmazero_proton_kplus.exp",
           }
        
        # Load data during initialization
        self.data = {key: self.load_data(key) for key in self.pp_exps}

    def load_exp(self, exp_path,reaction_key):
        with open(exp_path) as f:
            lines = f.readlines()

        lines =[i.split() for i in lines if i[0] != "#"]
        lines = [remove_comment(line) for line in lines if line]
      
        lines = np.array(lines)
        p_lab = lines[1:, 0].astype(float)
        sigma = lines[1:, 1].astype(float)
        yerror = lines[1:, 2].astype(float)
   
        if(reaction_key=="K⁺+Σ⁺" or reaction_key == "K⁺+Σ⁻" or reaction_key == "K⁰+Λ"):

            mp = 0.140
            mt = 0.938270
        else:
            mp = 0.938270
            mt = 0.938270
            
        
        ssqrt = (((p_lab**2 + mp**2)**0.5 + mt)**2 - p_lab**2)**0.5
        dic = {"ssqrt": ssqrt, "sigma": sigma, "yerror": yerror}
        return pd.DataFrame(dic)

    def load_data(self, reaction_key):
        exp_path = f"{self.cross_dir}/{self.pp_exps[reaction_key]}"
        return self.load_exp(exp_path,reaction_key)
    
    def plot_data(self):
        for key, df in self.data.items():
            plt.figure()
            plt.title("$pp \\rightarrow {} $".format(key))
            plt.errorbar(df["ssqrt"], df["sigma"], yerr=df["yerror"], fmt='o', label=key)
            plt.xlabel('$\sqrt{s}\ [\mathrm{GeV}]$',size=14)
            plt.ylabel("$\sigma\ [\mathrm{mb}]$",size=14)
            plt.legend()
            plt.show()

