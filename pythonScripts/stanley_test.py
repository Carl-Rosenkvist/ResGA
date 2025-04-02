#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 22:20:30 2024

@author: carl
"""
import pandas as pd
import numpy as np
from genetic import *
from decayFILE import *
path = "../stanley_test"
chromo = chromosome("small",massGENE("../ref.csv"))
ref = pd.read_csv("../particle_ref.csv")

for gene in chromo.genes.keys():
    chromo.genes[gene] = 0

makePARTICLEFILE(ref,chromo,path,case="mass")
