import os
import glob
import sys
from timeit import default_timer as timer
import numpy as np
import csv
import math
import pandas as pd

n_cores = 5

#For WT + random metrics

def run_net_metrics(net):
    topo_files = glob.glob("/mnt/e/varun/node_deletions/" + net + "/*.topo")
    files = []
    for file in topo_files:
        files.append(file.split('/')[6].split('.')[0])
    fields = ["Net","maxFrust", "minFrust", "meanFrust", "maxNetFrust", "minNetFrust", "meanNetFrust", "maxCoh", 
              "minCoh", "meanCoh", "maxNetCoh", "minNetCoh", "meanNetCoh", "corFreqFrust", "pFreqFrust", "corFreqCoh", "pFreqCoh", "corFrustCoh", "pFrustCoh"]
    with open('/mnt/e/varun/node_deletions/' + net + '_ALL.csv', 'a') as file:
        writer = csv.writer(file)
        writer.writerow(fields)
    end = math.ceil(len(files)/n_cores)
    for i in range(0,end):
        num_files_iter = np.arange(i*n_cores, (i*n_cores) + n_cores, 1)
        to_files = files[i*n_cores:(i*n_cores) + n_cores]
        arg = " ".join([str(i) for i in to_files])
        command = "mpirun -n " + str(n_cores) +" python3 -Wignore  /mnt/e/varun/coh_vs_p_n.py " + arg 
        print(command)
        os.system(command)  

nets = ["drosophila"]
[run_net_metrics(i) for i in nets]

#For n node pert hamming, coherence,score (WT)


'''
def run_multi_node_metrics(net):
    states_df = pd.read_csv("/mnt/e/varun/Results/Composition/EMT_MET_comp.csv")
    states = states_df['States'].tolist()
    field = [''] + states
    with open('/mnt/e/varun/try/EMT_MET_hamming.csv', 'a') as file:
        writer = csv.writer(file)
        writer.writerow(field)
    with open('/mnt/e/varun/try/EMT_MET_coherence.csv', 'a') as file:
        writer = csv.writer(file)
        writer.writerow(field)
    with open('/mnt/e/varun/try/EMT_MET_score.csv', 'a') as file:
        writer = csv.writer(file)
        writer.writerow(field)
    num_nodes = len(states[0]) - 2
    print(num_nodes)
    end = math.ceil(num_nodes/n_cores)
    for k in range(1):
        for i in range(0,end):
            num_n_iter = np.arange(i*n_cores, (i*n_cores) + n_cores, 1)
            arg = " ".join([str(i) for i in num_n_iter])
            command = "mpirun -n " + str(n_cores) +" python3 -W ignore n_node_pert_hamming_coh.py " + arg 
            print(command)
            os.system(command)  
'''
