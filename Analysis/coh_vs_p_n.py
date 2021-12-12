import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numba import jit
import os
from mpi4py import MPI
import sys
import networkx as nx
import csv
from scipy.stats import pearsonr

comm = MPI.COMM_WORLD # this is a communicators
rank = comm.Get_rank() # this tells you the id of the current instance

num_id = 1 + rank
file = str(sys.argv[num_id]) 
print("my rank is " + str(rank) + " and my file is " + file)
net = file.partition('_del')[0]

def compute_hamming(rand_net_file, wt_net_file, column):
    rand_net = pd.read_table(rand_net_file, sep = " ")
    wt_net = pd.read_table(wt_net_file, sep = " ")
    hamming = np.sum([0 if (rand_net.iloc[i,column] == wt_net.iloc[i,column]) else 1 for i in range(rand_net.shape[0])])
    return hamming

def topo2interaction(file, rule = 0, weight=100):
    table = pd.read_table(file, sep=" ")
    tmp = table[["Source", "Target"]].values.reshape(-1)
    node_labels = sorted(list(set(tmp)))
    N = len(node_labels)
    idx2label = dict(enumerate(node_labels))
    label2idx = {v: k for k, v in idx2label.items()}
    T2J_dict = {1: 1, 2: -1}
    if rule != 0:
        T2J_dict[rule] = T2J_dict[rule]*weight
    J = np.zeros((N, N), np.float64)
    for u, v, t in table.values:
        j = label2idx[u]
        i = label2idx[v]
        J[i, j] = T2J_dict[t]
    return J, np.array(node_labels)

@jit(nopython=True)
def _run(J, J_pseudo,initial_pert_ss,maxT,mode, can_be_updated):
    s = initial_pert_ss.copy()
    s_check = np.sign(J_pseudo@s)
    if (np.all(s_check==s)):
        convergence = True
        return convergence, s
    convergence = False
    CT = 0
    for ct in range(maxT):
        n = ct//100
        if mode == "async":
            k = np.random.choice(can_be_updated)
            sk_new = np.sign((J_pseudo@s)[k])
            if s[k] != sk_new:
                s[k] = sk_new
        else:
            s = np.sign(J_pseudo@s)
        s_check = np.sign(J_pseudo@s)
        if (np.all(s_check==s)):
            CT = n
            convergence = True
            return convergence, s
    return convergence, s

def state_run(net_file,states_file,rule,mode,maxT,lang):
    state_df = pd.read_csv(states_file)
    try:
        flag_state_df = state_df[state_df['flag'] == 1]
        flag_state_df = flag_state_df[flag_state_df['Avg0'].notna()]
        states = flag_state_df['states'].tolist()
        freq = flag_state_df['Avg0'].tolist()
    except:
        flag_state_df = state_df
        flag_state_df = flag_state_df[flag_state_df['async_0_Mean'].notna()]
        states = flag_state_df['State'].tolist()
        freq = flag_state_df['async_0_Mean'].tolist()
    J,node_labels = topo2interaction(net_file,rule=rule)
    num_nodes, M = J.shape
    can_be_updated = np.array([a for a, b in enumerate(node_labels)])
    J_pseudo = np.identity(num_nodes) + 2 * J
    w=0
    coherence_state = []
    for state in states:
        #print(w)
        state_s = state[1:len(state)-1]
        state_n = []
        for s in state_s:
            state_n.append(int(s))
        state_n = np.array(state_n, dtype=np.float64)
        indices_one = state_n == 1
        indices_zero = state_n == 0
        state_n[indices_one] = 1 
        state_n[indices_zero] = -1
        steady_state = state_n
        coherence_node = 0
        for i in range(len(state_n)):
            state_pert = state_n.copy()
            state_pert[i] = state_pert[i]*(-1)
            for j in range(100):
                convergence,s = _run(J, J_pseudo, state_pert,maxT, mode, can_be_updated)
                if convergence:
                    if np.all(s == steady_state):
                        coherence_node = coherence_node + 1
        coherence_state.append(coherence_node/(len(state_n)*100))
    return coherence_state,freq


def networkNcycles(name):
    G = nx.DiGraph()
    if '.topo' in name:
        f = open(str(name),'r')
        content = f.read()
        f.close()
        content = content.split('\n')
        content.remove(content[0])
        content = [content[i].split() for i in range(len(content))]
        cdict = {}
        content = [x for x in content if not x == []]
        for i in range(len(content)):
            G.add_edge(content[i][0], content[i][1])
            if content[i][2] == '2':
                cdict[str(content[i][0] + ','+content[i][1])] = -1
            elif content[i][2] == '1': 
                cdict[str(content[i][0] + ','+content[i][1])] = 1
        #print(cdict)
        cycles= list(nx.simple_cycles(G))
        return [G, cycles, cdict]
        # counting the number of positive and negative loop
    else:
        raise Exception("Input should be a topo file")


def loopNEdgeCounter(cycles, cdict):
    p_fbc = 0
    n_fbc = 0
    for i in range(len(cycles)):
        prod = 1
        for j in range(len(cycles[i])):
            prod = prod * cdict[cycles[i][j]+','+cycles[i][(j+1)%len(cycles[i])]]
        if prod == -1:
            print(cycles[i])
            p_fbc = p_fbc+1
        else:
            n_fbc = n_fbc + 1
    return [p_fbc, n_fbc]

def edgeCounter(cycles, cdict, G):
    edict = {'Cycles': [], 
    'Nature':[], 
    'Edge_count':[]};
    cyc = [",".join(x) for x in cycles]
    edict["Cycles"] = cyc
    for i in range(len(cycles)):
        prod = 1
        for j in range(len(cycles[i])):
            prod = prod * cdict[cycles[i][j]+','+cycles[i][(j+1)%len(cycles[i])]]
        if prod ==1:
            #print(cycles[i])
            edict["Nature"].append("P")
        else:
            edict["Nature"].append("N")
        edict["Edge_count"].append(len(cycles[i]))
        #print(edict)
    return(pd.DataFrame(edict))


def calc_loops(topofile):
    list1 = networkNcycles(topofile)
    cycles = list1[1]
    G = list1[0]
    cdict = list1[2]
    edict = edgeCounter(cycles, cdict, G)
    return edict

def calc_weighted_loops(loops_df, nature):
    n_df = loops_df[loops_df['Nature'] == nature]
    total = n_df.shape[0]
    li = np.array(n_df['Edge_count'].tolist())
    w_loops = np.sum(1/li)/total
    return w_loops

def get_frust(states_file):
    if '_finFlagFreq' in states_file:
        state_df = pd.read_csv(states_file)
        flag_state_df = state_df[state_df['flag'] == 1]
        flag_state_df = flag_state_df[flag_state_df['Avg0'].notna()]
        frust = flag_state_df['frust0'].tolist()
    return frust


topo_fl = "/mnt/e/varun/node_deletions/" + net + "/" + file + ".topo"
states_fl = "/mnt/e/varun/node_deletions/" + net + "/" + file + "_finFlagFreq.csv"

coherence, freq = state_run(topo_fl, states_fl, 0,"async",1000, "julia")
frust = get_frust(states_fl)
#loops_df = calc_loops(topo_fl)
#pos_wloops, neg_wloops = calc_weighted_loops(loops_df, "P"), calc_weighted_loops(loops_df, "N")
#hamming = compute_hamming(topo_fl, "/mnt/e/EMT_EXP/" + net + "_RANDOM_TOPO/" + net + ".topo", 1)

try:
    fields1 = [file,max(frust), min(frust), np.mean(frust), max(np.multiply(frust, freq)), min(np.multiply(frust, freq)), np.mean(np.multiply(frust, freq))]
except:
    fields1 = [file,hamming, None,None,None,None,None,None]

try:
    fields2 = [max(coherence), min(coherence), np.mean(coherence), max(np.multiply(coherence, freq)), min(np.multiply(coherence, freq)), np.mean(np.multiply(coherence, freq))]
except:
    fields2 = [None,None,None,None,None,None]

try:
    fields3 = [pearsonr(freq, frust)[0], pearsonr(freq, frust)[1], pearsonr(freq, coherence)[0], pearsonr(freq, coherence)[1], pearsonr(frust, coherence)[0], pearsonr(frust, coherence)[1]]
except:
    fields3 = [None, None, None]

with open('/mnt/e/varun/node_deletions/' + net + '_ALL.csv', 'a') as file:
    writer = csv.writer(file)
    writer.writerow(fields1 + fields2 + fields3)

