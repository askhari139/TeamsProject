import pandas as pd
import numpy as np
import random 

def get_random_nets(net,n_times):
	for n in range(n_times):
		random_pairs= random.sample(range(0, net.shape[0]), k =2)
		tmp = net.iloc[random_pairs[0],1]
		net.iloc[random_pairs[0],1] = net.iloc[random_pairs[1],1]
		net.iloc[random_pairs[1],1] = tmp
	return net

dest_file = "E:/EMT_EXP/EMT_MET_RANDOM_TOPO/"
is_not_500 = 1
i=0
gen_targets_full = []
while is_not_500:
	gen_net = get_random_nets(pd.read_csv("E:/varun/EMT_MET_reduced.topo", sep = " "), 10)
	gen_targets = gen_net['Target'].tolist()
	if gen_targets not in gen_targets_full:
		i = i + 1
		gen_net.to_csv(dest_file + "EMT_MET_reduced_rand_" + str(i) + ".topo", index = None, sep = " ")
		gen_targets_full.append(gen_targets)
	if len(gen_targets_full) == 500:
		is_not_500 = 0

