import numpy as np
import pandas as pd
from itertools import combinations,permutations
import random


def edge_del(df,inds):
  for i in inds:
    df = df.drop(i)
  return df

def random_combs(num_edges,num, total_edges):
  if num==1:
    poss_list = []
    poss = combinations(np.arange(0,total_edges,1),num)
    for i in poss:
      poss_list.append(list(i))
    return poss_list
  poss = []
  for i in range(num_edges):
    single = random.sample(list(np.arange(0,total_edges,1)), num)
    poss.append(tuple(single))
  return poss 

def delete_edges(topo_file):
  del_edges = []
  topo_df = pd.read_table(topo_file, sep=" ")
  num_edges = 100
  total_edges = topo_df.shape[0]
  num_to_del_list = np.arange(1,4,1)
  for num_to_del in num_to_del_list:
    poss_list = random_combs(num_edges,num_to_del, total_edges)
    i = 0
    for poss_ind in poss_list:
      table = edge_del(topo_df, poss_ind)
      if table.size > 0:
        fl = "EMT_RACIPE_edel_" + str(num_to_del) + "_" + str(i) + ".topo"
        table.to_csv("E:/varun/edge_deletions/EMT_RACIPE/" + fl, sep = " ", index = None)
        poss_ind_str = ",".join([str(i) for i in poss_ind])
        del_edges.append([fl,poss_ind_str])
        i = i+1
  pd.DataFrame(del_edges).to_csv("E:/varun/edge_deletions/EMT_RACIPE/edel_annot.csv", index = False, header = False)

delete_edges("E:/varun/EMT_RACIPE.topo")