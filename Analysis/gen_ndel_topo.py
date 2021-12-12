import numpy as np
import pandas as pd
from itertools import combinations,permutations
import random


def del_topo(table, signal, rule=0, weight=100):
    for node in (signal):
      table = table[table['Source'] != node ]
      table = table[table['Target'] != node ]
    return table

def get_node_labels(ind_list, nodes):
  node_list = []
  for ind in ind_list:
    node_list.append(nodes[ind])
  return node_list

def random_combs(num_nodes,num):
  if num==1 or num==num_nodes:
    poss_list = []
    poss = combinations(np.arange(0,num_nodes,1),num)
    for i in poss:
      poss_list.append(list(i))
    return poss_list
  poss = []
  for i in range(100):
    single = random.sample(list(np.arange(0,num_nodes,1)), num)
    poss.append(tuple(single))
  return poss 

def delete_nodes(topo_file):
  del_nodes_total = []
  topo_df = pd.read_table(topo_file, sep=" ")
  nodes = sorted(list(set(topo_df['Source'].tolist() + topo_df['Target'].tolist())))
  num_nodes = len(nodes)
  num_to_del_list = np.arange(1,2,1)
  for num_to_del in num_to_del_list:
    poss_list = random_combs(num_nodes,num_to_del)
    i = 0
    for poss_ind in poss_list:
      del_nodes = get_node_labels(poss_ind, nodes)
      table = del_topo(topo_df, signal = del_nodes)
      if table.size > 0:
        fl_name = "drosophila_del_" + str(num_to_del) + "_" + str(i) + ".topo"
        table.to_csv("E:/varun/node_deletions/drosophila/" + fl_name, sep = " ", index = None)
        del_nodes_total.append([fl_name, ",".join(del_nodes)])
        i = i+1
  pd.DataFrame(del_nodes_total).to_csv("E:/varun/node_deletions/drosophila_del_annot.csv", index = False, header = False)

#delete_nodes("E:/varun/drosophila.topo")


#Deleting specific nodes 
def nodes_to_delete2(df, net_G):
  df['G'] = df.apply(lambda x: np.mean([abs(x[i]) for i in range(1,5)]), axis=1)
  X = df['Node']
  Y = [abs(i) for i in np.log(df['G']/net_G)]
  sorted_nodes = [x for _, x in sorted(zip(Y, X), reverse=True)]
  return sorted_nodes[:5]


def random_combs2(nodes_to_del, n):
  poss_list = []
  poss = combinations(nodes_to_del,n)
  for i in poss:
    poss_list.append(list(i))
  return poss_list

def delete_nodes2(topo_file, net_G):
  del_nodes_total = []
  net = topo_file.split('/')[-1].replace(".topo", "")
  topo_df = pd.read_table(topo_file, sep=" ")
  nodes = sorted(list(set(topo_df['Source'].tolist() + topo_df['Target'].tolist())))
  num_nodes = len(nodes)
  G_df = pd.read_csv("E:/varun/node_deletions/" + net + "_GroupMetrics.csv", header=0)
  annot_df = pd.read_csv("E:/varun/node_deletions/" + net + "_del_annot.csv", header=0)
  annot_df['Net'] = [i.replace(".topo", "") for i in annot_df['Net']]
  df = pd.merge(G_df, annot_df,on='Net')
  to_del = nodes_to_delete2(df, net_G)
  num_to_del = np.arange(2,len(to_del)+1)
  for n in num_to_del:
    poss_list = random_combs2(to_del, n)
    i = 0
    for poss_ind in poss_list:
      table = del_topo(topo_df, signal=poss_ind)
      if table.size > 0:
        fl_name = net + "_del_" + str(n) + "_" + str(i) + ".topo"
        table.to_csv("E:/varun/node_deletions/" + net + "/" + fl_name, sep = " ", index = None)
        del_nodes_total.append([fl_name, ",".join(poss_ind)])
        i = i+1
  pd.DataFrame(del_nodes_total, columns =['Net', 'Nodes']).to_csv("E:/varun/node_deletions/" + net + "_del_annot2.csv", index = False, header = False)

delete_nodes2("E:/varun/drosophila.topo", net_G = np.mean([0.108, 0.002040816, 0.017142857, 0.038571429]))




