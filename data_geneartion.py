#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
from pin import ProteinGraph
import pandas as pd
from Bio.PDB import PDBList, PDBIO, PDBParser
import os
from networkx.readwrite import json_graph
import json
import numpy as np
import codecs 

# In[2]:


curated_data=pd.read_csv('curated_data.csv', error_bad_lines=False)
DataPdbs=list(set(curated_data['PDB'].to_list()))
DataPdbs.remove(DataPdbs[0])


# In[3]:


def aaName_dic(aa):
    d={'A':'ALA',
            'R':'ARG',
            'N':'ASN',
            'D':'ASP',
            'C':'CYS',
            'Q':'GLN',
            'E':'GLU',
            'G':'GLY',
            'H':'HIS',
            'I':'ILE',
            'L':'LEU',
            'K':'LYS',
            'M':'MET',
            'F':'PHE',
            'P':'PRO',
            'S':'SER',
            'T':'THR',
            'W':'TRP',
            'Y':'TYR',
            'V':'VAL'}
    return d[aa]


# In[219]:


#pdbl = PDBList()
#pdbl.download_pdb_files(DataPdbs[645:], pdir='msa_pdbs/',file_format='pdb')


# In[4]:


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def convert(o):
    if isinstance(o, np.int64): return int(o) 
        


# In[ ]:


from glob import glob

pdbs=glob('curated_pdbs/*.ent')
missing_wholePssms=[]
for pdb in pdbs:
    pdb_code=pdb.split('/')[-1][3:-4]
    print(pdb_code)
    try:
        graph=ProteinGraph(pdb)
    except:
        continue
        
    d=curated_data[curated_data['PDB'] == pdb_code].dropna()
    catalytic=[str(i[1])+str(i[2])+str(i[0]) for i in zip([x.upper() for x in d['PDB code'].to_list()], 
                                                          d['chain/kegg compound'].to_list(),
                                                          d['resid/chebi id'].to_list())]
    
    try:
        pssm=open(glob('proccessed_pssms/{pdb_code}*.pssm'.format(pdb_code=pdb_code))[0]).read().splitlines()[3:]
        pssm_nodes={}
    except:
        missing_wholePssms.append(pdb_code)
        continue 


    for p in pssm:
        allP=p.split()
        if len(allP[2])==1:
            node=allP[0]+allP[1]+aaName_dic(allP[2])
            pssmL=allP[3:23]
            pssm_nodes[node]= [int(score) for score in pssmL]

    for n in graph.nodes: #add pssm line to features
        graph.nodes[n]['features'][0].extend(pssm_nodes[n])
        break

    #for e in graph.edges:
    #    graph.edges[e]['features']=list(graph.edges[e]['features'])
        
    labeled_dic={key:0.0 for key in list(graph.nodes)}
    for cat in catalytic:
        if cat in graph.nodes:
            labeled_dic[cat]=1.
            
    graph = nx.node_link_graph(graph)
    g = dgl.DGLGraph()
      
      
    with open('json/'+pdb_code+'.json', 'w') as outfile1:
        G=json.dumps(json_graph.node_link_data(graph), cls=NumpyEncoder,default=convert)
        dic=json.dumps(labeled_dic)
        json.dump([G, dic],outfile1)
    #os.remove('pdbs/pdb{pdb}.ent'.format(pdb=pdb))
    break





