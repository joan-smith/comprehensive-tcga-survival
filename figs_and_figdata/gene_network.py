#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 17:46:34 2021

@author: Joan Smith

https://towardsdatascience.com/visualizing-protein-networks-in-python-58a9b51be9d5
"""

import networkx as nx
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import os


#%%

dropbox_dir = '~/Dropbox/comprehensive-tcga-survival/'
clinical = dropbox_dir + 'raw-data/TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx'
#%%

gene_list = pd.read_csv(os.path.join(dropbox_dir, 'Proliferation analysis', 'e2f targets.csv'), index_col=0)
#%%
proteins = '%0d'.join(gene_list.index.values)
url = 'https://string-db.org/api/tsv/network?identifiers=' + proteins + '&species=9606'
r = requests.get(url)

#%%
lines = r.text.split('\n') # pull the text from the response object and split based on new lines
data = [l.split('\t') for l in lines] # split each line into its components based on tabs
# convert to dataframe using the first row as the column names; drop empty, final row
df = pd.DataFrame(data[1:-1], columns = data[0])

# dataframe with the preferred names of the two proteins and the score of the interaction
interactions = df.loc[:,['preferredName_A', 'preferredName_B', 'score']]
interactions['score'] = interactions['score'].astype(float)
interactions = interactions[interactions['score'] >= interactions['score'].quantile(.75)]

#%%

G=nx.Graph(name='Protein Interaction Graph')

for j, interaction in interactions.iterrows():
    a = interaction['preferredName_A'] # protein a node
    b = interaction['preferredName_B'] # protein b node
    w = float(interaction['score']) # score as weighted edge where high scores = low weight
    G.add_weighted_edges_from([(a,b,w)]) # add weighted eG.dge to graph


#%%

graph_colormap = cm.get_cmap('bwr', 15)

norm = mpl.colors.Normalize(vmin=-15,vmax=15)
# node color varies with Degree
c = [gene_list.loc[g,'Z score'] for g in G]
colors = [graph_colormap(norm(i)) for i in c]


#%%
pos = nx.kamada_kawai_layout(G)
plt.figure(figsize=(19,9)) #,facecolor=[0.7,0.7,0.7,0.4])
nx.draw_networkx(G, pos=pos, with_labels=True, node_color=colors, node_size=1800, edge_color='grey',
                 font_color='black',font_weight='bold',font_size='11')
plt.axis('off')
plt.savefig(os.path.join(dropbox_dir, 'Proliferation analysis', 'rnaseq_e2f_targets.pdf'))

a = np.array([[-15,15]])
plt.figure(figsize=(1.5, 9))
img = plt.imshow(a, cmap='bwr')
plt.gca().set_visible(False)
cax = plt.axes([0.2, 0.6, 0.1, 0.2])
plt.colorbar(orientation="vertical", cax=cax, ticks=[-5, 0, 5, 10, 15])
plt.show()
plt.savefig(os.path.join(dropbox_dir, 'Proliferation analysis', 'rnaseq_e2f_targets.pdf'))


#%%%
#%%%
#%% Methylation network

gene_list = pd.read_csv(os.path.join(dropbox_dir, 'Gene Ontology', 'Methylation events', 'top50 - meth.csv'), index_col=0)
proteins = '%0d'.join(gene_list.index.values)
url = 'https://string-db.org/api/tsv/network?identifiers=' + proteins + '&species=9606'
r = requests.get(url)

#%%
lines = r.text.split('\n') # pull the text from the response object and split based on new lines
data = [l.split('\t') for l in lines] # split each line into its components based on tabs
# convert to dataframe using the first row as the column names; drop empty, final row
df = pd.DataFrame(data[1:-1], columns = data[0])

# dataframe with the preferred names of the two proteins and the score of the interaction
interactions = df.loc[:,['preferredName_A', 'preferredName_B', 'score']]
interactions['score'] = interactions['score'].astype(float)
#interactions = interactions[interactions['score'] >= interactions['score'].quantile(.75)]

#%%

G=nx.Graph(name='Protein Interaction Graph')

for j, interaction in interactions.iterrows():
    a = interaction['preferredName_A'] # protein a node
    b = interaction['preferredName_B'] # protein b node
    G.add_edges_from([(a,b)]) # add weighted eG.dge to graph


#%%
graph_colormap = cm.get_cmap('bwr', 18)

norm = mpl.colors.Normalize(vmin=-15,vmax=12)
# node color varies with Degree
c = [gene_list.loc[g,'Z score'] for g in G]
colors = [graph_colormap(norm(i)) for i in c]

#%%
pos = nx.kamada_kawai_layout(G)
plt.figure(figsize=(19,9)) #,facecolor=[0.7,0.7,0.7,0.4])
nx.draw_networkx(G, pos=pos, with_labels=True, node_color=colors, node_size=1800, edge_color='grey',
                 font_color='black',font_weight='bold',font_size='11')
plt.axis('off')
plt.savefig(os.path.join(dropbox_dir, 'Gene Ontology', 'Methylation events', 'methylation_top50_graph.pdf'))

a = np.array([[-5,15]])
plt.figure(figsize=(1.5, 9))
img = plt.imshow(a, cmap='bwr')
plt.gca().set_visible(False)
cax = plt.axes([0.2, 0.6, 0.1, 0.2])
plt.colorbar(orientation="vertical", cax=cax, ticks=[-5, 0, 5, 10, 15])
plt.show()
plt.savefig(os.path.join(dropbox_dir, 'Gene Ontology', 'Methylation events', 'methylation_top50_colorbar.pdf'))