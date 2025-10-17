#!/usr/bin/env python3

#Boilerplate for the imports needed for the scripts, and a few extras.
import sys

# from statsmodels.stats.proportion import proportion_confint
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import date
import datetime as dt

import csv
import baltic as bt

import math
import os

import glob

#import geopandas as gpd
#import geopandas
import matplotlib.patches as mpatches
#from cartopy import crs as ccrs
import warnings

warnings.filterwarnings('ignore')
#from shapely import Point

#import seaborn as sns
#import scipy
import collections
import numpy as np
#import pandas as pd

import datetime
from datetime import date
# import skbio

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as patches
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.collections import PatchCollection
# from epiweeks import Week, Year


font = {'family' : 'Helvetica',
        'weight' : 'bold',
        'size'   : 18}

from matplotlib.lines import Line2D

mpl.rcParams.update({'font.size': 18})

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}
mpl.rcParams.update(new_rc_params)

plt.rcParams['font.family'] = 'Helvetica'

#Define the tree building function
def tree_builder(treefile,refs,metadata,host_metadata,height=20,title=False):
    my_tree=bt.loadNewick(treefile, absoluteTime=False)
    ref_amplicons=open(refs,"r")
    metadata=open(metadata,"r")
    host_metadata=open(host_metadata,"r")

    save = False
    if not title:    
        title = treefile.split(".")[0] + " and " + refs.split("_")[0] + "_" + refs.split("_")[1]
        save = True
    
    ref_accessions = []
    for line in ref_amplicons:
        if line.split(",")[0] == "reference":
            continue
        else:
            ref_accessions.append(line.split(",")[0])

    ref_metadata = {}
    for line in metadata:
        line = line.lstrip(">")
        ref_metadata[line.split()[0]] = line.split(",")[0]

    #There's lots of ways to parse the host metadata, I've chosen to just pull out the ones that have "Homo sapiens" as host
    host_Hs =[]
    all_found_hosts = []
    no_host_data = []
    for line in host_metadata:
        all_found_hosts.append(line.split("\t")[0])
        if "Homo sapiens" in line.rstrip("\n").split("\t")[2]:
            host_Hs.append(line.split("\t")[0])
        elif "[]" in line.rstrip("\n").split("\t")[2]:
            no_host_data.append(line.split("\t")[0])
        elif "['']" in line.rstrip("\n").split("\t")[2]:
            no_host_data.append(line.split("\t")[0])

    print(host_Hs)
    #print(ref_accessions)
    #print(ref_metadata)

    width = 15
    height = height

    leafs = 0
    for k in my_tree.Objects:
        if k.branchType == 'leaf':
            leafs += 1
            if k.name in ref_accessions:
                if k.name.split(".")[0] in host_Hs:
                    k.traits["label"] = "hit man"
                elif k.name.split(".")[0] not in all_found_hosts:
                    k.traits["label"] = "hit miss"
                elif k.name.split(".")[0] in no_host_data:
                    k.traits["label"] = "hit miss"
                else:
                    k.traits["label"] = "hit"
            elif k.name.split(".")[0] in host_Hs:
                k.traits["label"] = "man"
            elif k.name.split(".")[0] not in all_found_hosts:
                k.traits["label"] = "miss"
            elif k.name.split(".")[0] in no_host_data:
                k.traits["label"] = "miss"
            else:
                k.traits["label"] = "Background"

    fig,ax = plt.subplots(figsize=(width,height),facecolor='w')

    x_attr=lambda k: k.height ## x coordinate of branches will be height attribute
    s_func=lambda k: 140 if 'man' in k.traits['label'] else (70 if 'miss' in k.traits['label'] else 20) ## size of tips
    c_func=lambda k: 'red' if 'man' in k.traits['label'] else ('green' if 'miss' in k.traits['label'] else
        'dimgrey') ## colour
    cb_func=lambda k: "dimgrey"


    # increment = my_tree.treeHeight/150
    my_tree.plotTree(ax,x_attr=x_attr,colour=cb_func) ## plot branches
    my_tree.plotPoints(ax,size=s_func,colour=c_func,x_attr=x_attr) ## plot circles at tips

    #Add text labels to the tips
    target_func=lambda k: k.is_leaf()   ## which branches will be annotated
    text_func=lambda k: ref_metadata[k.name] if 'hit' in k.traits['label'] else " " ## what text is plotted
    text_x_attr= lambda k: k.height+0.05 ## where x coordinate for text is
    kwargs={'va':'center','ha':'left','size': 10, 'color': 'blue'} ## kwargs for text
    my_tree.addText(ax,x_attr=text_x_attr,target=target_func,text=text_func,**kwargs)
    target_func=lambda k: k.is_leaf() ## which branches will be annotated
    text_func=lambda k: ref_metadata[k.name] if 'hit' not in k.traits['label'] else " "
    text_x_attr= lambda k: k.height+0.05 ## where x coordinate for text is
    kwargs={'va':'center','ha':'left','size': 10, 'color': 'grey'} ## kwargs for text
    my_tree.addText(ax,x_attr=text_x_attr,target=target_func,text=text_func,**kwargs)

    mpl.rcParams['font.family'] = 'sans-serif'

    print(len(set(ref_accessions)))        
    print(leafs)

    print("proportion: " + str((len(set(ref_accessions))/leafs)*100))

    [ax.spines[loc].set_visible(False) for loc in ['top','right','left','bottom']]
    ax.tick_params(axis='y',size=0)
    #ax.tick_params(axis='x',size=0)

    ax.set_yticklabels([])
    #ax.set_xticklabels([])

    plt.text(0.000,-1,title,size=30)

    plt.plot([0.000,0.5],[1.75,1.75],linewidth=3, color="dimgrey")
    plt.text(0.000,2.25,"0.5 subs/site",size=30)
    #plt.plot([0.0465,0.0465],[1900,3175],linewidth=3, color="steelblue")
    #plt.plot([0.0465,0.0465],[1200,1800],linewidth=3, color="orange")
    #plt.plot([0.0465,0.0465],[0,1000],linewidth=3, color="red")

    #plt.text(0.047,2437,"A.D.1",size=30, rotation=270, color="steelblue")
    #plt.text(0.047,1400,"A.D.3",size=30, rotation=270, color="orange")
    #plt.text(0.047,400,"A.D.5",size=30, rotation=270, color="red")
    #outfile = open("RSVB.tree", "w")
    #plt.savefig("tree.svg");
    #plt.savefig(f"./figures/{outfile}.png",bbox_inches='tight', 
    #               transparent=True);
    if not save:
        plt.savefig(title,bbox_inches='tight')
    else:
        plt.savefig(refs.split("_")[0] + "_" + refs.split("_")[1],bbox_inches='tight')
    #plt.show()

if __name__ == "__main__":
    # Example usage
    treefile = sys.argv[1]
    refs = sys.argv[2]
    metadata = sys.argv[3]
    host_metadata = sys.argv[4]
    title = sys.argv[5]
    tree_builder(treefile=treefile,refs=refs,metadata=metadata,host_metadata=host_metadata,title=title)