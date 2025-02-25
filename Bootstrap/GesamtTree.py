#!python
import os, glob, subprocess, sys
import numpy as np
from natsort import natsorted
from io import StringIO
import Bio
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix

import dendropy
from dendropy import Tree
from dendropy.calculate import treecompare

def MakeNJTreeData():
    # get all PDBs.
    f = open('resDB','a')
    _d_ = glob.glob('*.pdb')
    for i in _d_:
        for j in _d_:
            v = ""
            v = subprocess.Popen('/opt/xtal/ccp4-8.0/bin/gesamt ' + i + ' ' + j + ' | grep "Q-score"', shell = True, stdout = subprocess.PIPE).stdout.read().decode()
            if len(v) > 0:
                f.write('\t'.join([i.split('.')[0],j.split('.')[0],v.split('\n')[0].split()[2]]) + '\n')
            else:
                f.write('\t'.join([i.split('.')[0],j.split('.')[0],'0']) + '\n')
    f.close()
    # Generate distance matrix
    fl = 'resDB'
    d_0 = np.genfromtxt(fl,delimiter='\t',dtype=str,usecols=0,comments = "~")
    d_1 = np.genfromtxt(fl,delimiter='\t',dtype=str,usecols=1,comments = "~")
    d_2 = np.genfromtxt(fl,delimiter='\t',dtype=str,usecols=2,comments = "~")
    merged = []
    for i,j in zip(d_0,d_1):
        merged.append(i + '-' + j)
    merged = np.asarray(merged)
    ud_0 = np.unique(d_0)
    mat = np.zeros(shape=[len(ud_0), len(ud_0)])
    for i in range(len(ud_0)):
        for j in range(len(ud_0)):
            if i > j:
                v1,v2 = '',''
                pat1 = ud_0[i] + '-' + ud_0[j]
                pat2 = ud_0[j] + '-' + ud_0[i]
                idx1 = np.where(pat1 == merged)[0]
                idx2 = np.where(pat2 == merged)[0]
                if len(idx1) == 1:
	                v1 = float(d_2[idx1][0])
                if len(idx2) == 1:       
                        v2 = float(d_2[idx2][0])
                if len(idx1) == 1 and len(idx2) == 0:
                        v2 = v1
                if len(idx2) == 1 and len(idx1) == 0:
                        v1 = v2 
                if len(idx2) == 0 and len(idx1) == 0:
                    v1 = 0
                    v2 = 0
                val = (v1 + v2)/2
                _A_ = round(1 - val,5)
                mat[i,j] = _A_
                mat[j,i] = _A_
    # Generate phyloData
    phy = '#nexus\nBEGIN taxa;\nDIMENSIONS ntax=' + str(len(ud_0)) + ';\nTAXLABELS\n'
    c = 0
    ud_0L = []
    for i in ud_0:
    	c += 1
    	phy += '[' + str(c) + ']' + '_'.join(i.split('_')[0:2]) + '\n'
    	ud_0L.append('_'.join(i.split('_')[0:2]))
    phy += ';\nEND [taxa];\nBEGIN distances;\nDIMENSIONS ntax=' + str(len(ud_0)) + ';\nFORMAT labels diagonal triangle=both;\nMATRIX\n'
    c = 0
    for i in mat:
    	c += 1
    	phy += '[' + str(c) + ']' + '_'.join(ud_0[c-1].split('_')[0:2]) + '\t' + '\t'.join(list(map(str,i))) + '\n'
    dist1_ = np.tril(mat)
    dist2_ = []
    for x in range(0,np.shape(dist1_)[0]):
        dist2_.append(list(dist1_[x][0:x+1]))
    mm_ = _DistanceMatrix(ud_0L,dist2_)
    root_ = DistanceTreeConstructor()
    tree = root_.nj(mm_) 
    Phylo.write(tree,'tree.nex','newick')
    tree = ""
    os.system("sed -i 's,:-[0-9\.]\+,:0.0,g' tree.nex")
    with open('tree.nex') as f:
    	for line in f:
    		tree += line.split('\n')[0]
    with open('Data.nex','a') as f:
    	f.write(phy)
    with open('Taxon.txt','a') as f:
    	for i in ud_0L:
            f.write(i + '\n')
    with open('TreeFormatted.nex','a') as f:
        f.write("#NEXUS\nBEGIN TREES;\n  Tree tree1 = " + tree + '\nEND;\n')
    os.system("sed -i 's,:-[0-9\.]\+,:0.0,g' TreeFormatted.nex")
    os.system('rm resDB')

data = os.getcwd() + '/' + sys.argv[1]
path = os.getcwd()

os.chdir(data)
MakeNJTreeData()
os.chdir(path)

