#! /usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.PDB import *
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import normal
import pylab as P






def analyse2(name, f):
	parser=PDBParser()

	#structure=parser.get_structure('1NTI', '1NTI.pdb')
	structure=parser.get_structure(name, f)
	dbl=PDBList()
	dbl.retrieve_pdb_file(name)
	count = 1;
	ds = []
	temp = 0;
	for model in structure:
	 for chain in model:
	  for residue in chain:
	   for atom in residue:
	    if atom.get_name() == 'CA':
			if count > 1:
			 diff = abs(atom.get_coord() - temp);
			 square = diff*diff;
			 s = sum(square);
			 ds.append(s)
			 
			 
			 temp = atom.get_coord()
			else:
			 temp = atom.get_coord(); count += 1
	print ds	 
	print np.mean(ds)	
	print np.min(ds)
	print ds.index(np.min(ds))	
	print np.max(ds)
	print ds.index(np.max(ds))	
	
	plt.hist(ds, bins = 300, range = (0,20), histtype = 'step')
	plt.title("Distances between alpha carbons")
	plt.xlabel("Value")
	plt.ylabel("Frequency")
	plt.show()
	
def analyse(name, f):
	parser=PDBParser()

	#structure=parser.get_structure('1NTI', '1NTI.pdb')
	structure=parser.get_structure(name, f)
	dbl=PDBList()
	dbl.retrieve_pdb_file(name)
	count = 1;
	ds = []
	temp = 0;
	for model in structure:
	 for chain in model:
	  for residue in chain:
	   for atom in residue:
	    if atom.get_name() == 'CA':
			if count > 1:
			 dist = atom - temp
			 ds.append(dist)
			 temp = atom
			 count +=1
			else:
			 temp = atom; count += 1
			 		 
	print np.mean(ds)	
	print np.min(ds)
	print ds.index(np.min(ds))	
	print np.max(ds)
	print ds.index(np.max(ds))	
	
	plt.hist(ds, bins = 300, range = (2.5,4.5), histtype = 'step')
	plt.title("Distances between alpha carbons")
	plt.xlabel("Value")
	plt.ylabel("Frequency")
	plt.show()	
	
analyse('1NTI', '1NTI.pdb')
analyse('3SIL', '3SIL.pdb')
analyse('1E08', '1E08.pdb')
