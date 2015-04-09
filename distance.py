#! /usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.PDB import *
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import normal
import pylab as P

#Takes in a protein name and protein data bank file name and returns
#information on the distance distribution of the protein's C-alpha carbons.
#Returns the mean distance, min, max, as well as a histogram.
def analyse(name, f):
	parser=PDBParser()

	#structure=parser.get_structure('1NTI', '1NTI.pdb')
	structure=parser.get_structure(name, f)
	dbl=PDBList()
	dbl.retrieve_pdb_file(name)
	count = 1;
	ds = []
	temp = 0;
	#Iterate through each atom
	for model in structure:
	 for chain in model:
	  for residue in chain:
	   for atom in residue:
	    if atom.get_name() == 'CA': #Checks if the atom is C-alpha
			if count > 1: #If we are past the first atom, calculate the distance between the current and prior atom.
			 dist = atom - temp
			 ds.append(dist)
			 temp = atom
			 count +=1
			else: #Wait for second atom.
			 temp = atom; count += 1
			 		 
	print np.mean(ds)	
	print np.min(ds)
	print ds.index(np.min(ds))	
	print np.max(ds)
	print ds.index(np.max(ds))	
	
	#Creates histogram of the distribution
	plt.hist(ds, bins = 300, range = (2.5,4.5), histtype = 'step')
	plt.title("Distances between alpha carbons")
	plt.xlabel("Value")
	plt.ylabel("Frequency")
	plt.show()	
	
analyse('1NTI', '1NTI.pdb')
analyse('3SIL', '3SIL.pdb')
analyse('1E08', '1E08.pdb')
