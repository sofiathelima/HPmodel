# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 22:45:57 2022

@author: sofia, kevin, abe
"""
from aminoAcids import *
import copy
import math
import matplotlib.pyplot as plt
import numpy as np
from presentation import *
import pandas as pd
import proteomics
import seaborn as sns
import random
import time
import prody as prdy

dummyProtein5 = ["H","H","P","H","H"]    #Protein sequence being folded
dummyProtein10 = ["H","H","P","H","H","H","H","P","H","H"]
dummyProtein20 = ["H","H","P","H","H","H","H","P","H","H","H","H","P","H","H","H","H","P","H","H"]

M = 20
numstarts = 100      #number of evaluations to perform
directions = [0,1,2,3,4,5]			#0,1 x-axis, 2,3 y-axis, 4,5 z-axis one for each direction of pivot
f_QUARK = 'model1_qa15349.pdb'
f_1si4 = '1si4.pdb'

save_num = '6'

class AminoAcid:
	def __init__(self, dir, classification, r, c, z):
		self.dir = dir #direction is either 0,1,2,3,4,5. Corresponding directions listed above.
		self.type = classification
		self.r = r
		self.c = c
		self.z = z

# isValid checks all indexes prior to the current amino acid and compares the coordinates to
# itself. Returns False if non self avoiding; returns true otherwise.
def isValid(idx, aa, aaSeq):

	for n in aaSeq[:idx-1]:
		if n.r == aa.r and n.c == aa.c and n.z == aa.z:
			return False
	return True

# step Returns amino acid with coordinates incremented from its predescesor corresponding to the 
# given direction: 0 = x-axis-"up"; 1 = x-axis-"down"; 2 = y-axis-"up"; 3 = y-axis-"down"; 
# 4 = z-axis-"up"; 5 = z-axis-"down".
def step(dir, aaSeq, idx, protein):

	r = aaSeq[idx-1].r
	c = aaSeq[idx-1].c
	z = aaSeq[idx-1].z
	if dir == 0:
		r+=1.0
	if dir == 1:
		r-=1.0
	if dir == 2:
		c+=1.0
	if dir == 3:
		c-=1.0
	if dir == 4:
		z+=1.0
	if dir == 5:
		z-=1.0
	
	aa = AminoAcid(dir, protein[idx], r, c, z)

	return aa

# fold Maintains a self avoiding fold by iterating through the protein sequence and making a step 
# in a random direction at each amino acid. Returns a list of amino acids with updated coordinates.
def foldRandom(protein):
	# initialize first step
	dir = random.choice(directions)
	aa = AminoAcid(dir, protein[0], float(len(protein)), float(len(protein)), float(len(protein)))
	aaSeq = [aa]

	# while the iterator has not reached the length of the protein: 
	# check if the step results in a valid protein, otherwise remove the last step
	idx = 1
	while idx < len(protein):
		if isValid(idx, aa, aaSeq):
			dir = random.choice(directions)
			aa = step(dir,aaSeq, idx, protein)
			aaSeq.append(aa)
			idx+=1
		else:
			aaSeq = aaSeq[:-1]
			idx-=1
	

	# previous implementation: takes too long because it starts over when its non-self-avoiding
	# idx = 1
	# while isValid(idx, aa, aaSeq) and idx < len(protein):
	# 	print(f'i:{idx}')
	# 	aa = step(aaSeq, idx)
	# 	aaSeq.append(aa)
	# 	idx+=1


	return aaSeq

# inputs a sequence of amino acids and an H-H energy penalty, outputs total energy of chain
def energy(seq, penalty):
    energy = 0
    for i in range(len(seq)):
        # range through all other amino acids, determine if amino acid is close enough to affect energy
        for j in range(len(seq)):
            if i != j and i != j+1 and i != j-1:
                # see if both 'H'
                if seq[i].type == 'H' and seq[j].type == 'H':
                # see if adjacent
                # place points in numpy array
                    aa_1 = np.array([seq[i].r, seq[i].c, seq[i].z])
                    aa_2 = np.array([seq[j].r, seq[j].c, seq[j].z])
                    # calculate distance
                    # dist = np.sum((aa_1-aa_2)**2, axis=0)
                    # dist = np.sqrt(dist)
                    dist = np.linalg.norm(aa_1-aa_2)

                    # see if adjacent
                    if dist == 1: # adjacent
                        energy += penalty

    return energy

# randomFoldSampler Executes numstarts number of trials, where each trial is a randomly folded valid protein
# structure. It compares the energy of each structure, and returns the amino acid sequence with the 
# highest energy. This will provide a baseline to demonstrate the efficiency of Simulated Annealing.
def randomFoldSampler(protein):
	bestSeq = foldRandom(protein)
	best_energy = energy(bestSeq, -1)
	all_energies = [best_energy]
	
	for i in range(1, numstarts):
		aaSeq = foldRandom(protein)
		curr_energy = energy(aaSeq, -1)
		all_energies += [curr_energy]
		if curr_energy < best_energy:
			best_energy = curr_energy
			bestSeq = aaSeq
	
	return bestSeq, all_energies

# getFoldSeq Encodes a sequence of amino acids to a sequence of fold directions ([]integers).
def getFoldSeq(aaSeq):
	foldSeq = []
	for aa in aaSeq:
		foldSeq += [aa.dir]
	return foldSeq

# foldSeq_to_aaSeq Translates the sequence of fold directions {0,1,2,3,4,5} to a sequence of amino
# acids with appropriate coordinates, corresponding to a pivot at the ith position of the
# initialized fold structure. Truncates if non-self-avoiding.
def foldSeq_to_aaSeq(newFold, protein):
	aa = AminoAcid(newFold[0], protein[0], float(len(protein)), float(len(protein)), float(len(protein)))
	translated_aaSeq = [aa]

	idx = 1
	while idx < len(protein):
		if isValid(idx, aa, translated_aaSeq):
			dir = newFold[idx]
			aa = step(dir, translated_aaSeq, idx, protein)
			translated_aaSeq.append(aa)
			idx+=1
		else:
			return translated_aaSeq

	# for j in range(i,len(aaSeq)):
	# 	if isValid(j, aaSeq[j], aaSeq):
	# 		aa = step(newFold[j], aaSeq, j, protein)
	# 		translated_aaSeq += [aa]
	# 	else:
	# 		print(f'not valid here: {i}')
	# 		return translated_aaSeq
	
	return translated_aaSeq

# getNeighbors Searches for all valid "neighbor" of the given peptide structure, using a DP approach
#  to find and calculate coordinates of "neighboring structures"
# that differ from the initialized structure one position at a time.
def getNeighbors(aaSeq, protein):
	print("Entering getNeighbors")
	# Copy sequence of fold directions for each aa in initialized structure
	foldSeq = getFoldSeq(aaSeq)
	
	# for every position in the sequence, pivot the structure in each direction, and refold from i
	neighbors = []
	for i in range(0,len(aaSeq)):
		for d in directions:
			if d != aaSeq[i].dir:
				newFold = foldSeq
				newFold[i] = d
				newSeq = foldSeq_to_aaSeq(newFold, protein) # truncates if non self avoiding
				# print(f'len(newSeq)={len(newSeq)}')
				if len(newSeq) == len(aaSeq):
					neighbors += [newSeq]
	print(f'Added {len(neighbors)} neighbors to list.')

	return neighbors

# bestNeighbor Finds and calculate coordinates of "neighboring structures"
# that differ from the initialized structure one position at a time.
def bestNeighbor(aaSeq, protein):
	# get all possible neighboring amino acid structures (some shorter than initial aaSeq)
	print("\nEntering bestNeighbor")
	neighbors = getNeighbors(aaSeq, protein)
	bestNeighbor = aaSeq
	best_energy = energy(aaSeq, -1)
	all_energies = [best_energy]

	print("Re-Entering bestNeighbor. About to enter energy comparison...")
	# check each neighboring fold structure for validity, and if so compare energies
	for nSeq in neighbors:
		curr_energy = energy(nSeq, -1)
		if curr_energy < best_energy:
			bestNeighbor = nSeq
			best_energy = curr_energy
		all_energies += [curr_energy]
	
	# print(f'len(bestNeighbor)={len(bestNeighbor)}; should be {len(protein)}')
	print(f'# valid energies={len(all_energies)}')
	
	return bestNeighbor, best_energy, all_energies



# SA_Sampler Uses a simulated annealing technique to explores "neighboring structures" of peptides.
# It initializes a randomly folded valid protein, and pivots this initial structure at a random
# position in the sequence, accepting the pivoted structure with decreasing probability over time,
# until we have not changed our current structure for M iterations. It compares energies of viable 
# accepted structures and returns the peptide structure with the highest energy. It generates one 
# graph: Time-vs-iteration.
def SA_Sampler(protein):

	starttemp = 10 * len(protein) # 200
	temps = []
	all_final_energies = []
	all_all_energies = []
	longest_run = 0
	all_folds = []

	for i in range(numstarts):
		temp = starttemp
		rate = float(starttemp)/float(1000)
		boltzmann = 0.001987

		aaSeq = foldRandom(protein)
		bestSeq, best_energy, all_energies = bestNeighbor(aaSeq, protein)
		all_all_energies += [all_energies]
		if len(all_energies) > longest_run:
			longest_run = len(all_energies)
		nSeq = bestSeq

		all_accepted_energy = [best_energy]
		num_accepted = 1

		m = 0

		while m < M:
			temps += [temp]
			print(f'i: {i}, m: {m}, starttemp:{starttemp}, temp={temp}')

			neighbors = getNeighbors(bestSeq, protein)
			aaSeq = random.choice(neighbors)

			nSeq, curr_energy, all_energies = bestNeighbor(aaSeq, protein)

			if nSeq != aaSeq:

				if curr_energy < best_energy:
					all_accepted_energy += [curr_energy]
					num_accepted += 1
					best_energy = curr_energy
					bestSeq = nSeq
					aaSeq = nSeq

				else:
					d = curr_energy-best_energy
					accept = pow(math.e, -1*(d)/(boltzmann*temp)) # Calculate acceptance value using Metropolis algorithm
					
					rand_num = random.random() #random number between 0 and 1

					if rand_num < accept:
						all_accepted_energy += [curr_energy]
						num_accepted += 1
						aaSeq = nSeq

						if curr_energy <= best_energy:
							best_energy = curr_energy
							bestSeq = nSeq

				print(f'temp before: {temp}')
				temp -= rate
				print(f'temp -= rate: {temp}\n')
			else:
				m += 1
				print(f'\nm={m}')
		
		all_final_energies += [best_energy]
		all_folds.append(nSeq)

	
	x = np.linspace(1, len(temps), len(temps))
	proteomics.output_plot(x, temps,"iteration","temperature", "Line Plot of Temperature vs Algorithm Iteration (SA)", f'Temp-vs-Iteration({len(protein)}aa{numstarts}runs)({save_num})')

	all_runs_x = np.linspace(1, longest_run, longest_run)
	proteomics.output_overlay_graph(all_runs_x, all_all_energies,"iteration", "energy score", "All energies of accepted", f'accepted_all_energies({len(protein)}aa{numstarts}runs({save_num})')
	


	return bestSeq, all_final_energies, all_folds

if __name__ == "__main__":

	start = time.time()

	with open('seq.txt') as f:
		#read txt file line by line
		fin = f.read()
		protein_name, seq = fin.strip().split('\n')
	# print(seq)
	protein = HP_chain(amino_acid_dict,seq,True)
	# print(protein)


	aaSeq_mc, energies_mc = randomFoldSampler(protein)
	proteomics.output_distribution(energies_mc, f'randomFoldSampler-energy-distribution({len(aaSeq_mc)}aa{numstarts}trials)({save_num})')
	# proteomics.Plot_ProteinLength_vs_Runtime_MC(protein)
	hp_mc = proteomics.output_pdb(aaSeq_mc,f'HP-MC({len(aaSeq_mc)}aa{numstarts}runs)({save_num})')
	hp_mc_scaled = proteomics.scale(hp_mc, f'HP-MC({len(aaSeq_mc)}aa{numstarts}runs)({save_num})(scaled)')

	aaSeq_sa, energies_sa, allFolds_sa = SA_Sampler(protein)
	proteomics.output_distribution(energies_sa, f'SA_Sampler-energy-distribution({len(aaSeq_sa)}aa{numstarts}trials)({save_num})')
	# proteomics.Plot_ProteinLength_vs_Runtime_SA(protein)
	hp_sa = proteomics.output_pdb(aaSeq_sa, f'HP-SA({len(aaSeq_sa)}aa{numstarts}runs)({save_num})')
	hp_sa_scaled = proteomics.scale(hp_sa, f'HP-SA({len(aaSeq_sa)}aa{numstarts}runs)({save_num})(scaled)')

	proteomics.output_distributions(energies_mc, energies_sa)

	proteomics.CompareModels(hp_mc, hp_sa, hp_mc_scaled, hp_sa_scaled)
	proteomics.CompareScoring(energies_sa, allFolds_sa)


	print(f'\nTotal program time = {(time.time()-start)/60} min')