
import pytest

def read_test_ligands(csv):
	"""
	Returns a list of 10 ligands. The first 5 are ligand 1 from ligand_information.csv, and the next 5 are ligand 2 from the same file.

		Params:
			csv - a csv file containing ligand information in the same format as ligand_information.csv

		Returns:
			a list of Ligand() objects
	"""
	ligands = []
	i = 0 #gotta be a better way to do this
	index = 0
	with open(csv,'r') as f:
		for line in f:
			if i>0:
				for j in range(5):
					ligandID, score, smiles, onbits = line.split(',',3) #set max number of splits = 3
					ligands.append(Ligand(index, ligandID, score, smiles, onbits.replace('"', '').rstrip())) #remove ""'s and eol char
					index += 1
			i+=1
			if i == 3:
				break
	return ligands

def test_partitioning():
	ligands = read_test_ligands('ligand_information.csv')
	distanceMatrix = makeDistanceMatrix(ligands)
	pc = PartitionClustering()
	pcclusters = pc.cluster(ligands, distanceMatrix, 2)
	ligandIDs = []
	for cluster in pcclusters:
		for ligand in cluster.ligands:
			ligandIDs.append(ligand.ligandID)
	assert ligandIDs == [0,0,0,0,0,1,1,1,1,1] or ligandIDs == [1,1,1,1,1,0,0,0,0,0], 'Partition Clustering Test Failed'
	print('Partition Clustering Test Passed')

def test_hierarchical():
	ligands = read_test_ligands('ligand_information.csv')
	distanceMatrix = makeDistanceMatrix(ligands)
	hc = HierarchicalClustering()
	hcclusters = hc.cluster(ligands, distanceMatrix, 2)
	ligandIDs = []
	for cluster in hcclusters:
		for ligand in cluster.ligands:
			ligandIDs.append(ligand.ligandID)
	assert ligandIDs == [0,0,0,0,0,1,1,1,1,1], 'Hierarchical Clustering Test Failed :('
	print('Hierarchical Clustering Test Passed')
