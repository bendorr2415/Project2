
import pytest

from clusters import algs


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
	distanceMatrix = algs.makeDistanceMatrix(ligands)
	pc = algs.PartitionClustering()
	pcclusters = pc.cluster(ligands, distanceMatrix, 2)
	ligandIDs = []
	for cluster in pcclusters:
		for ligand in cluster.ligands:
			ligandIDs.append(ligand.ligandID)
	assert ligandIDs == [0,0,0,0,0,1,1,1,1,1] or ligandIDs == [1,1,1,1,1,0,0,0,0,0], 'Partition Clustering Test Failed'
	print('Partition Clustering Test Passed')

def test_hierarchical():
	ligands = read_test_ligands('ligand_information.csv')
	distanceMatrix = algs.makeDistanceMatrix(ligands)
	hc = algs.HierarchicalClustering()
	hcclusters = hc.cluster(ligands, distanceMatrix, 2)
	ligandIDs = []
	for cluster in hcclusters:
		for ligand in cluster.ligands:
			ligandIDs.append(ligand.ligandID)
	assert ligandIDs == [0,0,0,0,0,1,1,1,1,1], 'Hierarchical Clustering Test Failed :('
	print('Hierarchical Clustering Test Passed')
	
def test_make_distance_matrix():
	ligands = read_test_ligands('ligand_information.csv')
	distanceMatrix = algs.makeDistanceMatrix(ligands)
	firstRow = []
	for e in distanceMatrix[0]:
		firstRow.append(e)
	assert firstRow == [1,1,1,1,1,0,0,0,0,0], 'Make Distance Matrix Test Failed'
	print('Make Distance Matrix Test Passed')
	
def test_silhouette_coeff():
	ligands = read_test_ligands('ligand_information.csv')
	distanceMatrix = algs.makeDistanceMatrix(ligands)
	cluster1 = algs.Cluster(algs.Centroid(ligands[0].bitstring))
	cluster2 = algs.Cluster(algs.Centroid(ligands[5].bitstring))
	for i in range(5):
		cluster1.addLigand(ligands[i])
	for i in range(5):
		cluster2.addLigand(ligands[i+5])
	sc = algs.silhouetteCoeff([cluster1, cluster2], distanceMatrix)
	assert sc == 1, 'Silhouette Coefficient Test Failed'
	print('Silhouette Coefficient Test Passed')
	
def test_jaccard_index():
	ligands = read_test_ligands('ligand_information.csv')
	distanceMatrix = algs.makeDistanceMatrix(ligands)
	cluster1 = algs.Cluster(algs.Centroid(ligands[0].bitstring))
	cluster2 = algs.Cluster(algs.Centroid(ligands[5].bitstring))
	cluster3 = algs.Cluster(algs.Centroid(ligands[0].bitstring))
	cluster4 = algs.Cluster(algs.Centroid(ligands[5].bitstring))
	for i in range(5):
		cluster1.addLigand(ligands[i])
		cluster3.addLigand(ligands[i])
	for i in range(5):
		cluster2.addLigand(ligands[i+5])
		cluster4.addLigand(ligands[i+5])
	ji1 = algs.jaccardIndex([cluster1, cluster2], [cluster3, cluster4])
	ji2 = algs.jaccardIndex([cluster1, cluster3], [cluster2, cluster4])
	assert ji1 == 1 and ji2 == 0, 'Jaccard Index Test Failed'
	print('Jaccard Index Test Passed')

def test_csv_io():
	ligands = algs.read_ligand_csv('ligand_information.csv')
	assert ligands[0].ligandID == 0
	assert ligands[0].score == -1.3
	print("'CSV File I/O' Test Passed")
