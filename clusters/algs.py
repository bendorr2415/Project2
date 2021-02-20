import numpy as np
import random


class HierarchicalClustering():
	def __init__(self):
		pass

	def cluster(self, ligands, distanceMatrix, numClusters, verbose = False):
		"""
		Clusters the given ligands into numClusters clusters using the average linkage hierarchical clustering algorithm.  Treats each individual ligand
		as its own cluster from the start.  Finds the two most similar clusters using the given distance matrix and combines them into the same
		cluster. Calculates the clusters' average distance to all other clusters and adds this array of distances to the last row and last column of
		the distance matrix, then removes the rows and columns that correspond to the individual clusters that were combined.  Repeats this process until
		there are numClusters clusters remaining.

			Params:
				ligands - a list of Ligand() objects

				distanceMatrix - a len(ligands) x len(ligands) 2D array filled with Tanimoto Coefficients corresponding
				to each row and column ligand pair

				numClusters - an int dictating how many Cluster() objects will be returned

				verbose - when true, prints the iteration number of the algorithm as it's running

			Returns:
				a list of Cluster() objects
		"""

		# Make sure each dimension of the distance matrix is the same length as len(ligands), in case I use a subset of the ligands dataset
		# Made this extra complicated in case I ever want to use slices from the middle of the ligands dataset
		distanceMatrix = distanceMatrix[ligands[0].ligandIndex:ligands[-1].ligandIndex+1, ligands[0].ligandIndex:ligands[-1].ligandIndex+1]

		matrixHeader = list(map(lambda ligand:[ligand], ligands)) #returns a list of lists

		while len(matrixHeader) > numClusters:
			if verbose:
				print('Clustering step number %d' % (len(ligands) - len(matrixHeader)))

			# Determine the current max tanimoto coefficient in the distance matrix (which determines the next clustering event)
			maxVal = 0
			maxValIndeces = [0,0]
			for i in range(len(distanceMatrix)):
				for j in range(len(distanceMatrix[i])):
					if i != j:
						if maxVal < distanceMatrix[i][j]:
								maxVal = distanceMatrix[i][j]
								maxValIndeces = [i,j]

			# Populate the new row (and column) of the distance matrix with the minimum values from each of the combined row's indeces
			# Do not include -1's in the combined row (these correspond to a ligand's distance to itself)
			newRow = []
			for k in range(len(distanceMatrix[maxValIndeces[1]])):

				weightedVal1 = distanceMatrix[maxValIndeces[0]][k] * len(matrixHeader[maxValIndeces[0]])
				weightedVal2 = distanceMatrix[maxValIndeces[1]][k] * len(matrixHeader[maxValIndeces[1]])
				aveVal = (weightedVal1 + weightedVal2) / (len(matrixHeader[maxValIndeces[0]])+len(matrixHeader[maxValIndeces[1]]))
				newRow.append([aveVal])

			newRow = np.array(newRow)

			# Add the new row to the distance matrix, in both the rows and columns
			distanceMatrix = np.hstack((distanceMatrix, newRow))
			newRow = np.append(newRow, [1]) # 1 here because that's the tanimoto between the same ligand
			distanceMatrix = np.vstack((distanceMatrix, newRow))

			# Remove the rows and columns that were just clustered together from the distance matrix
			distanceMatrix = np.delete(distanceMatrix, [maxValIndeces[0], maxValIndeces[1]], axis = 0)
			distanceMatrix = np.delete(distanceMatrix, [maxValIndeces[0], maxValIndeces[1]], axis = 1)

			# Add a new column to the matrix header, which represents the combination of the two combined clusters
			matrixHeader.append(matrixHeader[maxValIndeces[0]] + matrixHeader[maxValIndeces[1]])

			# Remove the two individual clusters that were just combined from the matrix header
			matrixHeader.pop(max(maxValIndeces))
			matrixHeader.pop(min(maxValIndeces))

		clusters = []
		for cluster in matrixHeader:
			fake_centroid = Centroid(ligands[0].bitstring) #relic of partition clustering implementation
			new_cluster = Cluster(fake_centroid)
			for ligand in cluster:
				new_cluster.addLigand(ligand)
			clusters.append(new_cluster)

		return clusters



class PartitionClustering():
	def __init__(self):
		pass


	def cluster(self, ligands, distanceMatrix, numClusters, verbose=False):
		"""
		Clusters the given ligands into numClusters clusters using the K-Means++ partition clustering algorithm.  Initializes numClusters
		Centroid() objects, and uses these to initialize numClusters Cluster() objects.  For each ligand, calculates which cluster's centroid
		is closest to the ligand and adds that ligand to that cluster.  Once all ligands have been added to a cluster, ensures that no clusters
		have no ligands.  If a cluster is empty (has no ligands), then removes one ligand from a non-empty cluster and adds this ligand to the 
		empty cluster.  Then calculates a new centroid based on each clusters' new list of ligands, and determines if the centroids have moved from
		the last iteration.  If they have not moved, the clustering is complete.  If they have moved, then the clusters' centroids are updated, and
		the process of assigning ligands begins again with the new centroids.  This process repeats until either the centroids stop moving, or the
		algorithm has run for 100 iterations.

			Params:
				ligands - a list of Ligand() objects

				distanceMatrix - a len(ligands) x len(ligands) 2D array filled with Tanimoto Coefficients corresponding
				to each row and column ligand pair

				numClusters - an int dictating how many Cluster() objects will be returned

				verbose - when true, prints the iteration number of the algorithm as it's running, as well as the number of bits that changed
				across all the centroids' .bitstring attributes in that iteration (which is a measure of how much the centroids moved in that iteration)

			Returns:
				a list of Cluster() objects

		"""

		centroids = self.initializeCentroids(numClusters, ligands, distanceMatrix)

		clusters = self.initializeClusters(centroids) # clusters is a list of cluster objects with centroid attributes

		# centroidUpdate is the sum of each cluster's calcCentroidMoved, which is 0 when a centroid doesn't move
		prevUpdates = []
		centroidUpdate = 1
		iterations = 0

		while centroidUpdate != 0 and iterations < 100:
			centroidUpdate = 0
			iterations += 1

			# Empty the ligand list in each cluster. They'll be repopulated with the new centroids
			for cluster in clusters:
				cluster.resetLigands()

			# Determine which cluster each ligand should join
			for ligand in ligands:
				tanCo = ligand.tanimotoCoefficient(clusters[0].centroid) # this works because centroid has .onbits attribute
				#print('tanCo of ligand with clusters[0].centroid = %.2f' % tanCo)
				i, joinClusterIndex = 0, 0
				for cluster in clusters:
					if tanCo > ligand.tanimotoCoefficient(cluster.centroid):
						tanCo = ligand.tanimotoCoefficient(cluster.centroid)
						joinClusterIndex = i
					i += 1
				clusters[joinClusterIndex].addLigand(ligand)

			# If any clusters are empty, add one ligand to the cluster from a cluster with 2 or more ligands
			for cluster in clusters:
				if len(cluster.ligands) == 0:
					for cluster2 in clusters:
						if len(cluster2.ligands) > 1:
							cluster.addLigand(cluster2.ligands[-1])
							cluster2.removeLigand(cluster2.ligands[-1])

			# Determine whether the centroids moved
			for cluster in clusters:
				centroidUpdate += cluster.calcCentroidMoved()
				cluster.updateCentroid()
			prevUpdates.append(centroidUpdate)

			if verbose:
				print('Iteration number %d. Number of changed centroid bits = %d' % (iterations, centroidUpdate))

		return clusters
		

	def initializeClusters(self, centroids):
		"""
		For each Centroid() object in centroids, creates a new Cluster() object with the Centroid.  Returns a list of these Cluster objects

			Params:
				centroids - a list of Centroid() objects

			Returns:
				a list of Cluster() objects
		"""
		clusters = []
		for centroid in centroids:
			newCluster = Cluster(centroid)
			clusters.append(newCluster)
		return clusters


	def initializeCentroids(self, numClusters, ligands, distanceMatrix):
		"""
		Picks numClusters centroids from the ligands list. The probability of each ligand being chosen as a
		centroid is inversely proportional to its distance from the closest centroid that has already been
		chosen.

			Params:
				numClusters - an int dictating how many centroids will be returned

				ligands - a list of Ligand() objects

				distanceMatrix - a len(ligands) x len(ligands) 2D array filled with Tanimoto Coefficients corresponding
				to each row and column ligand pair

			Returns:
				a list of Centroid() objects
		"""
		indeces = []

		# Randomly choose one ligand as a first centroid
		i = random.randint(0,len(ligands)-1) # randint is inclusive of second number
		indeces.append(i)

		# The probability of choosing a second ligand is inversely proportional to the
		# second ligand's similarity to the first
		weights = distanceMatrix[i]
		weights = 1 - weights

		# Choose the rest of the centroids with their associated probabilities, updating these
		# probabilities after each new centroid is chosen
		for c in range(numClusters-1): # -1 because we've already picked one index (/centroid)
			i = random.choices(np.arange(len(distanceMatrix)), weights)
			newWeights = 1 - distanceMatrix[i]

			# update weights with the lowest similarity values for each of the already-chosen centroids
			weights = np.minimum(weights, newWeights)[0] #this function returns a one-element list
			indeces.append(i[0]) # i is a one-element list

		centroids = []
		for i in indeces:
			centroids.append(Centroid(ligands[i].bitstring))

		return centroids


class Centroid():
	def __init__(self, bitstring):
		self.bitstring = bitstring
		self.onbits = self.bitstringToOnbits(bitstring)


	def bitstringToOnbits(self, bitstring):
		"""
		Calculates and returns the densified representation of the bitstring by returning a list of indeces at
			which the bitstring had 1's.

			Params:
				bitstring - a 1024 element list of 1's and 0's

			Returns:
				a list of indeces at which the bitstring contained 1's
		"""
		onbits = []
		for i in range(len(bitstring)):
			if bitstring[i] == 1:
				onbits.append(i)
		return onbits



class Cluster():
	def __init__(self, centroid):
		self.centroid = centroid
		self.ligands = []

	def set_centroid(self, centroid):
		self.centroid = centroid

	def addLigand(self, ligand):
		self.ligands.append(ligand)

	def removeLigand(self, ligand):
		self.ligands.remove(ligand)

	def resetLigands(self):
		self.ligands = []

	def calcCentroid(self, ligands):
		"""
		Calculates and returns a bitstring that represents the average bitstring across all the ligands in the ligand list parameter. Any bit
		position that is 1 in most of the given ligand bitstrings will have a 1 in the final bitstring, and the same is true for 0s.
		Creates a new Centroid object with this bitstring and returns it.
		If the ligand list is empty, it returns the Cluster's current centroid object.

			Params:
				ligands - a list of ligand objects

			Returns:
				a centroid object.  Either a new centroid object with the new average bitstring, or the current centroid attribute of the Cluster object
		"""
		if len(ligands) > 0:
			centroidBitstring = np.zeros(1024, dtype=int)
			for ligand in ligands:
				centroidBitstring += ligand.bitstring
			centroidBitstring = centroidBitstring / len(ligands)
			centroidBitstring = np.rint(centroidBitstring)
			return Centroid(centroidBitstring)
		else:
			return self.centroid


	def calcCentroidMoved(self):
		"""
		Returns the number of bits that changed between the bitstrings of the Cluster() object's previous centroid and its new centroid.

			Params:
				None

			Returns:
				an int, representing the total number of bits that changed between the Cluster's old and new centroids' bitstrings
		"""
		j = 0
		newBitstring = self.calcCentroid(self.ligands).bitstring
		for i in range(len(self.centroid.bitstring)):
			if self.centroid.bitstring[i] != newBitstring[i]:
				j+=1
		return j


	def updateCentroid(self):
		"""
		Sets the Cluster() object's centroid attribute to what is returned from the Cluster() object's calcCentroid method.

			Params:
				None

			Returns:
				None
		"""
		self.centroid = self.calcCentroid(self.ligands)



class Ligand():
	def __init__(self, ligandIndex, ligandID, score, smiles, onbits):
		self.ligandIndex = ligandIndex
		self.ligandID = int(ligandID)
		self.score = float(score)
		self.smiles = smiles
		self.onbits = list(map(int, onbits.split(',')))
		self.bitstring = self.onbitsToBitstring(onbits)

	def __str__(self):
		return "%d" % self.ligandID

	def onbitsToBitstring(self, onbits):
		"""
		Converts a densified string of onbits to an expanded, 1024-bit-long bitstring.

			Params:
				onbits - a list of onbits for a 1024-bit-long bitstring.

			Returns:
				A 1024-element numpy array of 1s and 0s.
		"""
		onbits = onbits.split(',')
		bitstring = np.zeros(1024, dtype=int)
		for ob in onbits:
			bitstring[int(ob)] = 1
		return bitstring


	def tanimotoCoefficient(self, other):
		"""
		Calculates the Tanimoto Coefficient between two ligands' list of onbits.

			Params:
				An object that contains an .onbits attribute, such as a Ligand() or Centroid() object

			Returns:
				a float
		"""
		selfOnbits = set(self.onbits)
		otherOnbits = set(other.onbits)

		numerator = selfOnbits.intersection(otherOnbits)
		denomenator = selfOnbits.union(otherOnbits)

		return len(numerator) / len(denomenator)


def read_ligand_csv(csv):
	"""
	Creates a new Ligand() object for each line in the ligand_information.csv file, setting the Ligand() object's attributes
		to the appropriate values from the file.

		Params:
			csv - a csv file containing ligand information in the same format as ligand_information.csv

		Returns:
			a list of Ligand() objects
	"""
	ligands = []
	i = 0 #gotta be a better way to do this
	with open(csv,'r') as f:
		for line in f:
			if i>0:
				ligandID, score, smiles, onbits = line.split(',',3) #set max number of splits = 3
				ligands.append(Ligand(i-1, ligandID, score, smiles, onbits.replace('"', '').rstrip())) #remove ""'s and eol char
			i+=1
	return ligands


def silhouetteCoeff(clusters, distanceMatrix):
	"""
	Calcualtes and returns the silhouette coefficient for a set of clusters.

		Params:
			clusters - a list of Cluster() objects

			distanceMatrix - a len(ligands) x len(ligands) 2D array filled with Tanimoto Coefficients corresponding
				to each row and column ligand pair

		Returns:
			a float
	"""

	# the dM has Tanimoto Coeffs. Now change its values to be proportional to distances
	distanceMatrix = 1 - distanceMatrix

	s = []
	a, b = [], []
	for c in range(len(clusters)):
		for ligand in clusters[c].ligands:

			# For each ligand in the current cluster, calculate its 'a' term, its cohesion, to all the other
			#  ligands in its cluster
			da = []
			for ligand2 in clusters[c].ligands:
				if ligand.ligandIndex != ligand2.ligandIndex:
					da.append(distanceMatrix[ligand.ligandIndex][ligand2.ligandIndex])
				elif ligand.ligandIndex == ligand2.ligandIndex:
					da.append(0)
			if len(clusters[c].ligands) == 1:
				a.append(0)
			else:
				a.append((1 / (len(clusters[c].ligands)-1)) * sum(da))

			# Calculate the 'b'' term, the separation, between the current ligand and all the ligands in the other clusters
			bTemp = []
			db = []
			for c2 in range(len(clusters)):
				if c2 != c:
					for ligand2 in clusters[c2].ligands:
						if ligand.ligandIndex == ligand2.ligandIndex:
							db.append(0)
						else:
							db.append(distanceMatrix[ligand.ligandIndex][ligand2.ligandIndex])
					bTemp.append((1 / len(clusters[c2].ligands)) * sum(db))
			b.append(min(bTemp))

	# calc s for each ligand
	for i in range(len(a)):
		if max(a[i], b[i]) == 0:
			s.append(0)
		else:
			s.append((b[i] - a[i]) / max(a[i], b[i]))
	sc = sum(s) / len(s)
	return sc


def makeDistanceMatrix(ligands, verbose = False):
	"""
	Fills a len(ligands) x len(ligands) matrix with the Tanimoto Coefficients of each row and column pair of ligands.
	A TC of 1 means the ligands have identical fingerprints.

		Params:
			ligands - a list of Ligand() objects

		Returns:
			a 2D numpy array
	"""
	matrix = []
	for i in range(len(ligands)):
		if verbose:
			print('filling row %d' % i)
		row = []
		for j in range(len(ligands)):
			if i == j:
				row.append(1) # the Tanimoto Coefficient between a ligand and itself is 1
			else:
				row.append(ligands[i].tanimotoCoefficient(ligands[j]))
		matrix.append(row)
	return np.array(matrix)


def jaccardIndex(clusters1, clusters2):
	"""
	Calculates and returns the Jaccard Index of two sets of clusters.
	clusters1 must have a subset of or the same set of ligands as clusters2.

		Params:
			hcclusters - a list of Cluster() objects. A subset or the same set of ligands as in pcclusters.

			pcclusters - a list of Cluster() objects

		Returns:
			A float
	"""
	allLigandIDs = []
	clusters1Pairs = []
	clusters2Pairs = []

	# Record all the ligands and ligand pairs in clusters1, which may have a subset of the ligands in clusters2
	for cluster in clusters1:
		for l1 in range(len(cluster.ligands)):
			allLigandIDs.append(cluster.ligands[l1].ligandID)
			for l2 in range(len(cluster.ligands)-l1-1):
				clusters1Pairs.append((cluster.ligands[l1].ligandID, cluster.ligands[l1+l2+1].ligandID))

	# Record all the ligand pairs in clusters2 in which both ligands are also in clusters1
	for cluster in clusters2:
		for l1 in range(len(cluster.ligands)):
			for l2 in range(len(cluster.ligands)-l1-1):
				if cluster.ligands[l1].ligandID in allLigandIDs and cluster.ligands[l2].ligandID in allLigandIDs:
					clusters2Pairs.append((cluster.ligands[l1].ligandID, cluster.ligands[l1+l2+1].ligandID))

	clusters1Set = set(clusters1Pairs)
	clusters2Set = set(clusters2Pairs)

	numerator = clusters1Set.intersection(clusters2Set)
	denomenator = clusters1Set.union(clusters2Set)			

	return len(numerator) / len(denomenator)






#######
## UNIT TESTS
#######


def test_csv_io():
	ligands = read_ligand_csv('ligand_information.csv')
	assert ligands[0].ligandID == 0
	assert ligands[0].score == -1.3
	print("'CSV File I/O' Test Passed")


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


def test_silhouette_coeff():
	ligands = read_test_ligands('ligand_information.csv')
	distanceMatrix = makeDistanceMatrix(ligands)
	cluster1 = Cluster(Centroid(ligands[0].bitstring))
	cluster2 = Cluster(Centroid(ligands[5].bitstring))
	for i in range(5):
		cluster1.addLigand(ligands[i])
	for i in range(5):
		cluster2.addLigand(ligands[i+5])
	sc = silhouetteCoeff([cluster1, cluster2], distanceMatrix)
	assert sc == 1, 'Silhouette Coefficient Test Failed'
	print('Silhouette Coefficient Test Passed')


def test_partition_clustering():
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


def test_hierarchical_clustering():
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


def test_make_distance_matrix():
	ligands = read_test_ligands('ligand_information.csv')
	distanceMatrix = makeDistanceMatrix(ligands)
	firstRow = []
	for e in distanceMatrix[0]:
		firstRow.append(e)
	assert firstRow == [1,1,1,1,1,0,0,0,0,0], 'Make Distance Matrix Test Failed'
	print('Make Distance Matrix Test Passed')


def test_jaccard_index():
	ligands = read_test_ligands('ligand_information.csv')
	distanceMatrix = makeDistanceMatrix(ligands)
	cluster1 = Cluster(Centroid(ligands[0].bitstring))
	cluster2 = Cluster(Centroid(ligands[5].bitstring))
	cluster3 = Cluster(Centroid(ligands[0].bitstring))
	cluster4 = Cluster(Centroid(ligands[5].bitstring))
	for i in range(5):
		cluster1.addLigand(ligands[i])
		cluster3.addLigand(ligands[i])
	for i in range(5):
		cluster2.addLigand(ligands[i+5])
		cluster4.addLigand(ligands[i+5])
	ji1 = jaccardIndex([cluster1, cluster2], [cluster3, cluster4])
	ji2 = jaccardIndex([cluster1, cluster3], [cluster2, cluster4])
	assert ji1 == 1 and ji2 == 0, 'Jaccard Index Test Failed'
	print('Jaccard Index Test Passed')

	


