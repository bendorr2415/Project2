import numpy as np
import random

# sklearn has good pca, import that

class HierarchicalClustering():
	def __init__(self, ligands, numClusters):
		self.ligands = ligands
		self.numClusters = numClusters

	def cluster(self, ligands, numClusters):
		distanceMatrix = self.makeDistanceMatrix(ligands)

		matrixHeader = list(map(lambda ligand:[ligand], ligands)) #returns a list of lists

		while len(matrixHeader) > numClusters: # changing this from len(distanceMatrix) to len(matrixHeader) seems to have solved an error where the distance matrix stopped getting smaller once it was all zeros, but I hope there isn't an issue with the all zero distance matrix.
			
			# Determine the current max tanimoto coefficient in the distance matrix (which determines the next clustering event)
			maxVal = 0
			maxValIndeces = [0,0]
			for i in range(len(distanceMatrix)):
				for j in range(len(distanceMatrix[i])):
					if maxVal < distanceMatrix[i][j]:
							maxVal = distanceMatrix[i][j]
							maxValIndeces = [i,j]

			# Add a new column to the matrix header, which represents the combination of the two combined clusters
			matrixHeader.append(matrixHeader[maxValIndeces[0]] + matrixHeader[maxValIndeces[1]]) #Don't delete this line when deleting test code...

			# Remove the two individual clusters that were just combined from the matrix header
			matrixHeader.pop(max(maxValIndeces))
			matrixHeader.pop(min(maxValIndeces))

			# Populate an array with the min tanimoto values from each of the combined clusters (so comparing to least similar member of cluster - looks better with 40 examples)
			# I wonder if any -1's make it into the combined rows...

			# Populate the new row (and column) of the distance matrix with the minimum values from each of the combined row's indeces
			# Do not include -1's in the combined row (these correspond to a ligand's distance to itself)
			newRow = []
			for k in range(len(distanceMatrix[maxValIndeces[0]])):
				if distanceMatrix[maxValIndeces[0]][k] < distanceMatrix[maxValIndeces[1]][k] and distanceMatrix[maxValIndeces[0]][k] > 0: #so as not to include -1's
					newRow.append([distanceMatrix[maxValIndeces[0]][k]])
				else:
					newRow.append([distanceMatrix[maxValIndeces[1]][k]])
			newRow = np.array(newRow)

			# Add the new row to the distance matrix, in both the rows and columns
			distanceMatrix = np.hstack((distanceMatrix, newRow))
			newRow = np.append(newRow, [0])
			distanceMatrix = np.vstack((distanceMatrix, newRow))

			# Remove the rows and columns that were just clustered together from the distance matrix
			distanceMatrix = np.delete(distanceMatrix, [maxValIndeces[0], maxValIndeces[1]], axis = 0)
			distanceMatrix = np.delete(distanceMatrix, [maxValIndeces[0], maxValIndeces[1]], axis = 1)


		### TESTING CODE
		# Let's see what clusters are in the matrix header
		print('Lets see whats in the matrix Header;;;;;;;;;;;;;;;;;;;;;;;;;;')
		for cluster in matrixHeader:
			lst = []
			print('length of cluster = %d' % len(cluster))
			for ligand in cluster:
				lst.append(ligand.ligandID)
			print(lst)
		###END TESTING CODE


	def makeDistanceMatrix(self, ligands):
		matrix = []
		for i in range(len(ligands)):
			row = []
			for j in range(len(ligands)):
				if i == j:
					row.append(-1) # Ideally, I'd skip these values
				else:
					row.append(ligands[i].tanimotoCoefficient(ligands[j]))
			matrix.append(row)
		return np.array(matrix)



class PartitionClustering():
	def __init__(self, ligands, numClusters):
		self.ligands = ligands
		self.numClusters = numClusters

	def cluster(self, ligands, numClusters):
		centroids = self.initializeCentroids(numClusters, ligands)

		for centroid in centroids:
			print(centroid.bitstring)

		clusters = self.initializeClusters(centroids) # clusters is a list of cluster objects with centroid attributes

		# centroidUpdate is the sum of each cluster's calcCentroidMoved, which is 0 when a centroid doesn't move
		centroidUpdate = 1
		iterations = 0
		while centroidUpdate != 0 and iterations < 30:
			centroidUpdate = 0
			iterations += 1

			# Empty the ligand list in each cluster. They'll be repopulated with the new centroids
			for cluster in clusters:
				cluster.resetLigands()

			# Determine which cluster each ligand should join
			for ligand in ligands:
				tanCo = ligand.tanimotoCoefficient(clusters[0].centroid) # this works because centroid has .bitstring attribute
				i, joinClusterIndex = 0, 0
				for cluster in clusters:
					if tanCo > ligand.tanimotoCoefficient(cluster.centroid):
						tanCo = ligand.tanimotoCoefficient(cluster.centroid)
						joinClusterIndex = i
					i += 1
				clusters[joinClusterIndex].addLigand(ligand)

			# Determine whether the centroids moved
			for cluster in clusters:
				centroidUpdate += cluster.calcCentroidMoved()
				cluster.updateCentroid()


			##TESTING CODE
			print('centroid update = %d' % centroidUpdate)
			for cluster in clusters:
				lst = []
				for ligand in cluster.ligands:
					lst.append(ligand.ligandID)
				print(cluster.centroid.bitstring)
				print('sum of centroid bitstring above = %d' % sum(cluster.centroid.bitstring))
				print(lst)
				print('end of cluster')
			print('centroid update = %d' % centroidUpdate)
			##END TESTING CODE

		
		for cluster in clusters:
			lst = []
			for ligand in cluster.ligands:
				lst.append(ligand.ligandID)
			print(cluster.centroid.bitstring)
			print(lst)
			print('end of cluster')
			# 2.7.21 2:42pm:  It looks like this clustered resonably.  This code in the cluster function is pretty ugly though. There might be a prettier way to break it up
			# We'll see how reasonable this clustering looks with the real distance function, too
		

	def initializeClusters(self, centroids):
		clusters = []
		for centroid in centroids:
			newCluster = Cluster(centroid)
			clusters.append(newCluster)
		return clusters


	def initializeCentroids(self, numClusters, ligands):
		centroids = []
		indeces = random.sample(range(0,len(ligands)), numClusters)
		for i in indeces:
			centroids.append(Centroid(ligands[i].bitstring))

			print('ligand ID of centroid = %d' % ligands[i].ligandID) #this works

		return centroids


class Centroid():
	def __init__(self, bitstring):
		self.bitstring = bitstring


class Cluster():
	def __init__(self, centroid):
		self.centroid = centroid
		self.ligands = []

	def set_centroid(self, centroid):
		self.centroid = centroid

	def addLigand(self, ligand):
		self.ligands.append(ligand)

	def resetLigands(self):
		self.ligands = []

	def calcCentroid(self, ligands):
		"""
		Returns a bitstring that represents the average bitstring across all ligands in the cluster.
		Any position whose average is >= 0.5 will have a 1 in the final bitstring.
		"""
		if len(ligands) > 0:
			centroidBitstring = np.zeros(1024, dtype=int)
			for ligand in ligands:
				centroidBitstring += ligand.bitstring
			centroidBitstring = centroidBitstring / len(ligands) + 0.1 # adding this 0.1 helped prevent just two large clusters
			centroidBitstring = np.rint(centroidBitstring)
			#print('centroid bitstring = ')
			#print(centroidBitstring)
			return Centroid(centroidBitstring)
		else:
			return self.centroid


	def calcCentroidMoved(self):
		"""
		Returns the distance the centroid moved
		"""
		j = 0
		newBitstring = self.calcCentroid(self.ligands).bitstring
		for i in range(len(self.centroid.bitstring)):
			if self.centroid.bitstring[i] != newBitstring[i]:
				j+=1
		return j


	def updateCentroid(self):
		self.centroid = self.calcCentroid(self.ligands)


class Ligand():
	def __init__(self, ligandID, score, smiles, onbits):
		self.ligandID = int(ligandID)
		self.score = float(score)
		self.smiles = smiles
		self.onbits = onbits
		self.bitstring = self.onbitsToBitstring(onbits)

		#fake distance metric
		self.x = int(ligandID)
		self.y = int(ligandID) - 1

	def __str__(self):
		return "%d" % self.ligandID

	def onbitsToBitstring(self, onbits):
		onbits = onbits.split(',')
		bitstring = np.zeros(1024, dtype=int)
		for ob in onbits:
			bitstring[int(ob)] = 1
		return bitstring


	### 2.9.21 There's an issue with how I'm calculating the distance between a ligand and the centroid in partition clustering.
	### I'm assuming both will have ones and zeros, but the centroid is full of FLOATS!!!
	def tanimotoCoefficient(self, other):
		numerator, denomenator = 0, 0
		for i in range(len(self.bitstring)):
			if self.bitstring[i] == 1 and self.bitstring[i] == other.bitstring[i]:
				numerator += 1
				denomenator += 1
			elif self.bitstring[i] == 1 or other.bitstring[i] == 1:
				denomenator += 1
		return numerator / denomenator



def read_ligand_csv(csv):
	ligands = []
	i = 0 #gotta be a better way to do this
	with open(csv,'r') as f:
		for line in f:
			if i>0:
				ligandID, score, smiles, onbits = line.split(',',3) #set max number of splits = 3
				ligands.append(Ligand(ligandID, score, smiles, onbits.replace('"', '').rstrip())) #remove ""'s and eol char
			i+=1
	return ligands



csv = '../ligand_information.csv'

ligands = read_ligand_csv(csv)
#print(ligands[0].score)
#print(len(ligands))


# I need to figure out how to represent the centroids using this method
pc = PartitionClustering(ligands, 5)
pc.cluster(ligands, 5)


"""
hc = HierarchicalClustering(ligands, 5)
hc.cluster(ligands[:50], 5) #Need to limit num ligands here, since this function takes in ligands as a param
"""





#######
## UNIT TESTS
#######


def test_csv_io():
	ligands = read_ligand_csv('../ligand_information.csv')
	assert ligands[0].ligandID == 0
	assert ligands[0].score == -1.3
	print("'CSV File I/O' Test Passed")

#test_csv_io()
