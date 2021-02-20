# Project 2 - Clustering and Drug Discovery
## Due 02/17/2021

![BuildStatus](https://github.com/ucsf-bmi-203-2021/Project2/workflows/HW2/badge.svg?event=push)

In this assignment, you will evaluate results from a high-throughput virtual screen against the SARS-CoV2 Spike protein / Human ACE2 interface.  There are two parts to this assignment and Part 2 requires completion of Part 1. We recommend reading through both Part 1 and Part 2 before beginning this assignment. 

* Part 1 - API and implementation
* Part 2 - Evaluating clustering

The data we are considering comes from [Smith and Smith, 2020](https://chemrxiv.org/articles/preprint/Repurposing_Therapeutics_for_the_Wuhan_Coronavirus_nCov-2019_Supercomputer-Based_Docking_to_the_Viral_S_Protein_and_Human_ACE2_Interface/11871402). In this study, they generated 6 Spike-Ace2 interface poses using MD simulations. They then docked ~10k small molecules against each protein conformation. Provided for you is the top (#1) pose for each ligand docked against one Spike-ACE2 interface conformation, as well as the corresponding SMILES string, AutoDock Vina score, and the “On” bits in the Extended Connectivity Fingerprint for that compound. These can all be found in ligand\_information.csv.


### main
Runs all code in align/\_\_main\_\_.py, useful for part 2
```
python -m algs
```

### testing
Testing is as simple as running
```
python -m pytest test/*
```
from the root directory of this project.


### Layout of Repo

```
-Ben_Orr_BMI203_HW2.ipynb: a Jupyter notebook containing written answers to assignment questions in prose, as well as function definitions and function
calls that were used in answering assignment questions.

-clusters/
  -algs.py: contains HierarchicalClustering, ParitionClustering, Centroid, Cluster, and Ligand class definitions, as well as read_ligand_csv, silhouetteCoeff,
  makeDistanceMatrix, and jaccardIndex function definitions.  Lastly, contains unit test definitions for each of the functions above, as well as for the .cluster
  methods for the two Clustering classes.

-test/
  -test_clusters.py: contains function definitions for unit tests
  
-ligand_information.csv: a .csv file containing information about ligands that were docked against the SARS-CoV2 Spike protein / Human ACE2 Interface
```

### API

```
Class HierarchicalClustering():

  - Attributes:
  
    - None


  - Methods:
  
    - cluster(self, ligands, distanceMatrix, numClusters, verbose = False):
    
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
```

```
          
Class PartitionClustering():

  - Attributes:
  
    -  None


  - Methods:

    - cluster(self, ligands, distanceMatrix, numClusters, verbose=False):

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


    - initializeClusters(self, centroids):

      For each Centroid() object in centroids, creates a new Cluster() object with the Centroid.  Returns a list of these Cluster objects

        Params:
          centroids - a list of Centroid() objects

        Returns:
          a list of Cluster() objects
          
          
    - initializeCentroids(self, numClusters, ligands, distanceMatrix):
		
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
```

```
class Centroid():

  - Attributes:
    
    - bitstring - a 1024-element array of 1's and 0's
    - onbits - a list of the indeces in the bitstring that contain a '1'


  - Methods:
    
    - bitstringToOnbits(self, bitstring):

      Calculates and returns the densified representation of the bitstring by returning a list of indeces at
        which the bitstring had 1's.

        Params:
          bitstring - a 1024 element list of 1's and 0's

        Returns:
          a list of indeces at which the bitstring contained 1's

```

