package ranjit.lrising;

import java.util.Random;


public class lrClusters{
  static private final int NONE = Integer.MIN_VALUE;  	
  int L,R,N;                        // linear dimension of lattice
 
  double bondProbability;       // probability of a bond between sites
  
  int nuclei_size,nuclei_cm;				// Root of largest cluster
  double[] nuclei_profile;
  
  public int[] numClusters;     // number of clusters of size s (n_s)
  public int numSitesOccupied;         // number of occupied lattice sites 

  // secondClusterMoment stores the quantity sum{s^2 n_s}, where the sum is over
  // all s (not counting the spanning cluster)
  // n_s is the number of clusters of size s.  The first
  // cluster moment, sum{s n_s} is equal to numSitesOccupied.
  // A quantity of interest is the "mean cluster size" S, which is defined as
  // S = secondClusterMoment / numSitesOccupied
  public int secondClusterMoment;
  public int massLargestCluster=1;
  
  // the parent[] array serves two purposes: It stores the cluster size
  // when the site is the (unique) root.  Otherwise, it stores the index of
  // the site's "parent." We can find the root from an arbitrary occupied site
  // by recursively following the parent array. This recursion terminates
  // when we encounter a negative value in the parent array, which indicates we have found
  // the unique cluster root.
  //    if (parent[s] >= 0)          parent[s] is the parent site index
  //    if (0 > parent[s] > NONE)    s is the root of size -parent[s]
  //    if (parent[s] == NONE)       site s is vacant
  private int[] parent;
  Random r= new Random();
  
  public lrClusters(int L, double bondProbability) {
    this.L = L;
    this.N=L*L;
    this.bondProbability = bondProbability;
    numClusters = new int[N+1];
    parent      = new int[N];
    
    nuclei_profile = new double[N];
    
  }
  
  public void newLattice() {
    // initially all sites are empty, and there are no clusters
    numSitesOccupied = secondClusterMoment = 0;
    for (int s = 0; s < N; s++) {
      numClusters[s] = 0;
      parent[s] = NONE;
    }
    numClusters[N] = 0;
    massLargestCluster=1;
  }
  
  public void newNuclei(){
  	for(int i = 0; i < N; i++)nuclei_profile[i] = 0.0;
  }
  
  // add a site to the lattice, and update the clusters.
  public void addSite(int newSite) {
    if (parent[newSite] != NONE) return;
    // create a new cluster containing only the site newSite.
    numClusters[1]++;
    secondClusterMoment++;
    numSitesOccupied++;
    
    // store the new cluster's size in parent[].  The negative sign
    // distinguishes newSite as a root, with a size value.  (Positive values
    // correspond to non-root sites with index pointers).
    parent[newSite] = -1;
    // Merge newSite with occupied neighbors.  'root' is the 
    // index of the merged cluster root at each step
    int root = newSite;
    // n[j] is the jth site neighboring the newly added site newSite
    //pick random sites to bond to
    for(int i=0;i<parent.length;i++){
      if (i!=root && parent[i] != NONE && r.nextDouble() < bondProbability){
        root = mergeRoots(root, findRoot(i));	
      }
    }
  }
  
  // get size of cluster to which site s belongs.
  public int getClusterSize(int s) {
    return parent[s] == NONE ? 0 : -parent[findRoot(s)];
  }
  
  // returns S ("mean cluster size")
  public double getMeanClusterSize() {
  	int numclusters;
  	numclusters = 0;
  for(int i = 1; i < N + 1; i++)numclusters += numClusters[i];
  	
    if (numSitesOccupied > 0)
      return (double)numSitesOccupied/(double)numclusters;
    else
      return 0;
  }
  
  public int getMaxSize(){	/*Size that has most # of clusters*/
  	int max = 0;
  	for(int i = 0; i < N + 1; i++){
  		if(numClusters[i] > max)max = numClusters[i];
  	}
  	return max;
  }
  
  public int getNucleiStat(){
  	int max = 0;
  	int nuclei_root = -1;
  	
  	for(int i = 0; i < N; i++){
  		if(-parent[i] > max){
  			max = -parent[i];
  			nuclei_root = i;
  		}
  	}
  	
  	/*Getting density profile for nuclei*/
  	for(int i = 0; i < N; i++){
  		if(parent[i] != NONE){
  			if(findRoot(i) == nuclei_root){
  				nuclei_profile[i] += 1.0;

  				//System.out.println("Adding site\t" + i + "\t" + nuclei_profile[i]);
  			}
  		}
  	}
  	
  	nuclei_size += max;
  	return max;
  	
  }
  
  public void getCM(){
  	double mass = 0;
  	double cm = 0;
  	for(int i = 0; i < N; i++){
  		mass += nuclei_profile[i];
  		cm += i*nuclei_profile[i];
  	}
  	cm /= mass;
  	nuclei_cm = (int)cm;
  }
  
  public void clear_nuclei(){
  	newNuclei();
  	nuclei_size = 0;
  	nuclei_cm = 0;
  }
  // given a site index s, return the site index representing the root
  // of the cluster to which s belongs.
  private int findRoot(int s) {
    if (parent[s] < 0)
      return s; // i is the root site (with size -parent[s])
    else
      // first link parent[s] to the cluster's root to improve performance
      // (path compression.)  then return this value.
      return parent[s] = findRoot(parent[s]);
  }
  
  // utility method to square an integer
  private int sqr(int x) {return x*x;}
  
  public void setBondProbability(double _prob){
  	bondProbability = _prob;
  }
  
  // merge two root sites into one.  this represents cluster merging.
  // use the heuristic that the root of the smaller cluster points to
  // the root of the larger cluster, in order to improve performance.
  // remember that parent[root] stores negative cluster size.
  private int mergeRoots(int r1, int r2) {
    // clusters are uniquely identified by their root sites.  if they
    // are the same, then the clusters are already merged, and we need
    // do nothing
    if (r1 == r2)
      return r1;
    // if r1 has smaller cluster size than r2, reverse (r1,r2) labels
    else if (-parent[r1] < -parent[r2])
      return mergeRoots(r2, r1);
    else /* (-parent[r1] > -parent[r2]) */ {
      // update the cluster count, and second cluster moment to account for the
      // loss of two small clusters and the gain of one big cluster    	
      numClusters[-parent[r1]]--;
      numClusters[-parent[r2]]--;
      numClusters[-parent[r1]-parent[r2]]++;
      secondClusterMoment += sqr(parent[r1]+parent[r2]) - sqr(parent[r1]) - sqr(parent[r2]);
      
      // the cluster at r1 now includes the sites of old cluster at r2
      parent[r1] += parent[r2];
      if(-parent[r1]>massLargestCluster) massLargestCluster=-parent[r1];
      // make r1 the new parent of r2
      parent[r2] = r1;
      // return the new root site r1
      return r1;
    }
  }

}