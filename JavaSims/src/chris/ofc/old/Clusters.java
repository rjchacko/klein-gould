/*
 * Open Source Physics software is free software as described near the bottom of this code file.
 *
 * For additional information and documentation on Open Source Physics please see: 
 * <http://www.opensourcephysics.org/>
 */

package chris.ofc.old;


/**
 *  Clusters implements the Newman-Ziff algorithm for identifying clusters.
 *
 *  @author Jan Tobochnik, Wolfgang Christian, Harvey Gould
 *  @version 1.0 06/15/05
 *  
 *  @revised 04/09/08 
 *  @author Christopher A. Serino
 *  
 */
public class Clusters {
  private static final int EMPTY = Integer.MIN_VALUE; 
  private static final int BIG = Integer.MAX_VALUE;
  public int L;                                      
  public int N;                                     
  public int numSitesOccupied;                     
  public int[] numClusters; 
  
  private int[] cn, pcn; 
  private int ClCt;

  private int[] parent;
  
  private boolean PERIODIC; 
  

  public Clusters(int L, String BC) {
    this.L = L;
    N = L*L;
    numClusters = new int[N+1];
    parent = new int[N];

    ClCt = 1;
    cn  = new int[N];	// takes a site
    pcn = new int[N];	// takes a cn
    
    if (BC.equals("Periodic")){
    	PERIODIC = true;
    }
    else{
    	PERIODIC = false;
    }
    
    
    newLattice();
    
  }

  public void newLattice() {

    numSitesOccupied = 0;
    for(int s = 0;s<N;s++) {
      numClusters[s] = 0;
      parent[s] = EMPTY;
      cn[s] = BIG;
      pcn[s] = cn[s];
    }
  }

  public boolean addSite(int site){
	  
	  	int[] nbs = new int[4];
	  
	    numClusters[1]++;
	    parent[site] = -1;
	    int root = site;
	    for(int j = 0;j<4;j++) {		
	      int neighborSite = getNeighbor(site, j);
	      nbs[j] = neighborSite;
	      if((neighborSite!=EMPTY)&&(parent[neighborSite]!=EMPTY)) {
	        root = mergeRoots(root, findRoot(neighborSite));
	      }
	    }
	    
	    setClusterNumber(site,nbs);
	    
	    // did adding this site cause percolation?
	    return false;
  }

  private void setClusterNumber(int site, int[] nbs){
	  
	  // get the cn of the nbs
	  
	  int[] NBScn = new int[4];
	  
	  for (int i = 0 ; i < 4 ; i++){
		  NBScn[i] = cn[nbs[i]];
	  }
	  
	  // sort the nbs cn
	  
	  int[] order = sortVector(NBScn,"ID");
	  
	  if (NBScn[order[0]] == BIG){	// all neighbor sites are EMPTY
		  // give the site a new cluster number and set it as a pcn
		  cn[site] = ClCt++;
		  pcn[cn[site]] = cn[site];
		  
	  }
	  else if (NBScn[order[1]] == BIG){	// only one neighbor belongs to a cluster
		  // give the site the cluster number of its neighbor
		  cn[site] = pcn[NBScn[order[0]]];
		  // the pcn for this cn has already been determined
	  }
	  else {	// more than one neighbor belongs to a cluster

		  // determine which of the cn has the smallest pcn
		  int Sorg = nbs[order[0]];
		  for (int jj = 0 ; jj < 4 ; jj++){
			  if (cn[nbs[jj]] != BIG){
				  if (pcn[cn[nbs[jj]]] < pcn[cn[Sorg]]) Sorg = nbs[jj];
			  }
		  }
		  // so Sorg belongs to the cluster with the smallest pcn
		  // set up the new failed site
		  int SorgPCN = pcn[cn[Sorg]];
		  cn[site] = SorgPCN;
		  
		  // Initialize array to fix cn and pcn
		  
		  int[] temp = new int[3];
		  int tempindex = 0;
		  
		  // loop over the neighbors and fix their pcn
		  
		  for (int jj = 0 ; jj < 4 ; jj++){
			  if ( (nbs[jj] != Sorg) && (cn[nbs[jj]] != BIG) ){
				  temp[tempindex++] = pcn[cn[nbs[jj]]];	// record the old pcn as it is no longer a pcn
				  pcn[cn[nbs[jj]]] = SorgPCN;
			  }
		  }
		  
		  // loop over the array and fix all sites with pcn equal to one of the 
		  // pcn adjusted in the above loop
		  
		  if(tempindex == 1){
			  for(int kk = 0 ; kk < N ; kk++){
				  if(cn[kk] != BIG){
					  if(pcn[cn[kk]] == temp[0]) pcn[cn[kk]] = SorgPCN;	// fix the pcn
				  }
			  }
		  }
		  else if(tempindex == 2){
			  for(int kk = 0 ; kk < N ; kk++){
				  if(cn[kk] != BIG){
					  if((pcn[cn[kk]] == temp[0]) || (pcn[cn[kk]] == temp[1]) ) pcn[cn[kk]] = SorgPCN;	// fix the pcn
				  }
			  }
		  }
		  else { //tempindex ==3
			  for(int kk = 0 ; kk < N ; kk++){
				  if(cn[kk] != BIG){
					  if((pcn[cn[kk]] == temp[0]) || (pcn[cn[kk]] == temp[1]) || (pcn[cn[kk]] == temp[2])) pcn[cn[kk]] = SorgPCN;	// fix the pcn
				  }
			  }
		  }		  
	  }
	  
	  return;
  }
  
  public int[] sortVector(int[] input, String idORval){

	  int LL = input.length;
	  int k0;
	  int[][] copy = new int[LL][2];
	  int[] ret    = new int[LL];
	  int[] sorted = new int[LL];
	  
	  for (int copyINDEX = 0 ; copyINDEX < LL ; copyINDEX++){
		  copy[copyINDEX][0] = input[copyINDEX];
		  copy[copyINDEX][1] = copyINDEX;
	  }
	  
	  for (int jj = 0 ; jj < LL ; jj++){
		  k0=jj;
		  for (int kk = jj ; kk < LL ; kk++){
			  if(copy[kk][0] < copy[k0][0]) k0 = kk;
		  }
		  int temp1 = copy[k0][0];
		  int temp2 = copy[k0][1];
		  copy[k0][0]=copy[jj][0];
		  copy[k0][1]=copy[jj][1];
		  copy[jj][0] = temp1;
		  copy[jj][1] = temp2;
	  }
	  
	  for (int retINDEX = 0 ; retINDEX < LL ; retINDEX++){
		  ret[retINDEX]    = copy[retINDEX][1];
		  sorted[retINDEX] = copy[retINDEX][0];
	  }
	  
	  if (idORval.equals("ID")){
		  return ret;
	  }
	  else if (idORval.equals("value")){
		  return sorted;
	  }
	  else{
		  return null;
	  }
	  
  }
  
  public int getClusterSize(int s) {
    return(parent[s]==EMPTY) ? 0 : -parent[findRoot(s)];
  }

  public int getClusterNumber(int s){
	return(cn[s] == BIG) ? 0 : pcn[cn[s]];
  }
  
  
  private int findRoot(int s) {
    if(parent[s]<0) {
      return s;
    } else {
      parent[s] = findRoot(parent[s]);
    }
    return parent[s];
  }

  private int getNeighbor(int s, int j) {
	  
	if (PERIODIC){
		
		  switch(j) {
		    case 0 :
		      return(s%L==0) ? s+L-1 : s-1;   // left
		    case 1 :
		      return(s%L==L-1) ? s-L+1 : s+1; // right
		    case 2 :
		      return(s/L==0) ? L*(L-1)+s : s-L;   // down
		    case 3 :
		      return(s/L==L-1) ? s-(L-1)*L : s+L; // above
		    default :
		      return EMPTY;
		    }
		  
	}
	else{
		
    switch(j) {
	    case 0 :
	      return(s%L==0) ? EMPTY : s-1;   // left
	    case 1 :
	      return(s%L==L-1) ? EMPTY : s+1; // right
	    case 2 :
	      return(s/L==0) ? EMPTY : s-L;   // down
	    case 3 :
	      return(s/L==L-1) ? EMPTY : s+L; // above
	    default :
	      return EMPTY;
	    }
	
	}
  }
  
  private int mergeRoots(int r1, int r2) {

    if(r1==r2) {
      return r1;
    } else if(-parent[r1]<-parent[r2]) {
      return mergeRoots(r2, r1);
    } else { // (-parent[r1] > -parent[r2])
      
      numClusters[-parent[r1]]--;
      numClusters[-parent[r2]]--;
      numClusters[-parent[r1]-parent[r2]]++;
      parent[r1] += parent[r2];
      parent[r2] = r1;

      return r1;
    }
  }
  

  
}