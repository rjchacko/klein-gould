/*
 * Open Source Physics software is free software as described near the bottom of this code file.
 *
 * For additional information and documentation on Open Source Physics please see: 
 * <http://www.opensourcephysics.org/>
 */

package chris.ofc;


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
public class ClustersV3 {
  private static final int EMPTY = Integer.MIN_VALUE; 
  private static final int BIG = Integer.MAX_VALUE;
  public int L;                                      
  public int N;                                     
  public int numSitesOccupied;                     
  public int[] numClusters; 
  
  private int[] cnARRAY, pcn, cDist, refSite;
  private int ClCt;

  private int[] parent;
  
  private boolean PERIODIC; 
  

  public ClustersV3(int L, String BC) {
    this.L = L;
    N = L*L;
    numClusters = new int[N+1];
    parent = new int[N];

    ClCt = 1;
    cnARRAY  = new int[N];	// takes a site
    pcn = new int[N];	// takes a cn
    cDist = new int[2*N]; // takes a site ( CDist[site] = x ; CDist[site+N] = y )
    refSite = new int[N]; // takes a pcn
    
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
      cnARRAY[s] = BIG;
      pcn[s] = BIG;
      cDist[s] = 0;
      cDist[s+N] = 0;
      refSite[s] = BIG;
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
		  NBScn[i] = cnARRAY[nbs[i]];
	  }
	  
	  // sort the nbs cn
	  
	  int[] order = sortVector(NBScn,"ID");
	  
	  if (NBScn[order[0]] == BIG){	// all neighbor sites are EMPTY
		  // give the site a new cluster number and set it as a pcn
		  cnARRAY[site] = ClCt;
		  pcn[cnARRAY[site]] = ClCt;
		  
		  // now set this site as the reference site to measure distance from
		  // and set its x and y distance to 0
		  refSite[ClCt++] = site;
		  cDist[site] = 0;
		  cDist[site+N] = 0;
		  
	  }
	  else if (NBScn[order[1]] == BIG){	// only one neighbor belongs to a cluster
		  // give the site the cluster number of its neighbor
		  cnARRAY[site] = pcn[NBScn[order[0]]];
		  // the pcn for this cn has already been determined
		  
		  // set the distance this site is from the reference point
		  setCdist2(site,order[0]);
		  
	  }
	  else {	// more than one neighbor belongs to a cluster
		  
		  // call a function to deal with cluster distances
		  clusterDistance(site,nbs);
		  
		  
		  // determine which of the cn has the smallest pcn
		  int Sorg = nbs[order[0]];
		  for (int jj = 0 ; jj < 4 ; jj++){
			  if (pcn[cnMETHOD(nbs[jj])] < pcn[cnARRAY[Sorg]]) Sorg = nbs[jj];
		  }		  
		  
		  // so Sorg belongs to the cluster with the smallest pcn
		  // set up the new failed site
		  int SorgPCN = pcn[cnARRAY[Sorg]];
		  cnARRAY[site] = SorgPCN;
		  
		  // Initialize array to fix cn and pcn
		  
		  int[] temp = new int[3];
		  int tempindex = 0;

		  for (int jj = 0 ; jj < 4 ; jj++){
			  if ( (nbs[jj] != Sorg) && (cnARRAY[nbs[jj]] != BIG) ){
			  }
		  }
		  
		  
		  // loop over the neighbors and fix their pcn
		  
		  for (int jj = 0 ; jj < 4 ; jj++){
			  if ( (nbs[jj] != Sorg) && (cnARRAY[nbs[jj]] != BIG) ){
				  temp[tempindex++] = pcn[cnARRAY[nbs[jj]]];	// record the old pcn so it can be updated in the next step
				  pcn[cnARRAY[nbs[jj]]] = SorgPCN;
			  }
		  }
		  
		  // loop over the array and fix all sites with pcn equal to one of the 
		  // pcn adjusted in the above loop
		 
		  if(tempindex == 1){
			  for(int kk = 0 ; kk < N ; kk++){
				  if(pcn[cnMETHOD(kk)] == temp[0]) pcn[cnARRAY[kk]] = SorgPCN;	// fix the pcn
			  }
		  }
		  else if(tempindex == 2){
			  for(int kk = 0 ; kk < N ; kk++){
				  if((pcn[cnMETHOD(kk)] == temp[0]) || (pcn[cnMETHOD(kk)] == temp[1]) ) pcn[cnARRAY[kk]] = SorgPCN;	// fix the pcn
			  }
		  }
		  else { //tempindex ==3
			  for(int kk = 0 ; kk < N ; kk++){
				  if((pcn[cnMETHOD(kk)] == temp[0]) || (pcn[cnMETHOD(kk)] == temp[1]) || (pcn[cnMETHOD(kk)] == temp[2])) pcn[cnARRAY[kk]] = SorgPCN;	// fix the pcn
			  }
		  }
	  }
	  
	  return;
  }
 
  public void clusterDistance(int ns, int[] oss){
	  
	  // first find out how many clusters we are dealing with
	  int temp[] = new int[4];
	  int tempIndex = 0;
	  
	  for (int kk = 0 ; kk < 4 ; kk++){
		  if(oss[kk] != EMPTY) temp[tempIndex++]=oss[kk];
	  }
	  
	  // now check and see if any of the clusters are the same
	  // next, if any of the clusters are the same, check the 
	  // delta x and delta y
	  
	  // on the other hand, for clusters that are not the same
	  // it will be necessary to relabel the distances from reference
	  // as the reference site will change as the two clusters
	  // are merged
	  
	  
	  switch(tempIndex){
	  
	  	case 2:	// there are two clusters

	  		// check to see if the two clusters are the same
	  		if ( pcn[cnARRAY[temp[0]]] == pcn[cnARRAY[temp[1]]] ){
	  			// the two clusters are the same
	  		}
	  		else{
	  			// the two clusters are different
	  		}
	  			  		
	  		break;
	  		
	  	case 3:	// there are three clusters
	  		
	  		if ( (pcn[cnARRAY[temp[0]]] == pcn[cnARRAY[temp[1]]]) ){
	  			if ( (pcn[cnARRAY[temp[0]]] == pcn[cnARRAY[temp[2]]]) ) {
	  				// all three are the same
	  			}
	  			else{
	  				// 0 and 1 are the same and 2 is unique
	  			}
	  		}
	  		else{
	  			if ( (pcn[cnARRAY[temp[0]]] == pcn[cnARRAY[temp[2]]]) ) {
	  				// 0 and 2 are the same and 1 is unique
	  			}
	  			else if ( (pcn[cnARRAY[temp[1]]] == pcn[cnARRAY[temp[2]]]) ){
	  				// 1 and 2 are the same and 0 is unique
	  			}
	  			else {
	  				// all three are unique
	  			}
	  		}
	  		
	  		break;
	  		
	  	case 4: // there are four clusters

	  		// sort the clusters by pcn
	  		int[] pcnSORTindex = sortVector(oss,"ID");
	  		
	  		if( pcn[cnARRAY[oss[pcnSORTindex[0]]]]  ==  pcn[cnARRAY[oss[pcnSORTindex[1]]]] ){
	  			// 0 and 1 are the same cluster
	  			if( pcn[cnARRAY[oss[pcnSORTindex[1]]]]  ==  pcn[cnARRAY[oss[pcnSORTindex[2]]]] ){
	  				// 0 1 and 2 are the same 
	  				if( pcn[cnARRAY[oss[pcnSORTindex[2]]]]  ==  pcn[cnARRAY[oss[pcnSORTindex[3]]]] ){
	  					/* 
	  					 * all clusters are the same
	  					 */
	  				}
	  				else{
	  					/*
	  					 * 0 1 and 2 are the same but 3 is unique
	  					 */
	  				}
	  			}
	  			else{
	  				// 0 and 1 are the same and 2 and 3 are not part of ~this~ cluster
	  				if( pcn[cnARRAY[oss[pcnSORTindex[2]]]]  ==  pcn[cnARRAY[oss[pcnSORTindex[3]]]] ){
	  					/*
	  					 * 0 and 1 are the same and 2 and 3 are the same (but 0 and 2 are unique)
	  					 */
	  				}
	  				else{
	  					/*
	  					 * 0 and 1 are the same and 2 is unique and 3 is unique 
	  					 */
	  				}
	  			}
	  		}
	  		else{
	  			// 0 is a unique cluster
	  			if( pcn[cnARRAY[oss[pcnSORTindex[1]]]]  ==  pcn[cnARRAY[oss[pcnSORTindex[2]]]] ){
	  				// 0 is unique and 1 and 2 are the same
	  				if( pcn[cnARRAY[oss[pcnSORTindex[2]]]]  ==  pcn[cnARRAY[oss[pcnSORTindex[3]]]] ){
	  					/*
	  					 * 0 is unique and 1 2 and 3 are the same
	  					 */
	  				}
	  				else{	
	  					/*
	  					 * 0 is unique and 1 and 2 are the same and 3 is unique
	  					 */
	  				}
	  			}
	  			else{
	  				// 0 is unique and 1 is unique
	  				if( pcn[cnARRAY[oss[pcnSORTindex[2]]]]  ==  pcn[cnARRAY[oss[pcnSORTindex[3]]]] ){
	  					/*
	  					 * 0 is unique and 1 is unique and 3 and 3 are the same
	  					 */
	  				}
	  				else{
	  					/*
	  					 * all clusters are unique
	  					 */
	  				}
	  			}
	  		}
	  		
	  		break;
	  
	  	default:
	  		System.err.println("Internal Inconsistancy");
	  		break;
	  
	  }
	  
	  
	  
	  return;
  }
  
  private void setCdist2(int ns, int os){
	  
	  // determine where os is relative to ns
	  
	  int delta = os - ns;
	  
	  if ( (delta == -1) || (delta == 1-L) ){	// ns is to the left of os
		  cDist[ns] = cDist[os] - 1;
		  cDist[ns+N] = cDist[os+N];
	  }
	  else if ( (delta == 1) || (delta == L-1) ) { // ns is to the right of os
		  cDist[ns] = cDist[os] + 1;
		  cDist[ns+N] = cDist[os+N];
	  }
	  else if ( (delta == -L) || (delta == N-L-1) ){  // ns is above os
		  cDist[ns] = cDist[os];
		  cDist[ns+N] = cDist[os+N] + 1;
	  }
	  else if ( (delta == L) || (delta == L+1-N) ){  // ns is below os
		  cDist[ns] = cDist[os];
		  cDist[ns+N] = cDist[os+N] - 1;
	  }
	  else{
		 System.err.println("Internal Inconsistancy"); 
	  }
	  
	  return;
  }
  
  private int cnMETHOD(int st){
	  return (cnARRAY[st] == BIG) ? N-1 : cnARRAY[st];
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
	return(cnARRAY[s] == BIG) ? 0 : pcn[cnARRAY[s]];
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