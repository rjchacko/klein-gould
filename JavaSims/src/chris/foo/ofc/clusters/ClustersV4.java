/*
 * Open Source Physics software is free software as described near the bottom of this code file.
 *
 * For additional information and documentation on Open Source Physics please see: 
 * <http://www.opensourcephysics.org/>
 */

package chris.foo.ofc.clusters;


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
public class ClustersV4 {
  private static final int EMPTY = Integer.MIN_VALUE; 
  private static final int BIG = Integer.MAX_VALUE;
  public int L;                                      
  public int N;                                     
  public int numSitesOccupied;                     
  public int[] numClusters; 
  
  private int[] cnARRAY, pcn, cDist;
  private int ClCt, PercolationCluster;
  

  
  
  private int[] parent;
  
  private boolean PERIODIC, SystemPercolate; 
  

  public ClustersV4(int L, String BC) {

	  PseudoConstructor(L,BC);
    
  }
  
  private void PseudoConstructor(int L, String BC){
	  
	    this.L = L;
	    N = L*L;
	    numClusters = new int[N+1];
	    parent = new int[N];

	    ClCt = 1;
	    cnARRAY  = new int[N];	// takes a site
	    pcn = new int[N];	// takes a cn
	    
	    if (BC.equals("Periodic")){
	    	PERIODIC = true;
	    }
	    else{
	    	PERIODIC = false;
	    }
	    
	    cDist = new int[2*N]; // takes a site ( CDist[site] = x ; CDist[site+N] = y )
	    
	    PercolationCluster = EMPTY;
	    SystemPercolate = false;
	    
	    newLattice();
	  
	  return;
  }

  public void newLattice() {

    numSitesOccupied = 0;
    for(int s = 0;s<N;s++) {
      numClusters[s] = 0;
      parent[s] = EMPTY;
      cnARRAY[s] = BIG;
      pcn[s] = BIG;
      cDist[s] = EMPTY;
      cDist[s+N] = EMPTY;
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
	    
	    boolean tempRET = setCNandDIST(site,nbs);
	    
	    return tempRET;
  }
  
  private boolean setCNandDIST(int site, int[] nbs){
	  
	  // get the cn of the nbs
	  
	  int[] NBScn = new int[4];
	  
	  for (int i = 0 ; i < 4 ; i++){
		  NBScn[i] = cnARRAY[nbs[i]];
	  }
	  
	  // sort the nbs cn
	  
	  int[] order = sortVector(NBScn,"ID");
	  
	  // proceed according to how many of these sites are OCCUPIED
	  
	  if (NBScn[order[0]] == BIG){	// all neighbor sites are EMPTY
		  // give the site a new cluster number and set it as a pcn
		  cnARRAY[site] = ClCt;
		  pcn[cnARRAY[site]] = ClCt++;
		  // now set this site as the reference site to measure distance from
		  // and set its x and y distance to 0
		  cDist[site] = 0;
		  cDist[site+N] = 0;
		  
		  // lattice has not percolated
	  }
	  else if (NBScn[order[1]] == BIG){	// only one neighbor belongs to a cluster
		  // give the site the cluster number of its neighbor
		  cnARRAY[site] = pcn[NBScn[order[0]]];
		  // the pcn for this cn has already been determined
		  
		  // set the distance this site is from the reference point
		  setCdist2(site,nbs[order[0]]);
		  
		  // lattice has not percolated
	  }
	  else{ // more than one neighbor site is OCCUPIED and thus the lattice may have percolated
		  
		  // first find out how many clusters we are dealing with
		  int OCCnbs[] = new int[4];
		  int numNeighbors = 0;
		  
		  for (int kk = 0 ; kk < 4 ; kk++){
			  if(NBScn[kk] != BIG) OCCnbs[numNeighbors++]=nbs[kk];
		  }
		  
		  // now check and see if any of the clusters are the same
		  // next, if any of the clusters are the same, check the 
		  // delta x and delta y
		  
		  // on the other hand, for clusters that are not the same
		  // it will be necessary to relabel the distances from reference
		  // as the reference site will change as the two clusters
		  // are merged
		  
		  switch(numNeighbors){

		  case 2:	// there are two clusters

			  // check to see if the two clusters are the same
			  if ( pcn[cnARRAY[OCCnbs[0]]] == pcn[cnARRAY[OCCnbs[1]]] ){
				  /*
				   * the two clusters are the same
				   */
				  if( look4percolation(site, new int[]{OCCnbs[0],OCCnbs[1]}) ){
					  //  system has percolated
					  
					  // set the cn to PercolationCluster and use one of the OCC's to set the distance
					  cnARRAY[site] = PercolationCluster;
					  setCdist2(site,OCCnbs[0]);
				  }
				  else{
					  //  system has NOT percolated
					  
					  // since the two clusters are the same but the system has not 
					  // percolated, we can simply add the site to the cluster of 
					  // either neighbor and use that neighbor to set the distance
					  
					  cnARRAY[site] = pcn[cnARRAY[OCCnbs[0]]];
					  setCdist2(site,OCCnbs[0]);
				  }
			  }
			  else{
				  /*
				   * the two clusters are different
				   */
				  
				  // add site to the cluster with the smallest pcn
				  // merge the two clusters
				  // correct all pcn's and distances
				  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1]});
				  
			  }

			  break;

		  case 3:	// there are three clusters

			  if ( (pcn[cnARRAY[OCCnbs[0]]] == pcn[cnARRAY[OCCnbs[1]]]) ){
				  if ( (pcn[cnARRAY[OCCnbs[0]]] == pcn[cnARRAY[OCCnbs[2]]]) ) {
					  /*
					   * all three are the same
					   */
					  if( look4percolation(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2]}) ){
						  //  system has percolated
						  
						  // set the cn to PercolationCluster and use one of the OCC's to set the distance
						  cnARRAY[site] = PercolationCluster;
						  setCdist2(site,OCCnbs[0]);
						  
					  }
					  else{
						  //  system has NOT percolated
						  
						  // since the three clusters are the same but the system has not 
						  // percolated, we can simply add the site to the cluster of 
						  // either neighbor and use that neighbor to set the distance
						  
						  cnARRAY[site] = pcn[cnARRAY[OCCnbs[0]]];
						  setCdist2(site,OCCnbs[0]);
					  }
				  }
				  else{
					  /*
					   * 0 and 1 are the same and 2 is unique
					   */
					  if( look4percolation(site, new int[]{OCCnbs[0],OCCnbs[1]}) ){
						  //  system has percolated
						 
						  cnARRAY[site] = PercolationCluster;
						  setCdist2(site,OCCnbs[0]);
						  // now merge the clusters
						  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2]});
						  PercolationCluster = pcn[cnARRAY[site]];
					  }
					  else{
						  //  system has NOT percolated
						  
						  //merge the clusters
						  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2]});
					  }
				  }
			  }
			  else{
				  if ( (pcn[cnARRAY[OCCnbs[0]]] == pcn[cnARRAY[OCCnbs[2]]]) ) {
					  /* 
					   * 0 and 2 are the same and 1 is unique
					   */
					  if( look4percolation(site, new int[]{OCCnbs[0],OCCnbs[2]}) ){
						  //  system has percolated
						  
						  cnARRAY[site] = PercolationCluster;
						  setCdist2(site,OCCnbs[0]);
						  // now merge the clusters
						  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2]});
						  PercolationCluster = pcn[cnARRAY[site]];
						  
					  }
					  else{
						  //  system has NOT percolated
						  
						  //merge the clusters
						  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2]});
					  }
				  }
				  else if ( (pcn[cnARRAY[OCCnbs[1]]] == pcn[cnARRAY[OCCnbs[2]]]) ){
					  /*
					   * 1 and 2 are the same and 0 is unique
					   */ 
					  if( look4percolation(site, new int[]{OCCnbs[1],OCCnbs[2]}) ){
						  //  system has percolated
						  
						  cnARRAY[site] = PercolationCluster;
						  setCdist2(site,OCCnbs[1]);
						  // now merge the clusters
						  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2]});
						  PercolationCluster = pcn[cnARRAY[site]];
					  }
					  else{
						  //  system has NOT percolated
						  
						  //merge the clusters
						  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2]});  
					  }
				  }
				  else {
					  /*
					   * all three are unique
					   */ 
					  
					  //merge the clusters
					  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2]});
					  
				  }
			  }

			  break;

		  case 4: // there are four clusters

			  
			  // sort the clusters by pcn
			  int[] pcnSORTindex = sortVector(nbs,"ID");

			  if( pcn[cnMETHOD(nbs[pcnSORTindex[0]])]  ==  pcn[cnMETHOD(nbs[pcnSORTindex[1]])] ){
				  // 0 and 1 are the same cluster
				  if( pcn[cnMETHOD(nbs[pcnSORTindex[1]])]  ==  pcn[cnMETHOD(nbs[pcnSORTindex[2]])] ){
					  // 0 1 and 2 are the same 
					  if( pcn[cnMETHOD(nbs[pcnSORTindex[2]])]  ==  pcn[cnMETHOD(nbs[pcnSORTindex[3]])] ){
						  /* 
						   * all clusters are the same
						   */
						  if( look4percolation(site, nbs) ){
							  //  system has percolated
							  
							  // set the cn to PercolationCluster and use one of the OCC's to set the distance
							  cnARRAY[site] = PercolationCluster;
							  setCdist2(site,OCCnbs[0]);
						  }
						  else{
							  //  system has NOT percolated
							  
							  // since the four clusters are the same but the system has not 
							  // percolated, we can simply add the site to the cluster of 
							  // either neighbor and use that neighbor to set the distance
							  
							  cnARRAY[site] = pcn[cnMETHOD(OCCnbs[0])];
							  setCdist2(site,OCCnbs[0]);
						  }
					  }
					  else{
						  /*
						   * 0 1 and 2 are the same but 3 is unique
						   */
						  if( look4percolation(site, new int[]{nbs[pcnSORTindex[0]],nbs[pcnSORTindex[1]],nbs[pcnSORTindex[2]]}) ){
							  //  system has percolated
							  
							  cnARRAY[site] = PercolationCluster;
							  setCdist2(site,OCCnbs[0]);
							  // now merge the clusters
							  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2],OCCnbs[3]});
							  PercolationCluster = pcn[cnMETHOD(site)];
						  }
						  else{
							  //  system has NOT percolated
							  
							  //merge the clusters
							  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2],OCCnbs[3]});
						  }
					  }
				  }
				  else{
					  // 0 and 1 are the same and 2 and 3 are not part of ~this~ cluster
					  if( pcn[cnMETHOD(nbs[pcnSORTindex[2]])]  ==  pcn[cnMETHOD(nbs[pcnSORTindex[3]])] ){
						  /*
						   * 0 and 1 are the same and 2 and 3 are the same (but 0 and 2 are unique)
						   * 
						   * THIS is the odd-ball (so to speak)
						   * 
						   */
						  if( look4percolation(site, new int[]{nbs[pcnSORTindex[0]],nbs[pcnSORTindex[1]]}) || look4percolation(site, new int[]{nbs[pcnSORTindex[2]],nbs[pcnSORTindex[3]]}) ){
							  //  system has percolated
							  cnARRAY[site] = PercolationCluster;
							  setCdist2(site,OCCnbs[0]);
							  // now merge the clusters
							  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2],OCCnbs[3]});
							  PercolationCluster = pcn[cnMETHOD(site)];
						  }
						  else{
							  //  system has NOT percolated
							  
							  //merge the clusters
							  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2],OCCnbs[3]});
						  }

					  }
					  else{
						  /*
						   * 0 and 1 are the same and 2 is unique and 3 is unique 
						   */
						  if( look4percolation(site, new int[]{nbs[pcnSORTindex[0]],nbs[pcnSORTindex[1]]}) ){
							  //  system has percolated
							  
							  cnARRAY[site] = PercolationCluster;
							  setCdist2(site,OCCnbs[0]);
							  // now merge the clusters
							  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2],OCCnbs[3]});
							  PercolationCluster = pcn[cnMETHOD(site)];
						  }
						  else{
							  //  system has NOT percolated
							  
							  //merge the clusters
							  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2],OCCnbs[3]});
						  }
					  }
				  }
			  }
			  else{
				  // 0 is a unique cluster
				  if( pcn[cnMETHOD(nbs[pcnSORTindex[1]])]  ==  pcn[cnMETHOD(nbs[pcnSORTindex[2]])] ){
					  // 0 is unique and 1 and 2 are the same
					  if( pcn[cnMETHOD(nbs[pcnSORTindex[2]])]  ==  pcn[cnMETHOD(nbs[pcnSORTindex[3]])] ){
						  /*
						   * 0 is unique and 1 2 and 3 are the same
						   */
						  if( look4percolation(site, new int[]{nbs[pcnSORTindex[1]],nbs[pcnSORTindex[2]],nbs[pcnSORTindex[3]]}) ){
							  //  system has percolated
							  
							  cnARRAY[site] = PercolationCluster;
							  setCdist2(site,OCCnbs[1]);
							  // now merge the clusters
							  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2],OCCnbs[3]});
							  PercolationCluster = pcn[cnMETHOD(site)];
							  
						  }
						  else{
							  //  system has NOT percolated
							  
							  //merge the clusters
							  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2],OCCnbs[3]});
						  }

					  }
					  else{	
						  /*
						   * 0 is unique and 1 and 2 are the same and 3 is unique
						   */
						  if( look4percolation(site, new int[]{nbs[pcnSORTindex[1]],nbs[pcnSORTindex[2]]}) ){
							  //  system has percolated
							  
							  cnARRAY[site] = PercolationCluster;
							  setCdist2(site,OCCnbs[1]);
							  // now merge the clusters
							  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2],OCCnbs[3]});
							  PercolationCluster = pcn[cnMETHOD(site)];
						  }
						  else{
							  //  system has NOT percolated
							  
							  //merge the clusters
							  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2],OCCnbs[3]});
						  }
					  }
				  }
				  else{
					  // 0 is unique and 1 is unique
					  if( pcn[cnMETHOD(nbs[pcnSORTindex[2]])]  ==  pcn[cnMETHOD(nbs[pcnSORTindex[3]])] ){
						  /*
						   * 0 is unique and 1 is unique and 2 and 3 are the same
						   */
						  if( look4percolation(site, new int[]{nbs[pcnSORTindex[2]],nbs[pcnSORTindex[3]]}) ){
							  //  system has percolated
							  
							  cnARRAY[site] = PercolationCluster;
							  setCdist2(site,OCCnbs[2]);
							  // now merge the clusters
							  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2],OCCnbs[3]});
							  PercolationCluster = pcn[cnMETHOD(site)];
						  }
						  else{
							  //  system has NOT percolated
							  
							  //merge the clusters
							  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2],OCCnbs[3]});
						  }

					  }
					  else{
						  /*
						   * all clusters are unique
						   */
						  
						  //merge the clusters
						  resetPCNandDIST(site, new int[]{OCCnbs[0],OCCnbs[1],OCCnbs[2],OCCnbs[3]});
					  }
				  }
			  }

			  break;

		  default:
			  System.err.println("Internal Inconsistancy");
		  break;

		  }
	  }
	  
	  return SystemPercolate;
  }
  
private void resetPCNandDIST(int ns, int[] os){
	  
	  // determine which cluster has the smallest pcn
	  int Sorg = os[0];
	  for (int jj = 0 ; jj < os.length ; jj++){
		  if (pcn[cnMETHOD(os[jj])] < pcn[cnMETHOD(Sorg)]) Sorg = os[jj];
	  }
	  int SorgPCN = pcn[cnARRAY[Sorg]];
	  // so Sorg is a site index which belongs to the cluster with the smallest pcn
	  // and SorgPCN is that cluster's PCN
	  
	  // get the pcn that need to be reset and
	  // set the distances of the neighbors wrt the new site
	  
	  int miniN = os.length; 
	  int[] resetPCN = new int[miniN];
	  int[] newDIST = new int[2*miniN];
	  int[] oldDIST = new int[2*miniN];
	  
	  for (int jj = 0 ; jj < miniN ; jj++){
		  resetPCN[jj] = pcn[cnMETHOD(os[jj])];

		  oldDIST[jj] = cDist[os[jj]];
		  oldDIST[jj+miniN] = cDist[os[jj]+N];
		  
		  int delta = os[jj] - ns;
		  if ( (delta == 1) || (delta == 1-L) ){	// ns is to the left of os
			  newDIST[jj] = 1;
			  newDIST[jj+miniN] = 0;
		  }
		  else if ( (delta == -1) || (delta == L-1) ) { // ns is to the right of os
			  newDIST[jj] = -1;
			  newDIST[jj+miniN] = 0;
		  }
		  else if ( (delta == -L) || (delta == N-L) ){  // ns is above os
			  newDIST[jj] = 0;
			  newDIST[jj+miniN] = -1;
		  }
		  else if ( (delta == L) || (delta == L-N) ){  // ns is below os
			  newDIST[jj] = 0;
			  newDIST[jj+miniN] = 1;
		  }
		  else{
			 System.err.println("Internal Inconsistancy in Finding Location of Neighbor"); 
		  }
	  }
	  
	  // loop over the lattice and find sites that need to be reset
	  int[] Relabel = new int[2*N];
	  int RelabelCounter = 0;
	  for (int jj = 0 ; jj < N ; jj++){
		  for (int kk = 0 ; kk < miniN ; kk++){
			  if(pcn[cnMETHOD(jj)] == resetPCN[kk]){
				  Relabel[RelabelCounter] = jj;
				  Relabel[RelabelCounter+N] = kk;
				  RelabelCounter++;
				  break;
			  }
		  }
	  }

	  // now, make the changes to the sites found above
	  for (int jj = 0 ; jj < RelabelCounter ; jj++){
		  pcn[cnMETHOD(Relabel[jj])] = SorgPCN;
		  cDist[Relabel[jj]] += (newDIST[Relabel[jj+N]] - oldDIST[Relabel[jj+N]]);
		  //System.out.println("y_old <cr> y_new <cr> newC <cr> oldC");
		  cDist[Relabel[jj]+N] += (newDIST[Relabel[jj+N]+miniN] - oldDIST[Relabel[jj+N]+miniN]);
	  }
	  
	  // setup the ns and set the ns as the index for pcn = SorgPCN  
	  cnARRAY[ns] = SorgPCN;
	  cDist[ns] = 0;
	  cDist[ns+N]=0;
	  
	  return;
  }
  
  
  private boolean look4percolation(int ns, int[] comp){	// comp means sites to compare
	  
	  int Length = comp.length;
	  int[][] compCheck = new int [Length][2];
	  
	  for (int qq = 0 ; qq < Length ; qq++){
		  int os = comp[qq];
		  int delta = os - ns;

		  if ( (delta == 1) || (delta == 1-L) ){	// ns is to the left of os
			  compCheck[qq][0] = cDist[os] - 1;
			  compCheck[qq][1] = cDist[os+N];
		  }
		  else if ( (delta == -1) || (delta == L-1) ) { // ns is to the right of os
			  compCheck[qq][0] = cDist[os] + 1;
			  compCheck[qq][1] = cDist[os+N];
		  }
		  else if ( (delta == -L) || (delta == N-L) ){  // ns is above os
			  compCheck[qq][0] = cDist[os];
			  compCheck[qq][1] = cDist[os+N] + 1;
		  }
		  else if ( (delta == L) || (delta == L-N) ){  // ns is below os
			  compCheck[qq][0] = cDist[os];
			  compCheck[qq][1] = cDist[os+N] - 1;
		  }
		  else{
			 System.err.println("Internal Inconsistancy in Finding Location of Neighbor"); 
		  }
  
	  }
	  
	  for (int qq = 0 ; qq < Length ; qq++){
		  for (int rr = qq+1 ; rr < Length ; rr++){
			  if( ( (compCheck[qq][0] - compCheck[rr][0]) != 0 )  || ( (compCheck[qq][1] - compCheck[rr][1]) != 0 ) ){
				  PercolationCluster = pcn[cnARRAY[comp[qq]]];
			  }
		  }
		  
	  }
	  
	  if(PercolationCluster != EMPTY) SystemPercolate = true;
	  
	  return SystemPercolate;

	  // IF SystemPercolate == true THEN
	  //
	  // PercolationCluster holds the pcn of the percolating cluster 
	  // 
	  // ELSE
	  //
	  // no changes are made
	  //
	  
  }
  
  private void setCdist2(int ns, int os){
	  
	  // determine where os is relative to ns
	  
	  int delta = os - ns;
	 
	  if ( (delta == 1) || (delta == 1-L) ){	// ns is to the left of os
		  cDist[ns] = cDist[os] - 1;
		  cDist[ns+N] = cDist[os+N];
	  }
	  else if ( (delta == -1) || (delta == L-1) ) { // ns is to the right of os
		  cDist[ns] = cDist[os] + 1;
		  cDist[ns+N] = cDist[os+N];
	  }
	  else if ( (delta == -L) || (delta == N-L) ){  // ns is above os
		  cDist[ns] = cDist[os];
		  cDist[ns+N] = cDist[os+N] + 1;
	  }
	  else if ( (delta == L) || (delta == L-N) ){  // ns is below os
		  cDist[ns] = cDist[os];
		  cDist[ns+N] = cDist[os+N] - 1;
	  }
	  else{
		 System.err.println("Internal Inconsistancy in Setting Position Relative to Neighbor"); 
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
  
  public int getClusterSizeV2(int s) {
	  
	  int counter = 0;
	  int pcns = pcn[s];
	  
	    for(int jj = 0 ; jj < N ; jj++){
	    	if(pcn[jj] == pcns) counter++;
	    }
	    return counter;
  }

  public int getLargestCluster() {
	  int max = 1;
	  for (int j = 1 ; j < N+1 ; j++){
		  if(numClusters[j]>0) max = j;
	  }
	  
	  return max;
  }
  
  public int getLargestCluster(int old) {
	  int max = old;
	  for (int j = old ; j < N+1 ; j++){
		  if(numClusters[j]>0) max = j;
	  }
	  
	  return max;
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
  
  public int getDist(int site){
	  
	  if(site < 2*N){
		  if(cDist[site]!=EMPTY){
			  return cDist[site];
		  }
		  else{
			  return N-1;
		  }
	  }
	  else{
		  return -100;
	  }
  }
  
  public int whichCluster(){
	  
	  return PercolationCluster;
  }
  
  
}