package rachele.ofcNetworks;

import kip.util.Random;

public class AbstractOFC_Lattice {

	public int L;				// Linear lattice size
	public int N;				// Total no of lattice sites = LxL
	public double [] stress;
	public double alpha;			// dissipation paramater, aka alpha	
	public int plateUpdates;
	public int [][] nbor;
	public int epicenterSite; 	// beginning site of current avalanch
//	public int avSize; 			// avalanch size
	double tStress = 1.0;				// threshold stress
	double rStress = 0;				// base residual stress		
	double maxResNoise = 0.0; 		// actual residual stress = rStress + residual noise with maxResNoise = max residual noise
	int [] sitesToFail;
	int noNbors;
	public int failCount;
	
	protected Random random = new Random();
	

	
	/**
	 * Performs one OFC step:
	 * Brings lattice to failure.
	 * Fails site with max stress.
	 * Scans lattice for new sites above threshold.
	 * Fails these sites.
	 * Repeats scan until no sites above failure.
	 * Advances plateUpdates by one step.
	 */
	public void ofcStep(){
		failCount = 0;
		bringToFailure();
		fail(epicenterSite);
		scanLattice();
		//int it = 1;
		while(sitesToFail[0] >= 0){
			int i = 0;
			while(sitesToFail[i] >= 0){
				fail(sitesToFail[i]);
				i ++;
			}
			scanLattice();
			//System.out.println("Scan " + it + " has " + scanFails + " fails");
			//it ++;
		}
		plateUpdates += 1;
	}
	
	/**
	 * Fails a site:
	 * Calculates excess stress and distributes stress among
	 * neighbors
	 */
	void fail(int failingSite){
		failCount ++;
		double stressPerNbor = (stress[failingSite] - rStress)*(1.0 - alpha)/(double)noNbors;
		for(int i = 0; i < noNbors; i++){
			int neighbor = nbor[failingSite][i];
			if(neighbor >= 0) stress[neighbor] += stressPerNbor; 
		}
		stress[failingSite] = rStress;
	}
	
	/**
	 * Finds site with highest stress and brings lattice to failure
	 */
	void bringToFailure(){
		// find the epicenter
		epicenterSite = siteWithMaxStress();	
		// bring lattice to failure
		double stressAdd = tStress - stress[epicenterSite];
		if(stressAdd < 0) System.out.println("Error: stress already above failure");
		for (int i = 0; i < N; i++){
			stress[i] += stressAdd;
		}
		stress[epicenterSite] -= stressAdd;
	}
	
	/**
	 * Returns the site with the maximum stress
	 */
	public int siteWithMaxStress(){
		double max = stress[0];
		int maxSite = 0;
		for (int i = 1; i < N; i++){
			if (stress[i] > max){
				max = stress[i];
				maxSite = i;
			}
		}
		return maxSite;
	}

	/**
	 * Scans the lattice and stores all sites that are above
	 * threshold stress into array sitesToFail.
	 */
	void scanLattice(){
		for (int i = 0; i < N; i++){
			sitesToFail[i]=-1;
		}
		int index = 0;
		for (int i = 0; i < N; i++){
			if(stress[i] >= tStress){
				sitesToFail[index] = i;
				index ++;
			}
		}
	}
	
	void initStress(){
		for(int i = 0; i < N; i++){
			double rand = random.nextDouble();
			stress[i] = rand;
		}
	}
	
	void initSitesToFail(){
		for(int i = 0; i < N; i++){
			sitesToFail[i] = -1;
		}
	}
	
}
