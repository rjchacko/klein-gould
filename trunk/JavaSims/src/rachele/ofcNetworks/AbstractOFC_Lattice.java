/**
 * AbstractOCF_Lattice is a long stress transfer range OFC app
 * 
 * @author Rachele Dominguez
 * @version 2.0 August 3, 2012
 * 
 * Changed update procedure to update all sites simultaneously for a given lattice scan.  
 * Now the failure procedure is cut into two steps: fail and distribute stress.
 * First all sites above the threshold are failed.
 * Then all stress is distributed.
 * Then all sites are scanned again to see if they are above threshold.
 * Requires a new array stressAdd to keep track of what stress should be added at the end of the failures.
 */

package rachele.ofcNetworks;

import kip.util.Random;

public class AbstractOFC_Lattice {

	public int L;				// Linear lattice size
	public int N;				// Total no of lattice sites = LxL
	public double [] stress;
	public double [] stressAdd;	//Keeps track of what stress is to be added to what sites during an update cycle.
	public double alpha;			// dissipation paramater, aka alpha	
	public int plateUpdates;
	public int [][] nbor;
	public int epicenterSite; 	// beginning site of current avalanch
//	public int avSize; 			// avalanch size
	double tStress = 1.0;				// threshold stress
	double rStress = 0;				// base residual stress		
	double maxResNoise = 0.0; 		// actual residual stress = rStress + residual noise with maxResNoise = max residual noise
	int [] sitesToFail;
	public int noNbors;
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
		initialScanLattice();
		while(sitesToFail[0] >= 0){
			clearStressAdd();
			int i = 0;
			while(sitesToFail[i] >= 0){
				fail(sitesToFail[i]);
				i ++;
			}
			distributeStress();
			scanLattice();
			//System.out.println("Scan " + it + " has " + scanFails + " fails");
			//it ++;
			//double totalStress = 0;
			//for (int j = 0; j < N; j++) totalStress += stress[j];
			//System.out.println("Total stress = " + totalStress);
		}
		plateUpdates += 1;
		//System.out.println("PU = " + plateUpdates);
	}
	
	/**
	 * Fails a site:
	 * Calculates excess stress and collects it in the stressAdd array
	 */
	void fail(int failingSite){
		failCount ++;
		double stressPerNbor = (stress[failingSite] - rStress)*(1.0 - alpha)/(double)noNbors;
		for(int i = 0; i < noNbors; i++){
			int neighbor = nbor[failingSite][i];
			if(neighbor >= 0) stressAdd[neighbor] += stressPerNbor; 
		}
		stress[failingSite] = rStress;
	}
	
	/**
	 * Clears the array stressAdd.  Should be called at the begining of each 
	 * fail cycle.
	 */
	void clearStressAdd(){
		for(int i = 0; i < N; i++){
			stressAdd[i] = 0.0;
		}
	}
	
	
	/**
	 * Distributes all stress from failures of previous cycle. Should be called after 
	 * each complete cycle of failures.
	 */
	void distributeStress(){
		for(int i = 0; i < N; i++){
			stress[i] +=stressAdd[i];
		}
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
		stress[epicenterSite] -= stressAdd;			//????
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
		//System.out.println("sites to fail index = " + index);
	}

	/**
	 * Scans the lattice to find one site with max stress, then goes back to 
	 * collect all sites with that same max stress.  This step is to 
	 * take special care with the no noise condition, since many sites
	 * may be reset to exactly the same stress during a large event.
	 */
	void initialScanLattice(){
		for (int i = 0; i <= N; i++){				//Clear sitesToFail array
			sitesToFail[i]=-1;
		}
		int maxSite = siteWithMaxStress();			//Find max stress
		double maxStress = stress[maxSite];
	//	sitesToFail[0] = maxSite;
		int index = 0;
		for (int i = 0; i < N; i++){
			if(stress[i] >= maxStress){				//Identify all sites with max stress
				sitesToFail[index] = i;				//Store these sites to be failed
				index ++;
			}
		}
		double stressBump = tStress - maxStress;
		if(stressBump < 0) System.out.println("Error: stress already above failure");
		for (int i = 0; i < N; i++){				//Bring lattice to failure
			stress[i] += stressBump;
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
	
	/**
	 * Boundary conditions for now are OPEN only!
	 */
	void setNearNbors(){
		for(int i = 0; i < N; i++){
			int x = i % L;
			int y = i / L;
			
			// left neighbor
			if (x == 0){
				nbor[i][0] = -1;
			}else{
				nbor[i][0] = i - 1;
			}
			// bottom neighbor
			if (y == 0){
				nbor[i][1] = -1;
			}else{
				nbor[i][1] = i - L;
			}
			// right neighbor
			if (x == L-1){
				nbor[i][2] = -1;
			}else{
				nbor[i][2] = i + 1;
			}
			// top neighbor
			if (y == L-1){
				nbor[i][3] = -1;
			}else{
				nbor[i][3] = i + L;
			}
		}
	}
	
}
