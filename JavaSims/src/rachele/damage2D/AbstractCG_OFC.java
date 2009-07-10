package rachele.damage2D;

import kip.util.Random;

public class AbstractCG_OFC {

	public int L;				// Linear lattice size
	public int N;				// Total no of lattice sites = LxL
	public int R;				// interaction range- stress this distance, R = 0 -> fully connected
	public int Lp;				// no of coarse grained blocks (linear)
	public int Np;				// total no of coarse grained blocks = LpxLp
	public int dx;				// coarse grained size = L/Lp 
	public int noNbors;			// no of neighbors within given interaction range
	public int plateUpdates;	// no of times avalanch has been forced by increassing stress to all sites
	public int time;			// time (interpreted either as no of updates or proportional to size of stress update ??)
	public int dt;			// Either 1 or proportional to stress update (see above)
	public int epicenterSite; 	// beginning site of current avalanch
	public int avSize; 			// avalanch size
	double tStress;				// threshold stress
	double rStress;				// base residual stress
	double dissParam;			// dissipation paramater, aka alpha			
	double maxResNoise; 		// actual residual stress = rStress + residual noise with maxResNoise = max residual noise
	boolean fullyConnected;
	public int [] epicenterCount;
	public double [] stress;
	Random random = new Random();
	
	int findCG_site(int s){
		int x = s % L;
		int y = s / L;
		int xp = x / dx;
		int yp = y / dx;
		int cgSite = yp *Lp + xp;
		return cgSite;
	}

	int findCircleNbors(int r){
		int no = 0;
		for(int y = -r; y <= r; y++){
			for(int x = -r; x <= r; x++){
				double distance = Math.sqrt(x*x + y*y);
				if (distance <= r) no += 1;
			}
		}
		no -= 1;  //do not count the site itself 
		return no;
	}

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
	
	double calcResNoise(int s){
		double resNoise = maxResNoise*(random.nextFloat()*2.0-1.0) + rStress;
//		System.out.println("Res noise = " + resNoise);
		return resNoise;
	}
	
	double calcStressPerNbor(int s, double rStressWithNoise){
		double stressPerNbor = ((1-dissParam)*(stress[s] - rStressWithNoise))/(double)noNbors;
		return stressPerNbor;
	}

	
	double  calcMetric0(){
		double spaceSum = 0;
		for (int i = 0; i < N; i++) spaceSum += stress[i];
		double spaceStressAve = (spaceSum)/(double)(N);
		double metricSum = 0;
		for (int i = 0; i < N; i++) metricSum += Math.pow(stress[i] - spaceStressAve, 2);
		double metric0 = metricSum / (double)N;
		return metric0;
	}

}
