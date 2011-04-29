package rachele.damage2D.multidx;

import kip.util.Random;

public class AbstractOFC_Multidx {

	public int pow;				// System length = 2**pow
	public int L;				// Linear lattice size
	public int N;				// Total no of lattice sites = LxL
	public int R;				// interaction range- stress this distance, R = 0 -> fully connected
	public int noNbors;			// no of neighbors within given interaction range
	public int plateUpdates;	// no of times avalanch has been forced by increassing stress to all sites
	public double cg_time;			// coarse grained time 
	public double dt;			// 1.0/coarse grained time chunk. (Coarse grained units of time)
	public int epicenterSite; 	// beginning site of current avalanch
	public int avSize; 			// avalanch size
	double tStress;				// threshold stress
	double rStress;				// base residual stress		
	double maxResNoise; 		// actual residual stress = rStress + residual noise with maxResNoise = max residual noise
	public boolean fullyConnected;
	public int [] epicenterCount;
	public double [] stress;
	double [] alpha;			// dissipation paramater, aka alpha	
	protected Random random = new Random();
	
	public int findCG_site(int s, int blockSize){
		int x = s % L;
		int y = s / L;
		int xp = x / blockSize;
		int yp = y / blockSize;
		int cgSite = yp *(L/blockSize) + xp;
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
		return resNoise;
	}
	
	double calcStressPerNbor(int s, double rStressWithNoise){
		double stressPerNbor = ((1-alpha[s])*(stress[s] - rStressWithNoise))/(double)noNbors;
		return stressPerNbor;
	}

	
	public int findCgArrayIndex(int cgSite, int ip){
		int ret = 0;
		for(int i = 0; i < ip; i++){
			int dx = (int)Math.pow(2, i);
			int Np = N/((int) Math.pow(dx, 2));
			ret += Np;
		}
		ret += cgSite;
		return ret;
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
	
	void write(String text){
		System.out.println(text);
	}
	
	void write(String text, boolean bool){
		System.out.println(text + " " + bool);
	}

}

