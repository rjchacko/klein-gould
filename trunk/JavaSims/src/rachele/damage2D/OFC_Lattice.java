package rachele.damage2D;

import kip.util.Random;
import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.jobs.params.Parameters;

public class OFC_Lattice extends AbstractCG_OFC{

	public int time, cgTime, epicenterSite;
	public double tStress, rStress, dissParam, maxResNoise, metric0;
	public Random random = new Random();
	public int [] cgCount;
	public double [] stress;
	public double [] stressTimeAve;
	boolean [] failList;
	public Accumulator iterAcc;
	public Histogram sizeHist;
	public Accumulator inverseMetricAcc;
	
	public OFC_Lattice(Parameters params) {
		setParameters(params);
		initLattice();
	}
	
	public void setParameters(Parameters params) {
		Lp = params.iget("CG size");
		dx = params.iget("dx");
		L = Lp*dx;
		N = L*L;
		Np = Lp*Lp;
		params.set("L", L);
		R = params.iget("R");
		if (R==0){
			fullyConnected = true;
			noNbors = N-1;
		}else{
			fullyConnected = false;
			noNbors = findCircleNbors(R);
		}
		System.out.println("Fully Connected = " + fullyConnected);
		
		dt = params.iget("dt");
		
		random.setSeed(params.iget("Random Seed"));
		
		tStress = 1.0;
		dissParam = params.fget("Dissipation Param");
		rStress = params.fget("Residual Stress");
		maxResNoise = params.fget("Res. Max Noise");
	}
	
	public void initLattice(){
		stress = new double [N];
		stressTimeAve = new double [N];
		cgCount = new int [Np];
		failList = new boolean [N];
		iterAcc = new Accumulator(dt);
		sizeHist = new Histogram(1);
		inverseMetricAcc = new Accumulator(dt);
		for (int i = 0; i < L*L; i++){
			stress[i] = random.nextFloat();
			stressTimeAve[i] = stress[i];
			failList[i] = false;
//			System.out.println("site " + i + " = " + stress[i]);
		}
		calcMetric0();
		time = 1;
		cgTime = 1;
	}
	
	public void prestep(){
		epicenterSite = siteWithMaxStress();
		int site = findCG_site(epicenterSite);
		cgCount[site] +=1;
	}
	
	public void step(){
		//Bring to failure
		double stressAdd = tStress - stress[epicenterSite];
		if(stressAdd < 0) System.out.println("Error: stress already above failure");
		for (int i = 0; i < N; i++) stress[i] += stressAdd;
		int iteration = 0;
		int size = 1;
	
		if(fullyConnected){
			failList[epicenterSite] = true;
			//Fail the failSite
			boolean failAgain = checkFailList();
			while(failAgain){
				size += fail();
				iteration += 1;
				failAgain = checkFailList();
			}
		}else{
			failSiteWithRange(epicenterSite);
			int nextSiteToFail = checkFail();
			while (nextSiteToFail >= 0){
				failSiteWithRange(nextSiteToFail);
				nextSiteToFail = checkFail();
				size += 1;
				iteration += 1;
//				System.out.println("Next Site to Fail "+ nextSiteToFail);
			}
		}
		iterAcc.accum(time, iteration);
		sizeHist.accum(size);
		for(int i=0; i < N; i++)
			stressTimeAve[i] = (stressTimeAve[i]*(double)time + stress[i])/(double)(time+1);
		time += 1;
		if(time%dt==0) cgTime += 1;
		calcMetric();
	}

	int checkFail(){
		int s = siteWithMaxStress();
		if (stress[s] < tStress) s = -1; 
//		System.out.println("Stress of max Site = " + stress[s] + " at site " + s);
		return s;
	}
	
	void failSiteWithRange(int s){
		double resStressWithNoise = calcResNoise(s);
		double stressPerNbor = calcStressPerNbor(s, resStressWithNoise);
		stress[s] = resStressWithNoise;
		int x = s%L;
		int y = s/L;
		 for(int dy = -R; dy <= R; dy++){
			 for(int dx = -R; dx <= R; dx++){
				 double distance = Math.sqrt(dx*dx + dy*dy);
				 if (distance <= R){
					 int xx = (x+dx+L)%L;
					 int yy = (y+dy+L)%L;
					 int nborSite = yy*L+xx;
					 stress[nborSite] += stressPerNbor;
				 }
			 }
		 }
		
	}
	
	void calcMetric(){
		double spaceSum = 0;
		for (int i = 0; i < N; i++) spaceSum += stressTimeAve[i];
		double spaceTimeStressAve = (spaceSum)/(double)(L*L);
		double 
		metricSum = 0;
		for (int i = 0; i < N; i++) metricSum += Math.pow(stressTimeAve[i] - spaceTimeStressAve, 2);
		double inverseMetric = (double)(N)*metric0/metricSum;
		inverseMetricAcc.accum(time, inverseMetric);
	}
	
	void calcMetric0(){
		double spaceSum = 0;
		for (int i = 0; i < N; i++) spaceSum += stress[i];
		double spaceStressAve = (spaceSum)/(double)(N);
		double metricSum = 0;
		for (int i = 0; i < N; i++) metricSum += Math.pow(stress[i] - spaceStressAve, 2);
		metric0 = metricSum / (double)N;
	}
	
	void calcCG_Metric(){
//		double spaceSum = 0;
//		for (int i = 0; i < L*L; i++) spaceSum += stressTimeAve[i];
//		double spaceTimeStressAve = (spaceSum)/(double)(L*L);
//		double 
//		metricSum = 0;
//		for (int i = 0; i < L*L; i++) metricSum += Math.pow(stressTimeAve[i] - spaceTimeStressAve, 2);
//		double inverseMetric = (double)(L*L)*metric0/metricSum;
//		inverseMetricAcc.accum(time, inverseMetric);
	}
	
	boolean checkFailList(){
		boolean failCheck = false;
		for(int i = 0; i < N; i++){
			if (stress[i] >= tStress){
				failList[i] = true;
				failCheck = true;
			}
		}
		return failCheck;
	}
	
	int fail(){
		//Bring to Residual
		//maybe need to generate a random list here??
		//Changed so that all stress is summed 1st before added
		// so I don't think a random list is needed, but
		// may be needed for finite ranges.
		double excessStressSum = 0.0;
		int noFailed = 0;
		for (int site = 0; site < N; site++){
			if(failList[site]){
				double rStressWithNoise = calcResNoise(site);
				double excessStressPerNbor = calcStressPerNbor(site, rStressWithNoise);
//				double resNoise = maxResNoise*(random.nextFloat()*2.0-1.0);
//				double excessStressPerNbor = ((1-dissParam)*(stress[site] - rStress) - resNoise)/(double)noNbors;
				excessStressSum =+ excessStressPerNbor;
				stress[site] = rStressWithNoise - excessStressPerNbor;
				failList[site] = false;
				noFailed += 1;
			}
		}
		for(int i = 0; i < N; i++){
			stress[i] += excessStressSum;
		}
		return noFailed;
	}
	
	double calcResNoise(int s){
		double resNoise = maxResNoise*(random.nextFloat()*2.0-1.0)+rStress;
		return resNoise;
	}
	
	double calcStressPerNbor(int s, double rStressWithNoise){
		double stressPerNbor = ((1-dissParam)*(stress[s] - rStressWithNoise))/(double)noNbors;
		return stressPerNbor;
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
	
}
