package rachele.damage2D;

import kip.util.Random;
import scikit.dataset.Accumulator;
import scikit.jobs.params.Parameters;

public class OFC_Lattice {

	public int L, R, Lp, dx, dt;
	public int time, cgTime, failSite;
	public double tStress, rStress, dissParam, maxResNoise, metric0;
	public Random random = new Random();
	public int [] cgCount;
	public double [] stress;
	public double [] stressTimeAve;
	boolean [] failList;
	public Accumulator iterAcc;
	public Accumulator inverseMetricAcc;
	
	public OFC_Lattice(Parameters params) {
		setParameters(params);
		initLattice();
	}
	
	public void setParameters(Parameters params) {
		Lp = params.iget("CG size");
		dx = params.iget("dx");
		L = Lp*dx;
		params.set("L", L);
		
		dt = params.iget("dt");
		
		random.setSeed(params.iget("Random Seed"));
		
		tStress = 1.0;
		dissParam = params.fget("Dissipation Param");
		rStress = params.fget("Residual Stress");
		maxResNoise = params.fget("Res. Max Noise");
	}
	
	public void initLattice(){
		stress = new double [L*L];
		stressTimeAve = new double [L*L];
		cgCount = new int [Lp*Lp];
		failList = new boolean [L*L];
		iterAcc = new Accumulator(dt);
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
		failSite = siteWithMaxStress();
		int site = findCG_site(failSite);
		cgCount[site] +=1;
	}
	
	public void step(){
		//Bring to failure
		double stressAdd = tStress - stress[failSite];
		if(stressAdd < 0) System.out.println("Error: stress already above failure");
		for (int i = 0; i < L*L; i++) stress[i] += stressAdd;
		failList[failSite] = true;
		//Fail the failSite
		boolean failAgain = checkFailList();
		int iteration = 0;
		while(failAgain){
			fail();
			iteration += 1;
			failAgain = checkFailList();
		}
		iterAcc.accum(time, iteration);
		for(int i=0; i < L*L; i++)
			stressTimeAve[i] = (stressTimeAve[i]*(double)time + stress[i])/(double)(time+1);
		time += 1;
		if(time%dt==0) cgTime += 1;
		calcMetric();
	}

	int findCG_site(int s){
		int x = s % L;
		int y = s / L;
		int xp = x / dx;
		int yp = y / dx;
		int cgSite = yp *Lp + xp;
		return cgSite;
	}
	
	void calcMetric(){
		double spaceSum = 0;
		for (int i = 0; i < L*L; i++) spaceSum += stressTimeAve[i];
		double spaceTimeStressAve = (spaceSum)/(double)(L*L);
		double 
		metricSum = 0;
		for (int i = 0; i < L*L; i++) metricSum += Math.pow(stressTimeAve[i] - spaceTimeStressAve, 2);
		double inverseMetric = (double)(L*L)*metric0/metricSum;
		inverseMetricAcc.accum(time, inverseMetric);
	}
	
	void calcMetric0(){
		double spaceSum = 0;
		for (int i = 0; i < L*L; i++) spaceSum += stress[i];
		double spaceStressAve = (spaceSum)/(double)(L*L);
		double metricSum = 0;
		for (int i = 0; i < L*L; i++) metricSum += Math.pow(stress[i] - spaceStressAve, 2);
		metric0 = metricSum / (double)(L*L);
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
		for(int i = 0; i < L*L; i++){
			if (stress[i] >= tStress){
				failList[i] = true;
				failCheck = true;
			}
		}
		return failCheck;
	}
	
	void fail(){
		//Bring to Residual
		int noNbors = L*L-1;
		//maybe need to generate a random list here??
		//Changed so that all stress is summed 1st before added
		// so I don't think a random list is needed, but
		// may be needed for finite ranges.
		double excessStressSum = 0.0;
		for (int site = 0; site < L*L; site++){
			if(failList[site]){
				double resNoise = maxResNoise*(random.nextFloat()*2.0-1.0);
				double excessStressPerNbor = ((1-dissParam)*(stress[site] - rStress) - resNoise)/(double)noNbors;
				excessStressSum =+ excessStressPerNbor;
				stress[site] = rStress - excessStressPerNbor;
				failList[site] = false;
			}
		}
		for(int i = 0; i < L*L; i++){
			stress[i] += excessStressSum;
		}
	}
	
	public int siteWithMaxStress(){
		double max = stress[0];
		int maxSite = 0;
		for (int i = 1; i < L*L; i++){
			if (stress[i] > max){
				max = stress[i];
				maxSite = i;
			}
		}
		return maxSite;
	}
	
}
