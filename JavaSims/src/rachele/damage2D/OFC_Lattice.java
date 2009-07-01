package rachele.damage2D;

import scikit.jobs.params.Parameters;
import scikit.util.DoubleArray;

public class OFC_Lattice extends AbstractCG_OFC{


	double metric0;					//metric measured at time t = 0;
	double lastCG_RecordTime;		//CG time parameter to keep track of last time CG info was recorded
	double lastRecordTime;
	public int lowerCutoff;
	
	public double [] stressTimeAve;
	public double [] CG_TimeAve;
	boolean [] failList;
	
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
		
//		CG_dt = params.fget("Coarse Grained dt");
//		dt = params.fget("dt");
		
		random.setSeed(params.iget("Random Seed"));
		
		tStress = 1.0;
		dissParam = params.fget("Dissipation Param");   // = alpha
		rStress = params.fget("Residual Stress");
		maxResNoise = params.fget("Res. Max Noise");
		lastCG_RecordTime = 0;
		lastRecordTime = 0;
		
		lowerCutoff = params.iget("Lower Cutoff");
	}
	
	public void initLattice(){
		stress = new double [N];
		stressTimeAve = new double [N];
		epicenterCount = new int [Np];
		CG_TimeAve = new double [Np];
		failList = new boolean [N];
		for (int i = 0; i < N; i++){
			stress[i] = random.nextFloat()*(tStress-rStress)+rStress;
			stressTimeAve[i] = stress[i];
			failList[i] = false;
		}
		for (int i = 0; i < Np; i++){
			epicenterCount[i] = 0;
			CG_TimeAve[i] = 0;
		}
		metric0 = calcMetric0();
		time = 0;
		plateUpdates = 0;
	}
	
	public void initEquilibrate(int maxPlateUpdates){
			plateUpdates = -maxPlateUpdates;
	}
	
	public void equilibrate(){

		epicenterSite = siteWithMaxStress();
		double stressAdd = tStress - stress[epicenterSite];
		for (int i = 0; i < N; i++) stress[i] += stressAdd;
		if(fullyConnected){
			failList[epicenterSite] = true;
			boolean failAgain = checkFailList();
			while(failAgain){
				avSize += fail();
				failAgain = checkFailList();
			}
		}else{
			failSiteWithRange(epicenterSite);
			int nextSiteToFail = checkFail();
			while (nextSiteToFail >= 0){
				failSiteWithRange(nextSiteToFail);
				nextSiteToFail = checkFail();
			}
		}
		plateUpdates += 1;

	}
	
	public void prestep(){
		epicenterSite = siteWithMaxStress();
		int site = findCG_site(epicenterSite);
		epicenterCount[site] +=1;
	}
	
	public void step(){
		//Bring to failure
		double stressAdd = tStress - stress[epicenterSite];
//		dt = stressAdd;
		dt=1.0;
		if(stressAdd < 0) System.out.println("Error: stress already above failure");
		for (int i = 0; i < N; i++) stress[i] += stressAdd;
		avSize = 1;
	
		if(fullyConnected){
			failList[epicenterSite] = true;
			//Fail the failSite
			boolean failAgain = checkFailList();
			while(failAgain){
				avSize += fail();
				failAgain = checkFailList();
			}
		}else{
			failSiteWithRange(epicenterSite);
			int nextSiteToFail = checkFail();
			while (nextSiteToFail >= 0){
				failSiteWithRange(nextSiteToFail);
				nextSiteToFail = checkFail();
				avSize += 1;
			}
		}

		time += dt;
		plateUpdates += 1;
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
	
	public double calcInverseMetric(){
		double del_t = time -lastRecordTime;
		for(int i=0; i < N; i++)
			stressTimeAve[i] = (stressTimeAve[i]*(lastRecordTime)+ stress[i]*del_t)/(time);
		double spaceSum = DoubleArray.sum(stressTimeAve);
		double spaceTimeStressAve = (spaceSum)/(double)(N);
		double metricSum = 0;
		for (int i = 0; i < N; i++) metricSum += Math.pow(stressTimeAve[i] - spaceTimeStressAve, 2);
		double inverseMetric = (double)(N)*metric0/metricSum;
		lastRecordTime = time;
		return inverseMetric;
	}
	

	
	public double calcCG_activity_Metric(){
		double del_t = time - lastCG_RecordTime;
		for (int i = 0; i < Np; i++){
			double countWithCutoff;
			if (epicenterCount[i] < lowerCutoff) countWithCutoff = 0;
			else countWithCutoff = epicenterCount[i];
			CG_TimeAve[i] = (lastCG_RecordTime*CG_TimeAve[i] + del_t*countWithCutoff) / (time);
//			CG_TimeAve[i] = (lastCG_RecordTime*CG_TimeAve[i] + del_t*cgCount[i]) / (time);
		}
		double CG_SpaceAve = DoubleArray.sum(CG_TimeAve)/(double)Np;
		double CG_metric = 0;
		for (int i = 0; i < Lp; i++){
			CG_metric += Math.pow(CG_TimeAve[i] - CG_SpaceAve,2);
		}
		CG_metric /= (double)Np;
		lastCG_RecordTime = time;
		for (int i = 0; i < Np; i++) epicenterCount[i] = 0;
		return CG_metric; 
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


	
}
