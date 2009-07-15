package rachele.damage2D;

import scikit.jobs.params.Parameters;
import scikit.util.DoubleArray;

public class OFC_Lattice extends AbstractCG_OFC{
 

	double metric0;					//metric measured at time t = 0;
	double lastCG_RecordTime;		//CG time parameter to keep track of last time CG info was recorded
	double lastCG_StressRecordTime; //Basically same as lasCG_RecordTime, but separate to keep calcInverseMeric subroutines independent in case one is not performed
	double lastRecordTime;
	public int lowerCutoff;
	
	public double [] stressTimeAve;
	public double [] CG_ActivityTimeAve;
	public double [] CG_StressTimeAve;
	public double [] CG_Stress;
	public int [] plateUpdateFailLocations;  // Record all fail sites for one plate update
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
		dt = 1.0/params.fget("Coarse Grained dt");
		
//		CG_dt = params.fget("Coarse Grained dt");
//		dt = params.fget("dt");
		
		random.setSeed(params.iget("Random Seed"));
		
		tStress = 1.0;
		dissParam = params.fget("Dissipation Param");   // = alpha
		rStress = params.fget("Residual Stress");
		maxResNoise = params.fget("Res. Max Noise");
		lastCG_RecordTime = 0;
		lastCG_StressRecordTime = 0;
		lastRecordTime = 0;
		
		lowerCutoff = params.iget("Lower Cutoff");
	}
	
	public void initLattice(){
		stress = new double [N];
		stressTimeAve = new double [N];
		epicenterCount = new int [Np];
		CG_ActivityTimeAve = new double [Np];
		CG_StressTimeAve = new double [Np];
		plateUpdateFailLocations = new int [N];
		CG_Stress = new double [Np];
		failList = new boolean [N];
		for (int i = 0; i < N; i++){
			stress[i] = random.nextFloat()*(tStress-rStress)+rStress;
			stressTimeAve[i] = stress[i];
			failList[i] = false;
		}
		for (int i = 0; i < Np; i++){
			epicenterCount[i] = 0;
			CG_ActivityTimeAve[i] = 0;
			CG_StressTimeAve[i] = 0;
		}
		metric0 = calcMetric0();
		cg_time = 0;
		plateUpdates = 0;

	}
	
	public void initEquilibrate(int maxPlateUpdates){
			plateUpdates = -maxPlateUpdates;
	}
	
	public void equilibrate(){
		//find the epicenter
		epicenterSite = siteWithMaxStress();
		//bring epicenter to failure by adding equal stress to all sites
		double stressAdd = tStress - stress[epicenterSite];
		for (int i = 0; i < N; i++) stress[i] += stressAdd;
		//find and fail all sites in avalanche
		if(fullyConnected){
			//failList records all sites that will fail in one iteration
			//1st iteration includes the epicenter only
			failList[epicenterSite] = true;
			//fail the epicenter
			fc_fail();			
			//is there another iteration of failures due to epicenter failure?  If so, make a new failList
			boolean failAgain = fc_checkFailList();
			// new iteration if there are more failures due to previous iterations
			while(failAgain){
				fc_fail();
				failAgain = fc_checkFailList();
			}
		}else{
			failSiteWithRange(epicenterSite);
			int nextSiteToFail = checkFailWithRange();
			while (nextSiteToFail >= 0){
				failSiteWithRange(nextSiteToFail);
				nextSiteToFail = checkFailWithRange();
			}
		}
		plateUpdates += 1;

	}
	
//	public void prestep(){
//		epicenterSite = siteWithMaxStress();
//		int site = findCG_site(epicenterSite);
//		epicenterCount[site] +=1;
//	}
	
	public void clearPlateUpdateFileLocs(){
		for(int i = 0; i < N; i++) plateUpdateFailLocations[i] = 0;
	}
	
	public void clearFailList(){
		for(int i = 0; i < N; i++) failList[i] = false;
	}
	
	public void step(){
		clearPlateUpdateFileLocs();
		clearFailList();
		//initiate avalanche size
		avSize = 0;
		//find the epicenter
		epicenterSite = siteWithMaxStress();
		//bring epicenter to failure by adding equal stress to all sites
		double stressAdd = tStress - stress[epicenterSite];
		if(stressAdd < 0) System.out.println("Error: stress already above failure");
		for (int i = 0; i < N; i++) stress[i] += stressAdd;
		//find and fail all sites in avalanche
		if(fullyConnected){
			//failList records all sites that will fail in one iteration
			//1st iteration includes the epicenter only
			failList[epicenterSite] = true;
			//fail the epicenter
			fc_failAndCount();			
			//is there another iteration of failures due to epicenter failure?  If so, make a new failList
			boolean failAgain = fc_checkFailList();
			// new iteration if there are more failures due to previous iterations
			while(failAgain){
				fc_failAndCount();
				failAgain = fc_checkFailList();


//				failList[epicenterSite] = true;
//				//			fc_fail(epicenter);
//				//Fail the failSite
//				boolean failAgain = fc_checkFailList();
//				while(failAgain){
//				avSize += fc_fail();
//				failAgain = fc_checkFailList();
			}
		}else{
			plateUpdateFailLocations[epicenterSite] += 1;			
			avSize = 1;
			failSiteWithRange(epicenterSite);
			int nextSiteToFail = checkFailWithRange();
			if(nextSiteToFail >= 0){
				while (nextSiteToFail >= 0){
					failSiteWithRange(nextSiteToFail);
//					System.out.println("sitetofail = " + nextSiteToFail);
					plateUpdateFailLocations[nextSiteToFail] += 1;
					nextSiteToFail = checkFailWithRange();
					avSize += 1;
				}
			}
		}
		if (avSize >= lowerCutoff) 	epicenterCount[findCG_site(epicenterSite)] +=1;
		cg_time += dt;
		plateUpdates += 1;
//		System.out.println("av size = " + avSize + " at " + plateUpdates);
	}

	int checkFailWithRange(){
		int s = siteWithMaxStress();
		if (stress[s] < tStress) s = -1; 
		//		System.out.println("Stress of max Site = " + stress[s] + " at site " + s);
		return s;
	}

	boolean fc_checkFailList(){
		boolean failCheck = false;
		for(int i = 0; i < N; i++){
			if (stress[i] >= tStress){
				failList[i] = true;
				failCheck = true;
			}else{
				failList[i] = false;
			}
		}
		return failCheck;
	}
	
	void fc_fail(){
		double excessStressSum = 0.0;
		for (int site = 0; site < N; site++){
			if (failList[site]){
				double rStressWithNoise = calcResNoise(site);
				double excessStressPerNbor = calcStressPerNbor(site, rStressWithNoise);
				excessStressSum += excessStressPerNbor;
				stress[site] = rStressWithNoise - excessStressPerNbor;
			}
		}
		for(int site = 0; site < N; site++)
			stress[site] += excessStressSum;		
	}
	
	void fc_failAndCount(){
		double excessStressSum = 0.0;
		for (int site = 0; site < N; site++){
			if (failList[site]){
//				System.out.println("fail site " + site);
				avSize += 1;
				plateUpdateFailLocations[site] += 1;
				double rStressWithNoise = calcResNoise(site);
				double excessStressPerNbor = calcStressPerNbor(site, rStressWithNoise);
				excessStressSum += excessStressPerNbor;
				stress[site] = rStressWithNoise - excessStressPerNbor;
			}
		}
		for(int site = 0; site < N; site++)
			stress[site] += excessStressSum;		
	}
	
//	int fc_fail(){
//		//Bring to Residual
//		//maybe need to generate a random list here??
//		//Changed so that all stress is summed 1st before added
//		// so I don't think a random list is needed, but
//		// may be needed for finite ranges.
//		double excessStressSum = 0.0;
//		int noFailed = 0;
//		for (int site = 0; site < N; site++){
//			if(failList[site]){
//				double rStressWithNoise = calcResNoise(site);
//				double excessStressPerNbor = calcStressPerNbor(site, rStressWithNoise);
////				double resNoise = maxResNoise*(random.nextFloat()*2.0-1.0);
////				double excessStressPerNbor = ((1-dissParam)*(stress[site] - rStress) - resNoise)/(double)noNbors;
//				excessStressSum += excessStressPerNbor;
//				stress[site] = rStressWithNoise - excessStressPerNbor;
//				failList[site] = false;
//				noFailed += 1;
//				plateUpdateFailLocations[site] += 1;
//			}
//		}
//		for(int i = 0; i < N; i++){
//			stress[i] += excessStressSum;
//		}
//		return noFailed;
//	}
	
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
		double del_t = cg_time -lastRecordTime;
		for(int i=0; i < N; i++)
			stressTimeAve[i] = (stressTimeAve[i]*(lastRecordTime)+ stress[i]*del_t)/(cg_time);
		double spaceSum = DoubleArray.sum(stressTimeAve);
		double spaceTimeStressAve = (spaceSum)/(double)(N);
		double metricSum = 0;
		for (int i = 0; i < N; i++) metricSum += Math.pow(stressTimeAve[i] - spaceTimeStressAve, 2);
		double inverseMetric = (double)(N)*metric0/metricSum;
		lastRecordTime = cg_time;
		return inverseMetric;
	}

	public double calcStressMetric(){
		double del_t = cg_time -lastRecordTime;
		for(int i=0; i < N; i++)
			stressTimeAve[i] = (stressTimeAve[i]*(lastRecordTime)+ stress[i]*del_t)/(cg_time);
		double spaceSum = DoubleArray.sum(stressTimeAve);
		double spaceTimeStressAve = (spaceSum)/(double)(N);
		double metricSum = 0;
		for (int i = 0; i < N; i++) metricSum += Math.pow(stressTimeAve[i] - spaceTimeStressAve, 2);
		double stressMetric = metricSum/(double)N;
		lastRecordTime = cg_time;
		return stressMetric;
	}

	//fixed to CG in space and time
	public double calcCG_stressMetric(){
		double del_t = cg_time - lastCG_StressRecordTime;
		//		double [] cgStressConfig = findCG_StressConfig();
		CG_Stress = findCG_StressConfig();
		// accumulate the most recent time step into the time ave
		for (int i = 0; i < Np; i++)
			CG_StressTimeAve[i] = (lastCG_StressRecordTime*CG_StressTimeAve[i] + del_t*CG_Stress[i]) / (cg_time);
		double CG_SpaceAve = DoubleArray.sum(CG_StressTimeAve)/(double)Np;
		double CG_metric = 0;
		//take the space ave
		for (int i = 0; i < Np; i++){
			CG_metric += Math.pow(CG_StressTimeAve[i] - CG_SpaceAve,2);
		}
		CG_metric /= (double)Np;
		lastCG_StressRecordTime = cg_time;
		return CG_metric; 
	}
	
//	//only CGs in time, not space
//	public double calcCG_stressMetricCGTime(){
//		double del_t = cg_time - lastCG_StressRecordTime;
//
//		// accumulate the most recent time step into the time ave
//		for (int i = 0; i < N; i++)
//			CG_StressTimeAve[i] = (lastCG_StressRecordTime*CG_StressTimeAve[i] + del_t*stress[i]) / (cg_time);
//		double CG_SpaceAve = DoubleArray.sum(CG_StressTimeAve)/(double)N;
//		double CG_metric = 0;
//		//take the space ave
//		for (int i = 0; i < N; i++){
//			CG_metric += Math.pow(CG_StressTimeAve[i] - CG_SpaceAve,2);
//		}
//		CG_metric /= (double)N;
//		lastCG_StressRecordTime = cg_time;
//		return CG_metric; 
//	}
	
	public double [] findCG_StressConfig(){
		double [] cgSpaceStressConfig = new double [Np];
		for(int i = 0; i < N; i++){
			int cgSite = findCG_site(i);
			cgSpaceStressConfig[cgSite] += stress[i];
		}
		int cgBoxSize = N/Np;
		for(int cgSite = 0; cgSite < Np; cgSite++){
			cgSpaceStressConfig[cgSite] /= (double) cgBoxSize;
		}
		return cgSpaceStressConfig;
	}
	
	public double calcCG_activityMetric(){
		double del_t = cg_time - lastCG_RecordTime;
//		System.out.println("del_t = " + del_t + " cg_time = " + cg_time + " lasCG_recordTime = " + lastCG_RecordTime);
		for (int i = 0; i < Np; i++){
//			double countWithCutoff;
//			if (epicenterCount[i] < lowerCutoff) countWithCutoff = 0;
//			else countWithCutoff = epicenterCount[i];
//			CG_ActivityTimeAve[i] = (lastCG_RecordTime*CG_ActivityTimeAve[i] + del_t*countWithCutoff) / (time);
			CG_ActivityTimeAve[i] = (lastCG_RecordTime*CG_ActivityTimeAve[i] + del_t*epicenterCount[i]) / (cg_time);
		}
		double CG_SpaceAve = DoubleArray.sum(CG_ActivityTimeAve)/(double)Np;
		double CG_metric = 0;
		for (int i = 0; i < Np; i++){
			CG_metric += Math.pow(CG_ActivityTimeAve[i] - CG_SpaceAve,2);
		}
		CG_metric /= (double)Np;
		lastCG_RecordTime = cg_time;
		for (int i = 0; i < Np; i++) epicenterCount[i] = 0;
		return CG_metric; 
	}
	



	
}
