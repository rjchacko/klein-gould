package rachele.damage2D;

import scikit.jobs.params.Parameters;
import scikit.util.DoubleArray;

/**
* OFC_Lattice is used for OFC only runs.  Probably will be replaced by OFC_DamageLattice-- 
* to get to OFC only, choose maxFail = 0 
*/
public class OFC_Lattice extends AbstractCG_OFC{
 

	double metric0;					//metric measured at time t = 0;
	double lastCG_RecordTime, lastCG_RecordTimePU;		//CG time parameter to keep track of last time CG info was recorded
	double lastCG_SizeActRecordTime, lastCG_SizeActRecordTimePU;
	double lastCG_StressRecordTime, lastCG_StressRecordTimePU; //Basically same as lasCG_RecordTime, but separate to keep calcInverseMeric subroutines independent in case one is not performed
	double lastRecordTime, lastRecordTimePU;
	public int lowerCutoff;
	
	public double [] stressTimeAve;
	public double [] CG_ActivityTimeAve;
	public double [] CG_SizeActTimeAve;
	public double [] CG_ActivityTimeAveDev;
	public double [] CG_SizeActTimeAveDev;
	public double [] CG_StressTimeAve;
	public double [] CG_Stress;
	public int [] plateUpdateFailLocations;  // Record all fail sites for one plate update
	public int [] cgFailCount; //Counts 1 in a cg block if there is any site that fails in that block.
	public int [] epicenterSize;
	boolean [] failList;
	public double am_sq_ave, am_ave_sq, sam_sq_ave, sam_ave_sq;
	
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
			System.out.println("No of neighbors = " + noNbors);
		}
		System.out.println("Fully Connected = " + fullyConnected);
		dt = 1.0/params.fget("Coarse Grained dt (PU)");
		
//		CG_dt = params.fget("Coarse Grained dt");
//		dt = params.fget("dt");
		
		random.setSeed(params.iget("Random Seed"));
		
		tStress = 1.0;
		dissParam = params.fget("Dissipation Param");   // = alpha
		rStress = params.fget("Residual Stress");
		maxResNoise = params.fget("Res. Max Noise");
		lastCG_RecordTime = 0;
		lastCG_SizeActRecordTime = 0;
		lastCG_StressRecordTime = 0;
		lastRecordTime = 0;
		
		lowerCutoff = params.iget("Lower Cutoff");
	}
	
	public void initLattice(){
		stress = new double [N];
		stressTimeAve = new double [N];
		epicenterCount = new int [Np];
		epicenterSize = new int [Np];
		CG_ActivityTimeAve = new double [Np];
		CG_SizeActTimeAve = new double [Np];
		CG_ActivityTimeAveDev = new double [Np];
		CG_SizeActTimeAveDev = new double [Np];
		CG_StressTimeAve = new double [Np];
		plateUpdateFailLocations = new int [N];
		cgFailCount = new int[Np];
		CG_Stress = new double [Np];
		failList = new boolean [N];
		for (int i = 0; i < N; i++){
			stress[i] = random.nextFloat()*(tStress-rStress)+rStress;
			stressTimeAve[i] = stress[i];
			failList[i] = false;
		}
		for (int i = 0; i < Np; i++){
			epicenterCount[i] = 0;
			epicenterSize[i] = 0;
			CG_ActivityTimeAve[i] = 0;
			CG_SizeActTimeAve[i] = 0;
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
//		System.out.println("epicenter = " + epicenterSite);
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
//			avSize += 1;
			int nextSiteToFail = checkFailWithRange();

			while (nextSiteToFail >= 0){
//				System.out.println("sitetofail = " + nextSiteToFail);
//				System.out.println("size = " + avSize); 
				failSiteWithRange(nextSiteToFail);
				avSize += 1;
				nextSiteToFail = checkFailWithRange();
			}
		}
//		System.out.println("plateUpdates = " + plateUpdates);
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
			}
		}else{
			plateUpdateFailLocations[epicenterSite] += 1;			
			avSize = 1;
			failSiteWithRange(epicenterSite);
			int nextSiteToFail = checkFailWithRange();
			if(nextSiteToFail >= 0){
				while (nextSiteToFail >= 0){
					failSiteWithRange(nextSiteToFail);
					plateUpdateFailLocations[nextSiteToFail] += 1;
					nextSiteToFail = checkFailWithRange();
					avSize += 1;
				}
			}
		}
		if (avSize >= lowerCutoff){
			int cgSite = findCG_site(epicenterSite);
			epicenterCount[cgSite] +=1;
			epicenterSize[cgSite] += avSize;
		}
		cg_time += dt;
		plateUpdates += 1;
	}

	int checkFailWithRange(){
		int s = siteWithMaxStress();
//		System.out.println("site with max stress = " + s);
		if (stress[s] < tStress) s = -1; 
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
//					 System.out.println("nborsite = " + nborSite);
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
	
	
	public void cgFailCountAdd(){
		for (int i = 0; i < Np; i++){
			int xp = i%Lp;
			int yp = i/Lp;
			
			outer:
			for(int x = xp*dx; x < xp*dx + dx; x ++){
				for(int y = yp*dx; y < yp*dx + dx; y++){
					int site = y*L+x;
					if(plateUpdateFailLocations[site] > 0){
						cgFailCount[i] += 1;
						break outer;
					}
				}
			}

		}
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
			CG_ActivityTimeAveDev[i] = CG_ActivityTimeAve[i] - CG_SpaceAve;
			CG_metric += Math.pow(CG_ActivityTimeAveDev[i],2);
		}
		CG_metric /= (double)Np;
		lastCG_RecordTime = cg_time;
		for (int i = 0; i < Np; i++) epicenterCount[i] = 0;
		return CG_metric; 
	}

	public double calcCG_sizeActMetric(){
		double del_t = cg_time - lastCG_SizeActRecordTime;
		for (int i = 0; i < Np; i++){
			CG_SizeActTimeAve[i] = (lastCG_SizeActRecordTime*CG_SizeActTimeAve[i] + del_t*epicenterSize[i]) / (cg_time);
		}
		double CG_SpaceAve = DoubleArray.sum(CG_SizeActTimeAve)/(double)Np;
		double CG_metric = 0;
		for (int i = 0; i < Np; i++){
			CG_SizeActTimeAveDev[i] = CG_SizeActTimeAve[i] - CG_SpaceAve;
			CG_metric += Math.pow(CG_SizeActTimeAveDev[i] ,2);
		}
		CG_metric /= (double)Np;
		lastCG_SizeActRecordTime = cg_time;
		for (int i = 0; i < Np; i++) epicenterSize[i] = 0;
		return CG_metric; 
	}
	
	public double calcStressMetricPU(){
		double del_t = plateUpdates -lastRecordTimePU;
		for(int i=0; i < N; i++)
			stressTimeAve[i] = (stressTimeAve[i]*(lastRecordTimePU)+ stress[i]*del_t)/(plateUpdates);
		double spaceSum = DoubleArray.sum(stressTimeAve);
		double spaceTimeStressAve = (spaceSum)/(double)(N);
		double metricSum = 0;
		for (int i = 0; i < N; i++) metricSum += Math.pow(stressTimeAve[i] - spaceTimeStressAve, 2);
		double stressMetric = metricSum/(double)N;
		lastRecordTimePU = plateUpdates;
		return stressMetric;
	}

	//fixed to CG in space and time
	public double calcCG_stressMetricPU(){
		double del_t = plateUpdates - lastCG_StressRecordTimePU;
		//		double [] cgStressConfig = findCG_StressConfig();
		CG_Stress = findCG_StressConfig();
		// accumulate the most recent time step into the time ave
		for (int i = 0; i < Np; i++)
			CG_StressTimeAve[i] = (lastCG_StressRecordTimePU*CG_StressTimeAve[i] + del_t*CG_Stress[i]) / (plateUpdates);
		double CG_SpaceAve = DoubleArray.sum(CG_StressTimeAve)/(double)Np;
		double CG_metric = 0;
		//take the space ave
		for (int i = 0; i < Np; i++){
			CG_metric += Math.pow(CG_StressTimeAve[i] - CG_SpaceAve,2);
		}
		CG_metric /= (double)Np;
		lastCG_StressRecordTimePU = plateUpdates;
		return CG_metric; 
	}
	
	public double calcCG_activityMetricPU(){
		double del_t = plateUpdates - lastCG_RecordTimePU;
		for (int i = 0; i < Np; i++){
			CG_ActivityTimeAve[i] = (lastCG_RecordTimePU*CG_ActivityTimeAve[i] + del_t*epicenterCount[i]) / (plateUpdates);
		}
		double CG_SpaceAve = DoubleArray.sum(CG_ActivityTimeAve)/(double)Np;
		double CG_metric = 0;
		for (int i = 0; i < Np; i++){
			CG_ActivityTimeAveDev[i] = CG_ActivityTimeAve[i] - CG_SpaceAve;
			CG_metric += Math.pow(CG_ActivityTimeAveDev[i],2);
		}
		CG_metric /= (double)Np;
		lastCG_RecordTimePU = plateUpdates;
		for (int i = 0; i < Np; i++) epicenterCount[i] = 0;
		
		//calc of breakdown of metric
		am_ave_sq = CG_SpaceAve*CG_SpaceAve*(double)(plateUpdates*plateUpdates);
		double sqSum = 0;
		for (int i = 0; i < Np; i++){
			sqSum += Math.pow(CG_ActivityTimeAve[i],2);
		}
		am_sq_ave = sqSum*(double)(plateUpdates*plateUpdates) / (double)Np;
		
		return CG_metric; 
	}

	public double calcCG_sizeActMetricPU(){
		double del_t = plateUpdates - lastCG_SizeActRecordTimePU;
		for (int i = 0; i < Np; i++){
			CG_SizeActTimeAve[i] = (lastCG_SizeActRecordTimePU*CG_SizeActTimeAve[i] + del_t*epicenterSize[i]) / (plateUpdates);
		}
		double CG_SpaceAve = DoubleArray.sum(CG_SizeActTimeAve)/(double)Np;
		double CG_metric = 0;
		for (int i = 0; i < Np; i++){
			CG_SizeActTimeAveDev[i] = CG_SizeActTimeAve[i] - CG_SpaceAve;
			CG_metric += Math.pow(CG_SizeActTimeAveDev[i] ,2);
		}
		CG_metric /= (double)Np;
		lastCG_SizeActRecordTimePU = plateUpdates;
		for (int i = 0; i < Np; i++) epicenterSize[i] = 0;
		
		//calc of breakdown of metric
		sam_ave_sq = CG_SpaceAve*CG_SpaceAve*(double)(plateUpdates*plateUpdates);
		double sqSum = 0;
		for (int i = 0; i < Np; i++){
			sqSum += Math.pow(CG_SizeActTimeAve[i],2);
		}
		sam_sq_ave = sqSum*(double)(plateUpdates*plateUpdates) / (double)Np;
		
		return CG_metric; 
	}
	
	
	/**
	* The CG fail metric is a metric that use a parameter I call the CG fail paramater.
	* For each plate update, we add one to the fail paramater in each CG block that
	* contains any sites that fail.  The parameter captures more spatial information that just
	* the epicenter locations and hopefully creates clustering that can be detected in the metric.
	*/
//	public double calcCG_failMetric(){
//		double del_t = cg_time - lastCG_RecordTime;
//		for (int i = 0; i < Np; i++){
////			double countWithCutoff;
////			if (epicenterCount[i] < lowerCutoff) countWithCutoff = 0;
////			else countWithCutoff = epicenterCount[i];
////			CG_ActivityTimeAve[i] = (lastCG_RecordTime*CG_ActivityTimeAve[i] + del_t*countWithCutoff) / (time);
//			CG_ActivityTimeAve[i] = (lastCG_RecordTime*CG_ActivityTimeAve[i] + del_t*epicenterCount[i]) / (cg_time);
//		}
//		double CG_SpaceAve = DoubleArray.sum(CG_ActivityTimeAve)/(double)Np;
//		double CG_metric = 0;
//		for (int i = 0; i < Np; i++){
//			CG_metric += Math.pow(CG_ActivityTimeAve[i] - CG_SpaceAve,2);
//		}
//		CG_metric /= (double)Np;
//		lastCG_RecordTime = cg_time;
//		for (int i = 0; i < Np; i++) epicenterCount[i] = 0;
//		return CG_metric; 
//	}

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
	
}
