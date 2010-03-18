package rachele.damage2D;

import java.io.File;

import scikit.jobs.params.Parameters;
import scikit.util.DoubleArray;
import rachele.util.FileUtil;

/**
* Extention of OFC_Lattice.  Probably will replace OFC_Lattice-- 
* to get to OFC only, choose maxFail = 0 
*/
public class OFC_DamageLattice extends AbstractCG_OFC{
 

	double metric0;					//metric measured at time t = 0;
	double lastCG_RecordTime;		//CG time parameter to keep track of last time CG info was recorded
	double lastCG_StressRecordTime; //Basically same as lasCG_RecordTime, but separate to keep calcInverseMeric subroutines independent in case one is not performed
	double lastRecordTime;
	double lastCG_SizeActRecordTime;
	public int lowerCutoff;
	public int meanMaxFail;  		// Mean value of the max no of fails before death
	public int noiseMaxFail; 		// Max fail value = meanMaxFail +- eta where max eta = maxFailNoise
	public int meanHealTime;			// average time before sites heal
	public int noiseHealTime;		// heal time for site = aveHealTime +- eta where max eta = noiseHealTime
	public boolean deadMode;
	public int noDeadSites;
	public int nextSiteToFail;
	public double initPercentDead;
	public double healProb; 		// Probability of healing a dead site when mode = healProb
	int maxNbors;
	
	public double [] stressTimeAve;
	public double [] CG_ActivityTimeAve;
	public double [] CG_SizeActTimeAve;
	public double [] CG_StressTimeAve;
	
	public double [] CG_Stress;
	boolean [] failList;  // Temporary list of failures
	public int [] noFails;  // No of failures at each site: once a site is dead -> = -1;
	public int [][]nborList;  // List of neighbors for each site: nborList[site = 0..N-1][nbor index = 0..maxNbors]  This is gonna be really big
	public int [] noNborsForSite;  // No of nbors at each site (for when different for each site.)
	public int [] maxNoFails; //  Max no of fails for each site
	public int [] healTime; // Time before site heals
	public String interactionTopology;
	public int [] plateUpdateFailLocations;  // Record all fail sites for one plate update
	public int [] cgFailCount; //Counts 1 in a cg block if there is any site that fails in that block.
	public int [] epicenterSize;
	
	public int noLiveSites; //For frozen damage mode only
	public int [] liveSiteLocs; // size = noLiveSites
	
	public String infoFileName;
	boolean deadStrip;
	boolean randomBlockDamage;
	String mode;
	String damageType;
	double maxPercentDamage;
	double [] cgBlockDamage;
	
	public OFC_DamageLattice(Parameters params, String mode) {
		infoFileName = params.sget("Data Dir") + File.separator + "info.txt";
		FileUtil.initFile(infoFileName, params);
		FileUtil.printlnToFile(infoFileName, "Runlog for OFC_DamageLattice");
		setParameters(params, mode);
		initLattice(mode);
	}
	
	public void setParameters(Parameters params, String mode) {
	
		FileUtil.printlnToFile(infoFileName, "Setting params");
		meanMaxFail = params.iget("Mean Max Failures");
		noiseMaxFail = params.iget("Failures Max Noise");
		if(mode == "healFrozen"){
			meanHealTime = params.iget("Mean Heal Time");
			FileUtil.printlnToFile(infoFileName, "mean heal time = " , meanHealTime);
			noiseHealTime = params.iget("Heal Time Noise");
		}else if(mode == "healProb"){
			healProb = params.fget("Heal Probability");
		}
		maxPercentDamage = params.fget("Max Percent Damage");
//		if (meanMaxFail == 0) deadMode = false;
//		else deadMode = true;
		deadMode = true;
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

		interactionTopology = params.sget("Interaction");
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
		initPercentDead = params.fget("Init Percent Dead");
		damageType = params.sget("Type of Damage");
	}
	
	public void initLattice(String importMode){
		System.out.println("Initializing lattice");
		stress = new double [N];
		stressTimeAve = new double [N];
		epicenterCount = new int [Np];
		epicenterSize = new int [Np];
		CG_ActivityTimeAve = new double [Np];
		CG_SizeActTimeAve = new double [Np];
		CG_StressTimeAve = new double [N];
		plateUpdateFailLocations = new int [N];
		cgFailCount = new int[Np];
		CG_Stress = new double [Np];
		failList = new boolean [N];
		noFails = new int [N];
		maxNoFails = new int [N];
		mode=importMode;
		if(mode == "healFrozen")
			healTime = new int [N];
		noNborsForSite = new int [N];
		maxNbors = findMaxNoNbors();	
		System.out.println("max no nbors = " + maxNbors);
		FileUtil.printlnToFile(infoFileName, "max no nbors = " , maxNbors);
		if(interactionTopology == "Circle")
			nborList = new int [N][maxNbors];						

		findNbors();
		noDeadSites = 0;
		noLiveSites = N;
		
		System.out.println("Setting max fails, heal times, and noFails");
		FileUtil.printlnToFile(infoFileName, "Setting max fails, heal times, and noFails");
		for (int i = 0; i < N; i++){
			stress[i] = random.nextFloat()*(tStress-rStress)+rStress;
			stressTimeAve[i] = stress[i];
			failList[i] = false;
			maxNoFails[i] = meanMaxFail;	
			CG_StressTimeAve[i] = 0;
//			FileUtil.printlnToFile(infoFileName, "deadmode = ", deadMode);
//			if(deadMode){
////				maxNoFails[i] = meanMaxFail;									
//				maxNoFails[i] = meanMaxFail + (int)(random.nextDouble()*(double)(2*noiseMaxFail+1))-noiseMaxFail;
//				if(mode == "healFrozen") healTime[i] = meanHealTime +  (int)(random.nextDouble()*(double)(2*noiseHealTime+1))-noiseHealTime;
////				if (i<L) System.out.println("no of fails = " + maxNoFails[i]);
//			}else{
//				maxNoFails[i] = meanMaxFail;									
//			}
////			noFails[i] = 0;
		}
		
		//Choose one type of damage
		if(damageType=="Random"){
			for(int i = 0; i < N; i++){
				if(random.nextDouble() > initPercentDead){
					setSite(i);
				}else{
					killSite(i);
				}
			}
		}else if(damageType=="Dead Strip"){
			int width = (int)(L*initPercentDead);
			for(int i  = 0; i < N; i++){
				int y = i%L;
				if(y < width){
					killSite(i);
				}else{
					setSite(i);
				}
			}
		}else if (damageType=="Random + Dead Strip"){
			System.out.println("Error!");
		}else if(damageType=="Random Blocks"){
			int damageBlockSize = 8;
			FileUtil.printlnToFile(infoFileName, "Damage block size", damageBlockSize);
			double [] damageBlockDamage = new double [damageBlockSize*damageBlockSize];
			for (int j = 0; j < damageBlockSize*damageBlockSize; j++){
				//assign a random amount of damage
				//half of blocks have no damage
				damageBlockDamage[j] = random.nextGaussian()*0.2;
				System.out.println("damage = " + j + " "+ damageBlockDamage[j]);
			}
			for (int i = 0; i < N; i++){
				int block = findCG_site(i, damageBlockSize);
				if (random.nextDouble() < damageBlockDamage[block]){
					killSite(i);
				}else{
					setSite(i);
				}
			}
			//find block damage for cged blocks
			cgBlockDamage = new double [Np];
			for (int i = 0; i < N; i++)
				if (noFails[i]<0) cgBlockDamage[findCG_site(i,dx)] +=1.0;	
			for (int j = 0; j < Np; j++)
				cgBlockDamage[j] /= (double)(dx*dx);
			
		}else{
			System.out.println("Error!");	
		}
			
		liveSiteLocs = new int [noLiveSites];

		for (int i = 0; i < Np; i++){
			epicenterCount[i] = 0;
			CG_ActivityTimeAve[i] = 0;
		}
		metric0 = calcMetric0();
		cg_time = 0;
		plateUpdates = 0;

	}
	
	void killSite(int site){
		noDeadSites+=1;
		noLiveSites-=1;
		//assign a heal time uniformly between -1 and -healtime+1
		if(mode == "healFrozen"){
			noFails[site] = -((int)(random.nextDouble()*(double)(healTime[site]-1))+1);
			System.out.println(site + " dead healTime = " + healTime[site] + " noFails = " + noFails[site]);
		}else{
			noFails[site] = -1;
		}
	}
	
	void setSite(int site){
		//the site is alive
		//assign no of fails uniformly between 0 and maxNoFails - 1
		noFails[site] = (int)(random.nextDouble()*((double)maxNoFails[site]+1.0));
//		System.out.println(i + " alive maxNoFails = " + maxNoFails[i] + " noFails = " + noFails[i]);
	}
	
	public void findNbors(){
		System.out.println("Finding neighbors ");
		if(fullyConnected){
			System.out.println("Fully connected");
		}else{
			if (interactionTopology == "Circle"){
				findCircleNbors();
			}else if(interactionTopology == "Square"){
				System.out.println("Need circle code for find SquareNbors");
				
			}
			//write the rest of this for other topologies	
			//}			
		}
	}
	
	public int findMaxNoNbors(){
		System.out.println(interactionTopology);
		int noNbors = 0;
		if (interactionTopology == "Circle"){
			noNbors = findNoCircleNbors();
			for (int i = 0; i < N; i++)	noNborsForSite[i] = noNbors; 
		}else if (interactionTopology == "Square"){	
			noNbors = (2*R+1)*(2*R+1)-1;
			for (int i = 0; i < N; i++)	noNborsForSite[i] = noNbors; 
		}else if (interactionTopology == "Fully Connected"){	
			noNbors = N-1;
			for (int i = 0; i < N; i++)	noNborsForSite[i] = noNbors; 
			System.out.println("nbors = " + noNbors);
		}else if (interactionTopology == "Small World"){	
			System.out.println("Need small world code for findMaxNbors");
		}		
		if (noNbors == 0) System.out.println("Error in findMaxNoNbors");
		return noNbors;
	}
	
	public int findNoCircleNbors(){
		int count = 0;
		 for(int dy = -R; dy <= R; dy++){
			 for(int dx = -R; dx <= R; dx++){
				 double distance = Math.sqrt(dx*dx + dy*dy);
				 if (distance <= R){
					 count += 1;
				 }
			 }
		 }
		 return count;
	}
	
	public void findCircleNbors(){
		int nborIndex = 0;
		for (int s = 0; s < N; s++){
			nborIndex = 0;
			int x = s%L;
			int y = s/L;
			for(int dy = -R; dy <= R; dy++){
				for(int dx = -R; dx <= R; dx++){
					double distance = Math.sqrt(dx*dx + dy*dy);
					if (distance <= R){
						int xx = (x+dx+L)%L;
						int yy = (y+dy+L)%L;
						int nborSite = yy*L+xx;
						nborList[s][nborIndex] = nborSite;
						nborIndex += 1;
					}
				}
			}
		}
		if (nborIndex != (maxNbors)) System.out.println("Problem in findCircleNbors");
	}
	
	public void initEquilibrate(int maxPlateUpdates){
			plateUpdates = -maxPlateUpdates;
	}
	
	public void equilibrate(){

		bringToFailure();
		
		failSite(epicenterSite);  //Just distributes stress to each neighbor
//		System.out.println("Site failed");
		int nextSiteToFail = checkFail();  //each checkFail loops through lattice once.
//		avSize = 1;
		while (nextSiteToFail >= 0){
			failSite(nextSiteToFail);
			nextSiteToFail = checkFail();
//			System.out.println("avSize = " + avSize + " nextSiteToFail = " + nextSiteToFail);
			avSize += 1;
		}
//		System.out.println("av size = " + avSize);
		
//		if(fullyConnected){
//			failList[epicenterSite] = true;
//			boolean failAgain = checkFailList();
//			while(failAgain){
//				avSize += fail();
//				failAgain = checkFailList();
//			}
//		}else{
//			failSiteWithList(epicenterSite);  //Just distributes stress to each neighbor
//			int nextSiteToFail = checkFail();  //each checkFail loops through lattice once.
//			avSize = 1;
//			while (nextSiteToFail >= 0){
//				failSiteWithList(nextSiteToFail);
//				nextSiteToFail = checkFail();
//				avSize += 1;
//			}
//		}
		plateUpdates += 1;

	}
	
	public void bringToFailure(){
		epicenterSite = siteWithMaxStress();	
		double stressAdd = tStress - stress[epicenterSite];
		if(stressAdd < 0) System.out.println("Error: stress already above failure");
		for (int i = 0; i < N; i++){
			if(noFails[i] >=0) stress[i] += stressAdd;
			plateUpdateFailLocations[i] = 0;
		}
		plateUpdateFailLocations[epicenterSite] = 1;
		avSize = 1;
	}
	
	/**
	* One step in OFC only mode- no damage, no fails, no healing
	*/
	public void ofcStep(){
		bringToFailure();
		dt=1;
	
		failSite(epicenterSite);
//		failSiteWithList(epicenterSite);
		int nextSiteToFail = checkFail();
		while(nextSiteToFail >= 0){
//			failSiteWithList(nextSiteToFail);
			failSite(nextSiteToFail);
			nextSiteToFail = checkFail();
			avSize += 1;
		}

		cgEpicenterSizeAccum();
		cg_time += dt;
		plateUpdates += 1;
	}
	
	
	public void cgEpicenterSizeAccum(){
		if (avSize >= lowerCutoff){
			int cgSite = findCG_site(epicenterSite, dx);
			epicenterCount[cgSite] +=1;
			epicenterSize[cgSite] += avSize;
		}
		
	}
	/**
	* One step in Damage OFC mode- accumulate no of fails, kill sites > maxFail
	*/
	public void damagePreStep(){
		dt = 1;
		bringToFailure();
		failSite(epicenterSite);
		countFailAndCheckForDeath(epicenterSite);
		nextSiteToFail = checkNextFailSite();
		plateUpdates += 1;	
		cg_time += dt;
	}

	public void failSite(int s){
//		System.out.println(fullyConnected);
		if(fullyConnected){
			failSiteForFC(s);
		}else{
			failSiteWithList(s);
		}
	}
	
	public void damageStepIter(){
		failSite(nextSiteToFail);
		countFailAndCheckForDeath(nextSiteToFail);
		avSize += 1;
		nextSiteToFail = checkNextFailSite();
	}
	
	/**
	* One step in Damage/Healing Mode- sites fail, die and are healed.
	*/
	public void healPreStep(){
		dt = 1;
		bringToFailure();
		failSite(epicenterSite);
		countFailAndCheckForDeath(epicenterSite);
		nextSiteToFail = checkNextFailSite();
		plateUpdates += 1;
		cg_time += dt;
	}
	
	public void healStepIter(){
		failSite(nextSiteToFail);
		countFailAndCheckForDeath(nextSiteToFail);
		avSize += 1;
		nextSiteToFail = checkNextFailSite();
	}
	
	public void checkForFrozenHealing(){
		int countHeal = 0;
		for (int s = 0; s < N; s++){
//			System.out.println("noFails = " + s + " " +noFails[s]);
			if(noFails[s] < 0){
//				int testNo=Math.abs(noFails[s] - (plateUpdates));
//				System.out.println("test no = " + testNo + " healTime[s] = " + healTime[s]);
//				if(testNo >= healTime[s]){
				//heal the site
				if(noFails[s] <= -healTime[s]){
					noFails[s] = 0;
					noDeadSites -= 1;
					for (int i = 0; i < maxNbors; i++){
						noNborsForSite[nborList[s][i]] += 1;
					}
					countHeal+=1;
				}else{
					noFails[s] -= 1;
				}
			}
		}
//		System.out.println("count heal = " + countHeal);
	}
	
	public void checkForProbHealing(){
		int countHeal = 0;
		for (int s = 0; s < N; s++){
			if(noFails[s] < 0){
				if(random.nextDouble() <= healProb){
					noFails[s] = 0;
					noDeadSites -= 1;
					for (int i = 0; i < maxNbors; i++){
						noNborsForSite[nborList[s][i]] += 1;
					}
					countHeal+=1;
				}
			}
		}
	}
	
	public void countFailAndCheckForDeath(int s){
		noFails[s] += 1;
		if (noFails[s] >= maxNoFails[s]){
//			noFails[s] = -Math.abs(plateUpdates);
			noFails[s] = -1;
			noDeadSites += 1;
//			stress[s] = 0;
			if(fullyConnected){
				noNbors -= 1;
			}else{
				for (int i = 0; i < maxNbors; i++){
					noNborsForSite[nborList[s][i]] -= 1;
				}
			}
		}
	}
	
	public int checkNextFailSite(){
		int s = siteWithMaxStress();
		if (stress[s] < tStress) s = -1; 
		return s;
	}
	
//	/**
//	 *  This method should not be used when equilibrating in earthquake mode.
//	 *  Do not combine this with failSiteWithRange.
//	 */
//	void countFail(int siteLocation){
//		if(noFails[siteLocation] >= maxNoFails[siteLocation]) noFails[siteLocation] = -1;
//		else noFails[siteLocation] += 1;
//	}
	
	int checkFail(){
		int s = siteWithMaxStress();
		if (stress[s] < tStress) s = -1; 
//		System.out.println("Stress of max Site = " + stress[s] + " at site " + s);
		return s;
	}
	
	public void failSiteForFC(int s){
		double resStressWithNoise = calcResNoise(s);
		double stressPerNbor = ((1-dissParam)*(stress[s] - resStressWithNoise))/(double)noNbors;
		stress[s] = resStressWithNoise;
		//distribute stress to the alive nbors
		int checkAliveCount = 0;
//		System.out.println(maxNbors);
		for (int i = 0; i <= maxNbors; i++){
			if(noFails[i] >= 0){
				stress[i] += stressPerNbor;
				checkAliveCount += 1;
			}
		}
		stress[s] -= stressPerNbor;
		checkAliveCount -= 1;
		plateUpdateFailLocations[s] += 1;
//		System.out.println("alive nbors " + noNbors + " " + checkAliveCount);
	}
	
	public void failSiteWithList(int s){
		double resStressWithNoise = calcResNoise(s);
		double stressPerNbor = ((1-dissParam)*(stress[s] - resStressWithNoise))/(double)noNborsForSite[s];
		stress[s] = resStressWithNoise;
		//distribute stress to the alive nbors
		int checkAliveCount = 0;
		for (int i = 0; i < maxNbors; i++){
			int nborSite = nborList[s][i];
			if(noFails[nborSite] >= 0){
				stress[nborSite] += stressPerNbor;
				checkAliveCount += 1;
			}
		}
		plateUpdateFailLocations[s] += 1;
//		System.out.println("alive nbors " + noNborsForSite[s] + " " + checkAliveCount);
	}
	
//	public double calcInverseMetric(){
//		double del_t = cg_time -lastRecordTime;
//		for(int i=0; i < N; i++)
//			stressTimeAve[i] = (stressTimeAve[i]*(lastRecordTime)+ stress[i]*del_t)/(cg_time);
//		double spaceSum = DoubleArray.sum(stressTimeAve);
//		double spaceTimeStressAve = (spaceSum)/(double)(N);
//		double metricSum = 0;
//		for (int i = 0; i < N; i++) metricSum += Math.pow(stressTimeAve[i] - spaceTimeStressAve, 2);
//		double inverseMetric = (double)(N)*metric0/metricSum;
//		lastRecordTime = cg_time;
//		return inverseMetric;
//	}

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

	/**
	* Calculate the coarse grained stress metric using all sites (dead or alive).  Fixed to CG in space and time.
	*/
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
	
	/**
	* Calculate the stress metric using only sites that are alive
	*/
	public double calcLiveStressMetric(){
		double del_t = cg_time -lastRecordTime;
		for(int i=0; i < N; i++)
			stressTimeAve[i] = (stressTimeAve[i]*(lastRecordTime)+ stress[i]*del_t)/(cg_time);
		double spaceSum = 0.0;
		int noAlive = 0;
		for(int i=0; i < N; i++){
			if(noFails[i] >= 0){
				noAlive += 1;
				spaceSum += stressTimeAve[i];
			}
		}
		double spaceTimeStressAve = (spaceSum)/(double)(noAlive);
		double metricSum = 0;
		for (int i = 0; i < N; i++){
			if(noFails[i] >= 0) metricSum += Math.pow(stressTimeAve[i] - spaceTimeStressAve, 2);
		}
		double stressMetric = metricSum/(double)noAlive;
//		System.out.println("no alive = " + noAlive);
		lastRecordTime = cg_time;
		return stressMetric;
	}
	
	/**
	* Calculate the coarse grained stress metric using only sites that are alive.
	*/
	public double calcLiveCG_stressMetric(){
		double del_t = cg_time - lastCG_StressRecordTime;
		CG_Stress = findCG_StressConfigNoDead();
		int noLiveCgedBlocks = 0;
		for (int i = 0; i < Np; i++){
			if(CG_Stress[i]!=0){
				CG_StressTimeAve[i] = (lastCG_StressRecordTime*CG_StressTimeAve[i] + del_t*CG_Stress[i]) / (cg_time);
				noLiveCgedBlocks+=1;
			}
		}
		double CG_SpaceAve = DoubleArray.sum(CG_StressTimeAve)/(double)noLiveCgedBlocks;
		double CG_metric = 0;
		for (int i = 0; i < Np; i++){
			if(CG_Stress[i]!=0){
				CG_metric += Math.pow(CG_StressTimeAve[i] - CG_SpaceAve,2);
			}
		}
		CG_metric /= (double)noLiveCgedBlocks;
		lastCG_StressRecordTime = cg_time;
		return CG_metric; 
	}

	/**
	* Find the average stress on each coarse grained block.  
	* If sites are dead, they do not contribute to the average.
	* If all sites in a coarse grained block are dead, then the 
	* coarse grained block itself is dead- the average stress is zero
	*  and it is not included in the 
	* final coarse grained stress metric calculation.  
	* 
	* This is pretty weird and needs to be thought about if used for something important:
	* why should one block with one live site be worth the same as one block with all live sites.
	* Unequal weighting because we take a average of stresses.
	*/
	public double [] findCG_StressConfigNoDead(){
		double [] cgSpaceStressConfig = new double [Np];
		int [] noLiveSitesInBlock = new int [Np];
		for(int i = 0; i < N; i++){
			if(noFails[i] >= 0){//only add stress for live sites
				int cgSite = findCG_site(i, dx);
				cgSpaceStressConfig[cgSite] += stress[i];
				noLiveSitesInBlock[cgSite] += 1;
			}
		}
		for(int cgSite = 0; cgSite < Np; cgSite++){
			if(noLiveSitesInBlock[cgSite]!=0)
				cgSpaceStressConfig[cgSite] /= (double) noLiveSitesInBlock[cgSite];
		}
		return cgSpaceStressConfig;
	}
	
	public double [] findCG_StressConfig(){
		double [] cgSpaceStressConfig = new double [Np];
		for(int i = 0; i < N; i++){
			int cgSite = findCG_site(i, dx);
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


	
	/**
	* Calculate the coarse grained activity metric using all sites (dead or alive).
	*/
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

	/**
	* Calculate the coarse grained activity metric using only alive coarse grained sites.
	* Add up all hits in each coarse grained box.
	* If there are boxes with no hits, do not include them in the 
	* metric calculation.
	* We calculate the activity time average at each step, but we do not include
	* boxes in the metric calculation that do not have any hits.
	* 
	* There is a problem for dt too small- get large no: divide by zero?
	*/
	public double calcLiveCG_activityMetric(){
		double del_t = cg_time - lastCG_RecordTime;
		for (int i = 0; i < Np; i++){
			CG_ActivityTimeAve[i] = (lastCG_RecordTime*CG_ActivityTimeAve[i] + del_t*epicenterCount[i]) / (cg_time);
		}
		int noCgBoxesWithHits = 0;
		double CG_SpaceAve = 0;
		for (int i = 0; i < Np; i++){
			if(epicenterCount[i] !=0){
				noCgBoxesWithHits += 1;
				CG_SpaceAve += CG_ActivityTimeAve[i];
			}
		}
//		System.out.println(" no of cg boxes with hits " + noCgBoxesWithHits);
		CG_SpaceAve /=(double) noCgBoxesWithHits;
//		double CG_SpaceAve = DoubleArray.sum(CG_ActivityTimeAve)/(double)Np;
		double CG_metric = 0;

		for (int i = 0; i < Np; i++){
			if(epicenterCount[i] !=0){
				CG_metric += Math.pow(CG_ActivityTimeAve[i] - CG_SpaceAve,2);
			}
		}
		CG_metric /= (double)noCgBoxesWithHits;
		lastCG_RecordTime = cg_time;
		for (int i = 0; i < Np; i++) epicenterCount[i] = 0;
		return CG_metric; 
	}
	
	/**
	* Calculate the coarse grained size-activity metric using alive sites only.
	*/
	public double calcLiveCG_sizeActMetric(){
		double del_t = cg_time - lastCG_SizeActRecordTime;
		for (int i = 0; i < Np; i++){
			CG_SizeActTimeAve[i] = (lastCG_SizeActRecordTime*CG_SizeActTimeAve[i] + del_t*epicenterSize[i]) / (cg_time);
		}
		int noCgBoxesWithHits = 0;
		double CG_SpaceAve = 0;
		for (int i = 0; i < Np; i++){
			if(epicenterSize[i] !=0){
				noCgBoxesWithHits += 1;
				CG_SpaceAve += CG_SizeActTimeAve[i];
			}
		}
		CG_SpaceAve /=(double) noCgBoxesWithHits;
		double CG_metric = 0;
		for (int i = 0; i < Np; i++){
			if(epicenterSize[i] !=0){
				CG_metric += Math.pow(CG_SizeActTimeAve[i] - CG_SpaceAve,2);
			}
		}
		CG_metric /= (double)noCgBoxesWithHits;
		lastCG_SizeActRecordTime = cg_time;
		for (int i = 0; i < Np; i++) epicenterSize[i] = 0;
		return CG_metric; 
	}
	
	public double calcUndamagedCG_activityMetric(){
		double blockDamageTolerance = maxPercentDamage;
		double del_t = cg_time - lastCG_RecordTime;
		for (int i = 0; i < Np; i++){
			CG_ActivityTimeAve[i] = (lastCG_RecordTime*CG_ActivityTimeAve[i] + del_t*epicenterCount[i]) / (cg_time);
		}
		int noCgBoxesWithHits = 0;
		double CG_SpaceAve = 0;
		for (int i = 0; i < Np; i++){
			if(cgBlockDamage[i] < blockDamageTolerance){
				noCgBoxesWithHits += 1;
				CG_SpaceAve += CG_ActivityTimeAve[i];
			}
		}
//		System.out.println(" no of cg boxes with hits " + noCgBoxesWithHits);
		CG_SpaceAve /=(double) noCgBoxesWithHits;
//		double CG_SpaceAve = DoubleArray.sum(CG_ActivityTimeAve)/(double)Np;
		double CG_metric = 0;

		for (int i = 0; i < Np; i++){
			if(cgBlockDamage[i] < blockDamageTolerance){
				CG_metric += Math.pow(CG_ActivityTimeAve[i] - CG_SpaceAve,2);
			}
		}
		CG_metric /= (double)noCgBoxesWithHits;
//		System.out.println("max no nbors = " + noCgBoxesWithHits);
		lastCG_RecordTime = cg_time;
		for (int i = 0; i < Np; i++) epicenterCount[i] = 0;
		return CG_metric; 

	}
	
	/**
	* Calculate the coarse grained size-activity metric using alive sites only.
	*/
	public double calcUndamagedCG_sizeActMetric(){
		double blockDamageTolerance = maxPercentDamage;
		double del_t = cg_time - lastCG_SizeActRecordTime;
		for (int i = 0; i < Np; i++){
			CG_SizeActTimeAve[i] = (lastCG_SizeActRecordTime*CG_SizeActTimeAve[i] + del_t*epicenterSize[i]) / (cg_time);
		}
		int noCgBoxesWithHits = 0;
		double CG_SpaceAve = 0;
		for (int i = 0; i < Np; i++){
			if(cgBlockDamage[i] < blockDamageTolerance){
				noCgBoxesWithHits += 1;
				CG_SpaceAve += CG_SizeActTimeAve[i];
			}
		}
		CG_SpaceAve /=(double) noCgBoxesWithHits;
		double CG_metric = 0;
		for (int i = 0; i < Np; i++){
			if(cgBlockDamage[i] < blockDamageTolerance){
				CG_metric += Math.pow(CG_SizeActTimeAve[i] - CG_SpaceAve,2);
			}
		}
		CG_metric /= (double)noCgBoxesWithHits;
		lastCG_SizeActRecordTime = cg_time;
		for (int i = 0; i < Np; i++) epicenterSize[i] = 0;
		if(plateUpdates <= 2*dt) FileUtil.printlnToFile(infoFileName, "No of coarse grained boxes with hits = ", noCgBoxesWithHits);
		return CG_metric; 
	}
	
	
	/**
	* Calculate the coarse grained size-activity metric using all sites (dead or alive).
	*/
	public double calcCG_sizeActMetric(){
		double del_t = cg_time - lastCG_SizeActRecordTime;
		for (int i = 0; i < Np; i++){
			CG_SizeActTimeAve[i] = (lastCG_SizeActRecordTime*CG_SizeActTimeAve[i] + del_t*epicenterSize[i]) / (cg_time);
		}
		double CG_SpaceAve = DoubleArray.sum(CG_SizeActTimeAve)/(double)Np;
		double CG_metric = 0;
		for (int i = 0; i < Np; i++){
			CG_metric += Math.pow(CG_SizeActTimeAve[i] - CG_SpaceAve,2);
		}
		CG_metric /= (double)Np;
		lastCG_SizeActRecordTime = cg_time;
		for (int i = 0; i < Np; i++) epicenterSize[i] = 0;
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
				noNbors -= 1;
			}
		}
		for(int i = 0; i < N; i++){
			stress[i] += excessStressSum;
		}
		return noFailed;
	}

	int failWithDamage(){
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

	double calcStressPerNborWithDamage(int s, double rStressWithNoise){
		noNbors = N - noDeadSites - 1;
		double stressPerNbor = ((1-dissParam)*(stress[s] - rStressWithNoise))/(double)noNbors;
		return stressPerNbor;
	}

	
	
	
}
