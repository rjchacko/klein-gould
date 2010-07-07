package rachele.damage2D.multidx;
import scikit.dataset.Histogram;

import rachele.util.FileUtil;
import rachele.util.MathTools;
import scikit.jobs.params.Parameters;

public class FrozenDamageLattice extends AbstractOFC_Multidx{

	public boolean [] aliveLattice;  // No of failures at each site: once a site is dead -> = -1;
	public int [] epicenterSize;
	public double alphaDissRate;
	public double deadDissRate;
	public double alphaDiss;
	public double deadDiss;
	public String boundaryConditions;
	public Histogram alpha_iHist = new Histogram(0.001);
	
	public int noLiveSites; 
	
	int noDeadSites;
	int nextSiteToFail;
	int maxNbors;
	int [] noNborsForSite;  // No of nbors at each site (for when different for each site.)
	int [][]nborList;  // List of neighbors for each site: nborList[site = 0..N-1][nbor index = 0..maxNbors]  This is gonna be really big
	
	double maxPercentDamage;
	double [] cgBlockDamage;
	double [] fracDeadNbors;
	double alphaDissAccum;
	double deadDissAccum;
	double rateCt;
	
	boolean [] setLattice;
	boolean deadStrip;
	boolean randomBlockDamage;
	boolean deadDissipation;

	String damage;
	String infoFileName;
	
	public FrozenDamageLattice(Parameters params, String name) {
		infoFileName = name;
		FileUtil.initFile(infoFileName, params);
		FileUtil.printlnToFile(infoFileName, "# Runlog for OFC_DamageLattice");
		setParameters(params);
		initLattice(params);
	}
	
	public void setParameters(Parameters params) {
	
		FileUtil.printlnToFile(infoFileName, "# Setting params");
		if(params.sget("Dead dissipation?")=="Yes") deadDissipation = true;
		else deadDissipation = false;
		System.out.println("Dead dissipation = " + deadDissipation);
		pow = params.iget("Size Power");
		L = (int) Math.pow(2, pow);
		N = L*L;
		R = params.iget("R");
		noNbors = findCircleNbors(R);
		random.setSeed(params.iget("Random Seed"));	
		tStress = 1.0;
		alpha = params.fget("Dissipation Param");   // = alpha
		rStress = params.fget("Residual Stress");
		maxResNoise = params.fget("Res. Max Noise");
		damage = params.sget("Type of Damage");
		deadDissRate = 0.0;
		alphaDissRate = 0.0;
		rateCt = 0;
		boundaryConditions = params.sget("Boundary Conditions");
	}
	
	void initLattice(Parameters params){
		
		p("Initializing lattice");
		stress = new double [N];
		epicenterCount = new int [N];
		epicenterSize = new int [N];
		aliveLattice = new boolean [N];
		setLattice = new boolean [N];
		noNborsForSite = new int [N];
		fracDeadNbors = new double [N];
		maxNbors = findNoCircleNbors();	
		System.out.println("max no nbors = " + maxNbors);
		FileUtil.printlnToFile(infoFileName, "# max no nbors = " , maxNbors);
		p("Allocating nbor list...");
		nborList = new int [N][maxNbors];	
		
		p("Finding neighbors...");
		findCircleNbors();
		noDeadSites = 0;
		noLiveSites = N;
		
		p("Finding stress...");
		for (int i = 0; i < N; i++){
			stress[i] = random.nextFloat()*(tStress-rStress)+rStress;	
		}
		
		p("Setting damage...");
		setDamage(damage, params.iget("Dead Parameter"), params.fget("Init Percent Dead"), params.iget("Number Dead"));
		
		p("Making neighbor lists...");
		makeNborLists();
		
		p("Calculating alpha prime and variance");
		double [] ap = calcAlphaP();
		System.out.println("Alpha prime = " + ap[0] + ", var = " + ap[1]);
		FileUtil.printlnToFile(infoFileName, "Alpha' = " , ap[0], " var = " , ap[1]);
		
		double percentSitesDamaged = (double)noDeadSites/(double)N;
		FileUtil.printlnToFile(infoFileName, "# No of sites damaged = " , noDeadSites);
		
//		params.set("Dead Parameter", noDeadSites);
		FileUtil.printlnToFile(infoFileName, "# Percent of sites damaged = ", percentSitesDamaged);
		System.out.println("Percent of sites damaged = "+ percentSitesDamaged);
		params.set("Percent Damage", percentSitesDamaged);

		for (int i = 0; i < N; i++){
			epicenterCount[i] = 0;
			epicenterSize[i] = 0;
		}
		cg_time = 0;
		plateUpdates = 0;

		p("Lattice init complete.");
	}
	
	void killSite(int site){
		noDeadSites+=1;
		noLiveSites-=1;
		aliveLattice[site] = false;
		stress[site] = 0.0;
	}
	
	void setSite(int site){
		aliveLattice[site] = true;
	}
	
	void setDamage(String damageType, int deadParam, double initPercentDead, int noDead){
//		int noDeadToPlace = noDead;
//		int deadScale = deadParam;
		//set alive lattice
		for(int i = 0; i < N; i ++){
			aliveLattice[i] = true;
		}
		if(damageType=="Random") setRandom(initPercentDead);
		else if(damageType == "Place Random Dead")	setPlaceDeadRandom(noDead);
		else if(damageType=="Dead Strip") setDeadStrip(noDead);
		else if(damageType=="Random Blocks") setRandomBlockDamage(deadParam);
		else if(damageType == "Dead Block") setDeadBlock(noDead);
		else if(damageType == "Cascade") setCascadeDamage(initPercentDead);
		else if(damageType == "Cascade Random") setCascadeRandom(initPercentDead);
		else if(damageType == "Dead Rectangle") setDeadRectangle(noDead);
		else if(damageType == "Dead Blocks") setDeadBlocks(deadParam, initPercentDead);
		else if(damageType == "Place Dead Blocks") placeDeadBlocks(deadParam, noDead);		
		else System.out.println("Error!");	
		
	}
	
	private void setDeadRectangle(int noDead) {
		int aspectRatio = 3;// width over height
		int w = (int)Math.sqrt(aspectRatio*noDead);
		if(w > L) System.out.println("Aspect Ratio error!");
		int h = noDead/w;
		int initialBlock = h*w;
		for (int x = 0; x < w; x++){
			for (int y = 0; y < h; y++){
				int s = y*L+x;
				killSite(s);
			}
		}
		int extra = noDead - initialBlock;
		if(extra > h){
			for(int i = 0; i < h; i++) killSite(i*L + w);
			for (int i = 0; i <= extra-h; i++) killSite(h*L + i);
		}else{
			for (int i = 0; i < extra; i++) killSite(i*L + w);
		}
		System.out.println(" dead = " + noDead + " no in block so far = "  + initialBlock);
		
	}

	private void setCascadeRandom(double p) {
		System.out.println("Setting Cascade Random Damage");
		int maxPow = pow-2;
		for (int i = 0; i < N; i++) setLattice[i] = false;

		for (int ip = maxPow; ip >= 0; ip--){
			int dx = (int)Math.pow(2 ,ip);
			int Lp = L/dx;
			int Np = Lp*Lp;
			System.out.println("dx = " + dx + " Lp = " + Lp);	
			for (int i = 0; i < Np; i++){
				// set blocks with probability if they are not already set
				if(setBlockTest(i, ip)==false){
					if (random.nextDouble()<p) 
						setBlock(dx, i);
				}
			}
				
		}
		
	}
	
	boolean setBlockTest(int blockNo, int ip){
		boolean ret = false;
		int dx = (int)Math.pow(2, ip);
		int Lp = L/dx;
		int xp = blockNo%Lp;
		int yp = blockNo/Lp;
		int x = dx*xp;
		int y = dx*yp;
		int site = y*L+x;
		if(setLattice[site]) ret = true;
		return ret;
	}
	
	void setBlock(int dx, int block){
		//all sites in block
		int Lp = L/dx;
		int xp = block%Lp;
		int yp = block/Lp;
		double p =random.nextDouble();
		for (int y = yp*dx; y < yp*dx+dx; y++){
			for (int x = xp*dx; x < xp*dx + dx; x++){
				int site = y*L+x;
				if(random.nextDouble()<p) killSite(site);
				setLattice[site] = true;
			}
		}
	}
	
	void setCascadeDamage(double p){
		System.out.println("Setting Cascading Damage");
		int maxPow = pow-2;
		//largest lattice size is 8 x 8

		for (int ip = maxPow; ip >= 0; ip--){
			int dx = (int)Math.pow(2 ,ip);
			int Lp = L/dx;
			int Np = Lp*Lp;
			System.out.println("dx = " + dx + " Lp = " + Lp);	
			for (int i = 0; i < Np; i++){
				// kill alive blocks with probability p
				if(liveBlockTest(i, ip)){
					if (random.nextDouble()<p) 
						killBlock(dx, i);
				}
			}
				
		}
	}
	
	boolean liveBlockTest(int blockNo, int ip){
		boolean ret = false;
		int dx = (int)Math.pow(2, ip);
		int Lp = L/dx;
		int xp = blockNo%Lp;
		int yp = blockNo/Lp;
		int x = dx*xp;
		int y = dx*yp;
		int site = y*L+x;
		if(aliveLattice[site]) ret = true;
		return ret;
	}
	
	void killBlock(int dx, int block){
		//kill sites in block
		int Lp = L/dx;
		int xp = block%Lp;
		int yp = block/Lp;
		for (int y = yp*dx; y < yp*dx+dx; y++){
			for (int x = xp*dx; x < xp*dx + dx; x++){
				int site = y*L+x; 
				killSite(site);
			}
		}
	}

	void setDeadStrip(int noDeadToPlace){
		for(int i  = 0; i < N; i++){
			if(i < noDeadToPlace) killSite(i);
			else setSite(i);
		}
	}
	
	void setPlaceDeadRandom(int noDeadToPlace){
		while(noDeadSites < noDeadToPlace){
			int randSite = (int)(random.nextDouble()*(double)N);
			if (aliveLattice[randSite]){
				killSite(randSite);
			}
		}
	}
	
	void setRandom(double initPercentDead){
		for(int i = 0; i < N; i++){
			if(random.nextDouble() > initPercentDead){
				setSite(i);
			}else{
				killSite(i);
			}
		}
	}
	
	void setDeadBlock(int noDeadToPlace){
		int sq = (int) Math.sqrt(noDeadToPlace);
		int initialBlock = sq*sq;
		for (int x = 0; x < sq; x++){
			for (int y = 0; y < sq; y++){
				int s = y*L+x;
				killSite(s);
			}
		}
		int extra = noDeadToPlace - initialBlock;
		if(extra > sq){
			for(int i = 0; i < sq; i++) killSite(i*L + sq);
			for (int i = 0; i <= extra-sq; i++) killSite(sq*L + i);
		}else{
			for (int i = 0; i < extra; i++) killSite(i*L + sq);
		}

	}
	
	void placeDeadBlocks(int blockSize, int noDeadBlocksToPlace){
		FileUtil.printlnToFile(infoFileName, "Damage block size", blockSize);
		int noDamageBlocks = N/(blockSize*blockSize);
		boolean [] blockAlive = new boolean [noDamageBlocks];
		for (int i = 0; i < noDamageBlocks; i++) blockAlive[i] = true;
		int noDeadBlocks = 0;
		while(noDeadBlocks < noDeadBlocksToPlace){
			int randBlock = (int)(random.nextDouble()*(double)noDamageBlocks);
			if (blockAlive[randBlock]){
				blockAlive[randBlock]=false;
				noDeadBlocks += 1;
			}
		}

		for (int i = 0; i < N; i++){
			int block = findCG_site(i, blockSize);
			if(blockAlive[block]) setSite(i);
			else killSite(i);
		}
	}
	
	void setDeadBlocks(int blockSize, double percentDeadBlocks){
		FileUtil.printlnToFile(infoFileName, "Damage block size", blockSize);
		int noDamageBlocks = N/(blockSize*blockSize);
		boolean [] blockAlive = new boolean [noDamageBlocks];
		for (int j = 0; j < noDamageBlocks; j++){
			if(random.nextDouble()<percentDeadBlocks){
				blockAlive[j] = false;
			}else{
				blockAlive[j] = true;
			}
		}
		for (int i = 0; i < N; i++){
			int block = findCG_site(i, blockSize);
			if(blockAlive[block]) setSite(i);
			else killSite(i);
		}
	}
	
	void setRandomBlockDamage(int damageBlockSize){
		FileUtil.printlnToFile(infoFileName, "Damage block size", damageBlockSize);
		int noDamageBlocks = N/(damageBlockSize*damageBlockSize);
		double [] damageBlockDamage = new double [noDamageBlocks];
		for (int j = 0; j < noDamageBlocks; j++){
			//assign a random amount of damage
			//half of blocks have no damage
//			damageBlockDamage[j] = random.nextGaussian()*0.07;   // this gives about 3 % damage
//			damageBlockDamage[j] = random.nextGaussian()*0.115;  // this gives about 5 % damage
			damageBlockDamage[j] = random.nextGaussian()*0.25;  // this gives about 10 % damage
//			damageBlockDamage[j] = random.nextGaussian()*0.59;  // this gives about 25 % damage
//			damageBlockDamage[j] = Math.abs(random.nextGaussian());  // this gives about ? % damage
		}
		for (int i = 0; i < N; i++){
			int block = findCG_site(i, damageBlockSize);
			if (random.nextDouble() < damageBlockDamage[block]){
				killSite(i);
			}else{
				setSite(i);
			}
		}
		
	}
	/**
	 * Calculate effective alpha (alpha' or ap) for each lattice site:
	 * phi_i*(1-alpha) = 1-alpha'_i
	 */
	double [] calcAlphaP(){
		double [] ap = new double [noLiveSites];
		int liveSite = 0;
		for (int i = 0; i < N; i ++){
			if(aliveLattice[i]){
				double phi_i = 1.0 - fracDeadNbors[i];
				ap[liveSite] = 1-phi_i*(1-alpha);
				alpha_iHist.accum(ap[liveSite]);
				liveSite += 1;
			}
		}
		double ave = MathTools.mean(ap);
		double var = MathTools.variance(ap);
		double [] ret = new double [2];
		ret[0] = ave;
		ret[1] = var;
		return ret;
	}
	
	void makeNborLists(){
		for (int i = 0; i < N; i++){
			int ct = 0;
			int boundsCt = 0;
			for (int j = 0; j < maxNbors; j++){
				int nborSite = nborList[i][j];
				if(nborSite>=1){
					if (aliveLattice[nborSite]){
						ct += 1;
					}
				}else{
					boundsCt += 1;
				}
				int noDead = maxNbors - ct;
				fracDeadNbors[i] = (double)noDead/maxNbors;
				if(deadDissipation==true) noNborsForSite[i] = noNbors; 
				else if (deadDissipation==false){
					if(boundaryConditions=="Open")
						noNborsForSite[i]=ct+boundsCt;
					else if(boundaryConditions == "Periodic")
						noNborsForSite[i]=ct;
				}
				
			}
		}
//		if(deadDissipation==true){
//			for (int i = 0; i < N; i++)	noNborsForSite[i] = noNbors; 
//			for (int i = 0; i < N; i++){
//				int ct = 0;
//				for (int j = 0; j < maxNbors; j++){
//					if (aliveLattice[nborList[i][j]]){
//						ct += 1;
//					}
//					int noDead = maxNbors - ct;
//					fracDeadNbors[i] = (double)noDead/maxNbors;
//				}
//			}
//		}else if (deadDissipation==false){
//			for (int i = 0; i < N; i++){
//				int ct = 0;
//				for (int j = 0; j < maxNbors; j++){
//					if (aliveLattice[nborList[i][j]]){
//						ct += 1;
//					}
//					int noDead = maxNbors - ct;
//					fracDeadNbors[i] = (double)noDead/maxNbors;
//					noNborsForSite[i]=ct;
//				}
//			}
//
//		}
	}
	
	int findNoCircleNbors(){
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
	
	void findCircleNbors(){
		if(boundaryConditions == "Periodic"){
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
		}else if(boundaryConditions == "Open"){
			int nborIndex = 0;
			for (int s = 0; s < N; s++){
				nborIndex = 0;
				int x = s%L;
				int y = s/L;
				for(int dy = -R; dy <= R; dy++){
					for(int dx = -R; dx <= R; dx++){
						double distance = Math.sqrt(dx*dx + dy*dy);
						if (distance <= R){
							int xx = (x+dx);
							int yy = (y+dy);
							if(xx>=0 & xx<L & yy>=0 & yy<L){
								int nborSite = yy*L+xx;
								nborList[s][nborIndex] = nborSite;
							}else{
								nborList[s][nborIndex] = -1;
							}
							nborIndex += 1;
						}
					}
				}
			}
		}
	}
	
	public void initEquilibrate(int maxPlateUpdates){
			plateUpdates = -maxPlateUpdates;
	}
	
	public void equilibrate(){

		bringToFailure();
		eqFail(epicenterSite);  //Just distributes stress to each neighbor
		int nextSiteToFail = checkFail();  //each checkFail loops through lattice once.
		while (nextSiteToFail >= 0){
			eqFail(nextSiteToFail);
			nextSiteToFail = checkFail();
			avSize += 1;
		}
		plateUpdates += 1;
	}
	
	void bringToFailure(){
		epicenterSite = siteWithMaxStress();	
		double stressAdd = tStress - stress[epicenterSite];
		if(stressAdd < 0) System.out.println("Error: stress already above failure");
		for (int i = 0; i < N; i++){
			if(aliveLattice[i]) stress[i] += stressAdd;
		}
		avSize = 1;
	}
	
	public void ofcStep(){
		bringToFailure();
		dt=1;
	
		fail(epicenterSite);
		int nextSiteToFail = checkFail();
		while(nextSiteToFail >= 0){
			fail(nextSiteToFail);
			nextSiteToFail = checkFail();
			avSize += 1;
			if (avSize%N==0){
				System.out.println("Error:  Av Size = " + avSize);
				nextSiteToFail = -1;
				FileUtil.printlnToFile(infoFileName, "Av Size = " , avSize);
			}
		}

		epicenterSizeAccum();
		cg_time += dt;
		plateUpdates += 1;
	}
	
	void epicenterSizeAccum(){
			epicenterCount[epicenterSite] +=1;
			epicenterSize[epicenterSite] += avSize;		
	}
	
	public void clearCounts(){
		for (int i = 0; i < N; i++){
			epicenterCount[i] = 0;
			epicenterSize[i] = 0;
		}
	}
		
	int checkFail(){
		int s = siteWithMaxStress();
		if (stress[s] < tStress) s = -1; 
		return s;
	}

	void eqFail(int s){
		double resStressWithNoise = calcResNoise(s);
		double stressDrop = stress[s] - resStressWithNoise;
		double stressPerNbor = ((1-alpha)*(stressDrop))/(double)noNborsForSite[s];
		stress[s] = resStressWithNoise;
		//distribute stress to the alive nbors
		int checkAliveCount = 0;
		for (int i = 0; i < maxNbors; i++){
			int nborSite = nborList[s][i];
			if (nborSite >= 0){
				if(aliveLattice[nborSite]){
					stress[nborSite] += stressPerNbor;
					checkAliveCount += 1;
				}
			}
		}
	}
	
	void fail(int s){
		double resStressWithNoise = calcResNoise(s);
		double stressDrop = stress[s] - resStressWithNoise;
		double stressPerNbor = ((1-alpha)*(stressDrop))/(double)noNborsForSite[s];
//		if(damage == "Cascade")
			if(deadDissipation) calcDissipation(stressDrop, s);
		stress[s] = resStressWithNoise;
		//distribute stress to the alive nbors
		int checkAliveCount = 0;
		for (int i = 0; i < maxNbors; i++){
			int nborSite = nborList[s][i];
			if(nborSite >= 0){
				if(aliveLattice[nborSite] ){
					stress[nborSite] += stressPerNbor;
					checkAliveCount += 1;
				}
			}
		}
	}
	

	void calcDissipation(double stressDrop, int site){
		rateCt += 1.0;
		//calc alpha dissipation = alpha x % alive
		double fracAlive = 1.0 - fracDeadNbors[site];
		alphaDiss = alpha*fracAlive;
		alphaDissAccum += alphaDiss;
		alphaDissRate = alphaDissAccum/rateCt;
//		System.out.println("aDiss = " + alphaDiss + " cgTime = " + cg_time + " rate = " + alphaDissRate);
		//calc deadDiss = fraction dead
		deadDiss = fracDeadNbors[site];
		deadDissAccum += deadDiss;
		deadDissRate = deadDissAccum/rateCt;
	}
	
	void p(String s){
		System.out.println(s);
	}

}
