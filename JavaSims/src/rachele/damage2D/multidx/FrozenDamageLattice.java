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
	public Histogram alpha_iHist = new Histogram(0.01);
	public double alphaPrime;
	public double [] alphaP;
	
	public int noLiveSites; 
	
	int noDeadSites;
	int nextSiteToFail;
	int maxNbors;
	int [] noNborsForSite;  // No of nbors at each site (for when different for each site.)
	int [][]nborList;  // List of neighbors for each site: nborList[site = 0..N-1][nbor index = 0..maxNbors]  This is gonna be really big
	
	double [] fracDeadNbors;
	
	double alphaDissAccum;
	double deadDissAccum;
	double rateCt;
	boolean deadDissipation;

	String damage;
	String infoFileName;
	
	public Damage dl;
	
	public FrozenDamageLattice(Parameters params, String name) {
		infoFileName = name;
		FileUtil.initFile(infoFileName, params);
		FileUtil.printlnToFile(infoFileName, "# Runlog for OFC_DamageLattice");
		setParameters(params);
		allocate();
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
//		double alphaParam = params.fget("Dissipation Param");   // = alpha
		rStress = params.fget("Residual Stress");
		maxResNoise = params.fget("Res. Max Noise");
		damage = params.sget("Type of Damage");
		deadDissRate = 0.0;
		alphaDissRate = 0.0;
		rateCt = 0;
		boundaryConditions = params.sget("Boundary Conditions");
	}
	
	void allocate(){
		stress = new double [N];
		epicenterCount = new int [N];
		epicenterSize = new int [N];
//		aliveLattice = new boolean [N];
//		setLattice = new boolean [N];
		noNborsForSite = new int [N];
		fracDeadNbors = new double [N];
		alpha = new double [N];
	}
	
	void initLattice(Parameters params){
		
		p("Initializing lattice");
		maxNbors = findNoCircleNbors();	
		System.out.println("max no nbors = " + maxNbors);
		FileUtil.printlnToFile(infoFileName, "# max no nbors = " , maxNbors);
//		p("Allocating nbor list...");
		nborList = new int [N][maxNbors];	
//		p("Finding neighbors...");
		findCircleNbors();

		
//		p("Finding stress...");
		for (int i = 0; i < N; i++){
			stress[i] = random.nextFloat()*(tStress-rStress)+rStress;	
		}
		
//		p("Setting damage...");
		dl = new Damage(pow, R, params.iget("Random Seed"), infoFileName);
		
		boolean [] dAlive = Damage.setDamage(damage, params.iget("Dead Parameter"), params.fget("Init Percent Dead"), params.iget("Number Dead"));
		aliveLattice = new boolean [N];
		int ct = 0;
		for (int i =0; i < N; i++){
			aliveLattice[i]=dAlive[i];
			if(aliveLattice[i]) ct += 1;
		}
		noLiveSites = ct;	
		noDeadSites = N-ct;

		for(int i = 0; i < N; i++) {
			if(aliveLattice[i]==false) stress[i] = 0;
		}
		
		p("Making neighbor lists...");
		makeNborLists();
		
		p("Setting alpha array");
		setAlphaArray(params.sget("Alpha Distribution"), params.fget("Dissipation Param"));
		
		p("Calculating alpha prime and variance");
		double [] ap = calcAlphaP();
		alphaPrime = ap[0];
		System.out.println("Alpha prime = " + ap[0] + ", var = " + ap[1]);
		FileUtil.printlnToFile(infoFileName, "Alpha' = " , ap[0], " var = " , ap[1]);
		double percentSitesDamaged = (double)noDeadSites/(double)N;
		FileUtil.printlnToFile(infoFileName, "# No of sites damaged = " , noDeadSites);
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
	
	void setAlphaArray(String alphaDistribution, double alphaParam){
		if(alphaDistribution=="Flat Random"){
			for (int i = 0; i < N; i++) alpha[i]=random.nextDouble();
		}else if(alphaDistribution=="Gaussian"){
			for (int i = 0; i < N; i++){
				double a =random.nextGaussian()/10 + alphaParam;
				if(a > 1.0) alpha[i] = 1.0;
				else if (a < 0) alpha[i] = 0.0;
				else alpha[i] = a;
			}
		}else if(alphaDistribution=="Constant"){
			for (int i = 0; i < N; i++) alpha[i]=alphaParam;
		}else if(alphaDistribution=="Gaussian about zero"){
			for (int i = 0; i < N; i++){
				double a =Math.abs(random.nextGaussian()/5.0);
				while(Math.abs(a)>1.0){
					a =(random.nextGaussian()/5.0);
				}
				alpha[i] = Math.abs(a);
			}
		}else if(alphaDistribution=="Gaussian split"){
			for (int i = 0; i < N; i++){
				double a =(random.nextGaussian()/10.0);
				while(Math.abs(a)>1.0){
					a =(random.nextGaussian()/10.0);
				}
				if(a < 0.0) a = 1.0 + a;
				alpha[i] = a;
			}
		}else if(alphaDistribution=="Gaussian about half"){
			for (int i = 0; i < N; i++){
				double a =(random.nextGaussian()/10.0);
				while(Math.abs(a)>0.5){
					a =(random.nextGaussian()/10.0);
				}
				a += 0.5;
				alpha[i] = a;
			}
		}else if(alphaDistribution=="Dead Blocks"){
			boolean [] aa = Damage.setDamage("Place Dead Blocks", 32,0.0, 32);
			double a;
			for(int i = 0; i < N; i++){
				a=(random.nextGaussian()/10.0);
				while(Math.abs(a)>1.0){
					a =(random.nextGaussian()/10.0);
				}
				if(aa[i]){
					alpha[i] = 1.0- Math.abs(a);
				}else{
					alpha[i] = Math.abs(a);
				}
			}	
		}else if(alphaDistribution=="Many Gaussians"){
			int dx = 128;
			int Lp = L/dx;
			int Np = Lp*Lp;
			boolean [] blockSet = new boolean [Np];
			for(int i = 0; i < Np; i ++) blockSet[i] = false;
			for(int i = 0; i < Np; i ++){
				//choose a block
				int block = (int)(random.nextDouble()*(double)(Np));
				while(blockSet[block]){
					block = (int)(random.nextDouble()*(double)(Np));
				}
				double center = ((double)i+0.5)/(double)(Np);
				setGaussianBlock(center, block, dx, Np);
				blockSet[block] = true;
			}
			
		}else if(alphaDistribution=="Fractal"){
			for (int i = 0; i < N; i++)
				alpha[i]=-1.0;
			
			int maxPow = pow-2;  //largest lattice size is 8 x 8
			int [] intervalCount = new int [maxPow+1];
			for (int ip = maxPow; ip >= 0; ip--){
				int dx = (int)Math.pow(2 ,ip);
				int Lp = L/dx;
				int Np = Lp*Lp;
				double center = 1.0-((double)ip+0.5)/(double)(maxPow+1);
				for (int i = 0; i < Np; i++){
					double pr = 0.25;  						//Set one quarter of blocks to alpha distributed about center
					if(alpha[getFirstBlockSite(i, ip)]<0){
						if (random.nextDouble()<pr){
							setGaussianBlock(center, i, dx, maxPow+1);
							intervalCount[ip]+=dx*dx;
						}
					}
				}	
			}
			int [] unsetList = randomizeUnsetSites();
			for (int i = 0; i < unsetList.length; i++){
				int minI = findMinInterval(intervalCount);
				double center = 1.0-((double)minI+0.5)/(double)(maxPow+1);
				setGaussianBlock(center, unsetList[i], 1, maxPow+1);
				intervalCount[minI]+=1;
			}
		}else if (alphaDistribution=="Quarters"){
			for (int i = 0; i < N; i++)
				alpha[i]=-1.0;
			
			//set one block 1/4 the lattice size
			int ip = pow - 1;
			int dx = (int)Math.pow(2 ,ip);
			double center = 0.5/4.0;
			setGaussianBlock(center, 1, dx, 4);
			//set 4 blocks 1/8 the size of lattice
			ip = pow - 2;
			dx = (int)Math.pow(2 ,ip);
			int Lp = L/dx;
			int Np = Lp*Lp;
			center = 1.5/4.0;
			int noSet = 0;
			while(noSet < 4){
				//pick a block of size dx
				int randBlock = (int)(random.nextDouble()*Np);
				if(alpha[getFirstBlockSite(randBlock, ip)]<0){
					setGaussianBlock(center, randBlock, dx, 4);
					noSet+=1;
					System.out.println("Set block of size " + dx);
				}
			}
			//set 16
			ip = pow - 3;
			dx = (int)Math.pow(2 ,ip);
			Lp = L/dx;
			Np = Lp*Lp;
			center = 2.5/4.0;
			noSet = 0;
			while(noSet < 16){
				//pick a block of size dx
				int randBlock = (int)(random.nextDouble()*Np);
				if(alpha[getFirstBlockSite(randBlock, ip)]<0){
					setGaussianBlock(center, randBlock, dx, 4);
					noSet+=1;
					System.out.println("Set block of size " + dx);
				}
			}
			//set 16 more blocks
			center = 3.5/4.0;
			noSet = 0;
			while(noSet < 16){
				//pick a block of size dx
				int randBlock = (int)(random.nextDouble()*Np);
				if(alpha[getFirstBlockSite(randBlock, ip)]<0){
					setGaussianBlock(center, randBlock, dx, 4);
					noSet+=1;
					System.out.println("Set block of size " + dx);
				}
			}
			int [] unsetList = randomizeUnsetSites();
			System.out.println("no unset = " + unsetList.length);
		}else if (alphaDistribution=="Eights"){
			for (int i = 0; i < N; i++)
				alpha[i]=-1.0;
			
			int maxPow = pow-2;  //largest lattice size is 4 X 4
			for (int ip = maxPow; ip > 0; ip--){
				int dx = (int)Math.pow(2 ,ip);
				int Lp = L/dx;
				int Np = Lp*Lp;
				double center = (7.5-(double)ip)/8.0;
					int noSet = 0;
					int maxToSet = Np/8;
					while(noSet < maxToSet){
						//pick a block of size dx
						int randBlock = (int)(random.nextDouble()*Np);
						if(alpha[getFirstBlockSite(randBlock, ip)]<0){
							setFlatBlock(center, randBlock, dx, 8);
							noSet+=1;
							System.out.println("Set block of size " + dx);
						}
					}

			}
			double center = 7.5/8.0;
			for (int i = 0; i < N; i++){
				if(alpha[i]< 0){
					setFlatBlock(center, i, 1, 8);
				}
			}
		}
	}
	
	int [] randomizeUnsetSites(){
		int ct = 0;
		for (int i = 0; i < N; i++){
			if(alpha[i]<0){
				ct += 1;
			}
		}
		int [] list = new int [ct];
		for (int i = 0; i < ct; i++){
			list[i] = -1;
		}
		for (int i = 0; i < N; i++){
			if(alpha[i]<0){
				int index = (int)(random.nextDouble()*ct);
				while(list[index]>0){
					index = (int)(random.nextDouble()*ct);
				}
				list[index]=i;
			}
		}
		return list;
	}
	
	int findMinInterval(int [] a){
		int minCt = N;
		int minI = -1;
		for (int i = 0; i < a.length; i++){
			if (a[i] < minCt){
				minCt = a[i];
				minI = i;
			}
		}
		return minI;
	}
	
	void setGaussianBlock(double center, int block, int dx, int noInts){
		int Lp = L/dx;
		int xp = block%Lp;
		int yp = block/Lp;
		for (int y = yp*dx; y < yp*dx+dx; y++){
			for (int x = xp*dx; x < xp*dx + dx; x++){
				int site = y*L+x; 
				double rand = random.nextGaussian();
				while(Math.abs(rand)>=1.0){
					rand = random.nextGaussian();
				}
				alpha[site] = rand/(double)(noInts*2)+center;
			}
		}
	}
	
	void setFlatBlock(double center, int block, int dx, int noInts){
		int Lp = L/dx;
		int xp = block%Lp;
		int yp = block/Lp;
		double intervalWidth = 1.0/noInts;
		for (int y = yp*dx; y < yp*dx+dx; y++){
			for (int x = xp*dx; x < xp*dx + dx; x++){
				int site = y*L+x; 
				double rand = random.nextDouble()*intervalWidth;

				alpha[site] = rand+center-intervalWidth/2.0;
			}
		}
	}
	
	int getFirstBlockSite(int blockNo, int ip){
		int dx = (int)Math.pow(2, ip);
		int Lp = L/dx;
		int xp = blockNo%Lp;
		int yp = blockNo/Lp;
		int x = dx*xp;
		int y = dx*yp;
		int site = y*L+x;
		return site;
	}
	
	/**
	 * Calculate effective alpha (alpha' or ap) for each lattice site:
	 * phi_i*(1-alpha) = 1-alpha'_i
	 */
	double [] calcAlphaP(){
		alphaP = new double [N];
		double [] ap = new double [noLiveSites];
		int liveSite = 0;
		for (int i = 0; i < N; i ++){
			alphaP[i] = 1.0;
			if(aliveLattice[i]){
				double phi_i = 1.0 - fracDeadNbors[i];
				ap[liveSite] = 1-phi_i*(1-alpha[i]);
				alpha_iHist.accum(ap[liveSite]);
				alphaP[i] = ap[liveSite];
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
		if(deadDissipation==true){
			for (int i = 0; i < N; i++)	noNborsForSite[i] = noNbors; 
			for (int i = 0; i < N; i++){
				int ct = 0;
				for (int j = 0; j < maxNbors; j++){
					if (aliveLattice[nborList[i][j]]){
						ct += 1;
					}
					int noDead = maxNbors - ct;
					fracDeadNbors[i] = (double)noDead/maxNbors;
				}
			}
		}else if (deadDissipation==false){
			for (int i = 0; i < N; i++){
				int ct = 0;
				for (int j = 0; j < maxNbors; j++){
					if (aliveLattice[nborList[i][j]]){
						ct += 1;
					}
					int noDead = maxNbors - ct;
					fracDeadNbors[i] = (double)noDead/maxNbors;
					noNborsForSite[i]=ct;
				}
			}

		}
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
		double stressPerNbor = ((1-alpha[s])*(stressDrop))/(double)noNborsForSite[s];
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
		double stressPerNbor = ((1-alpha[s])*(stressDrop))/(double)noNborsForSite[s];
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
		alphaDiss = alpha[site]*fracAlive;
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