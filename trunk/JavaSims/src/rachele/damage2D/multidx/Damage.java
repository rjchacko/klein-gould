package rachele.damage2D.multidx;

import kip.util.Random;
import rachele.util.FileUtil;
import rachele.util.MathTools;

public class Damage {

	static int pow;
	static int L;
	static int N;
	static int R;
//	static int noDeadSites;
//	static int noLiveSites;
	static boolean [] aliveLattice;
	static boolean [] setLattice;
	static Random random = new Random();
	static String infoFileName;
	int [] noNborsForSite;
	int [][]nborList;
	double [] fracDeadNbors;
	
	public Damage(int power, int range, int randomSeed, String fileName) {
		pow=power;
		L=(int) Math.pow(2, pow);
		N=L*L;
		System.out.println("N = " + N);
		R=range;
		aliveLattice = new boolean [N];
		setLattice = new boolean [N];
		setRandomSeed(randomSeed);
		infoFileName=fileName;
//		noLiveSites = N;
//		noDeadSites = 0;
	}
	
	static public boolean [] setDamage(String damageType, int deadParam, double initPercentDead, int noDead){
		for(int i = 0; i < N; i ++) aliveLattice[i] = true;
		if(damageType=="Random") setRandom(initPercentDead);
		else if(damageType == "PlaceRandomDead")	setPlaceDeadRandom(noDead);
		else if(damageType=="DeadStrip") setDeadStrip(noDead);
		else if(damageType=="RandomBlocks") setRandomBlockDamage(deadParam);
		else if(damageType == "DeadBlock") setDeadBlock(noDead);
		else if(damageType == "Cascade") setCascadeDamage(initPercentDead);
		else if(damageType == "CascadeRandom") setCascadeRandom(initPercentDead, deadParam);
		else if(damageType == "DeadRectangle") setDeadRectangle(noDead);
		else if(damageType == "DeadBlocks") setDeadBlocks(deadParam, initPercentDead);
		else if(damageType == "PlaceDeadBlocks") placeDeadBlocks(deadParam, noDead);
		else if(damageType == "AliveCascade") setCascadeLiveBlocks(initPercentDead);
		else System.out.println("Error!");
		System.out.println("site 0 = " + aliveLattice[0]);
		return aliveLattice;
	}
	
	static void setRandomSeed(int rs){
		random.setSeed(rs);	
	}
	
	
	static void killSite(int site){
//		noDeadSites+=1;
//		noLiveSites-=1;
		aliveLattice[site] = false;
	}
	
	static void setSite(int site){
		aliveLattice[site] = true;
	}
	
	private static void setDeadRectangle(int noDead) {
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

	private static void setCascadeRandom(double p, int d) {
		System.out.println("Setting Cascade Random Damage");
		int maxPow = pow-2;
		for (int i = 0; i < N; i++) setLattice[i] = false;

		double blockP=.25;
		for (int ip = maxPow; ip >= 0; ip--){
			int dx = (int)Math.pow(2 ,ip);
			int Lp = L/dx;
			int Np = Lp*Lp;
			for (int i = 0; i < Np; i++){
				// set blocks with probability if they are not already set
				if(setBlockTest(i, ip)==false){
					if (random.nextDouble()<blockP) 
//						setBlock(dx, i);
						setBlockWithp(dx, i, p, d);
				}
			}
				
		}
		
	}
	
	static boolean setBlockTest(int blockNo, int ip){
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

	static void setBlockWithp(int dx, int block, double p, int d){
		//all sites in block
		int Lp = L/dx;
		int xp = block%Lp;
		int yp = block/Lp;
		
//		double pr =random.nextDouble()+p;
		double pr =(100*random.nextGaussian()/(double)d+p);
		
		for (int y = yp*dx; y < yp*dx+dx; y++){
			for (int x = xp*dx; x < xp*dx + dx; x++){
				int site = y*L+x;
				if(random.nextDouble()<pr) killSite(site);
				setLattice[site] = true;
			}
		}
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

	static void setCascadeLiveBlocks(double p){
//		System.out.println("Setting Cascading Damage");
		int maxPow = pow-2;
		//largest lattice size is 8 x 8
		for (int i = 0; i < N; i++) aliveLattice[i] = false;
		for (int ip = maxPow; ip >= 0; ip--){
			int dx = (int)Math.pow(2 ,ip);
			int Lp = L/dx;
			int Np = Lp*Lp;
			System.out.println("dx = " + dx + " Lp = " + Lp);	
			for (int i = 0; i < Np; i++){
				double pr = p*Math.pow(Math.E,(-Math.pow(dx-R,2)/100));
				// kill alive blocks with probability p
				if(liveBlockTest(i, ip)==false){
					if (random.nextDouble()<pr) 
						setLiveBlock(dx, i);
				}
			}
				
		}
	}
	
	static void setCascadeDamage(double p){
//		System.out.println("Setting Cascading Damage");
		int maxPow = pow-2;
		//largest lattice size is 8 x 8

		for (int ip = maxPow; ip >= 0; ip--){
			int dx = (int)Math.pow(2 ,ip);
			int Lp = L/dx;
			int Np = Lp*Lp;
			System.out.println("dx = " + dx + " Lp = " + Lp);	
			for (int i = 0; i < Np; i++){
				double pr = p*Math.pow(Math.E,(-Math.pow(dx-R,2)/100));
				// kill alive blocks with probability p
				if(liveBlockTest(i, ip)){
					if (random.nextDouble()<pr) 
						killBlock(dx, i);
				}
			}
				
		}
	}
	
	static boolean liveBlockTest(int blockNo, int ip){
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
	
	static void setLiveBlock(int dx, int block){
		//kill sites in block
		int Lp = L/dx;
		int xp = block%Lp;
		int yp = block/Lp;
		for (int y = yp*dx; y < yp*dx+dx; y++){
			for (int x = xp*dx; x < xp*dx + dx; x++){
				int site = y*L+x; 
				setSite(site);
			}
		}
	}
	
	static void killBlock(int dx, int block){
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

	static void setDeadStrip(int noDeadToPlace){
		for(int i  = 0; i < N; i++){
			if(i < noDeadToPlace) killSite(i);
			else setSite(i);
		}
	}
	
	static void setPlaceDeadRandom(int noDeadToPlace){
		int noDead = 0;
		while(noDead < noDeadToPlace){
			int randSite = (int)(random.nextDouble()*(double)N);
			if (aliveLattice[randSite]){
				killSite(randSite);
				noDead+=1;
			}
		}
	}
	
	static void setRandom(double initPercentDead){
		for(int i = 0; i < N; i++){
			if(random.nextDouble() > initPercentDead){
				setSite(i);
			}else{
				killSite(i);
			}
		}
	}
	
	static void setDeadBlock(int noDeadToPlace){
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
	
	static void placeDeadBlocks(int blockSize, int noDeadBlocksToPlace){
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
	
	static void setDeadBlocks(int blockSize, double percentDeadBlocks){
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
	
	static void setRandomBlockDamage(int damageBlockSize){
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

	static int findCG_site(int s, int blockSize){
		int x = s % L;
		int y = s / L;
		int xp = x / blockSize;
		int yp = y / blockSize;
		int cgSite = yp *(L/blockSize) + xp;
		return cgSite;
	}
	
	/**
	 * Calculate effective alpha (alpha' or ap) for each lattice site:
	 * phi_i*(1-alpha) = 1-alpha'_i
	 */
	double [] calcAlphaP(double alpha, int R){
		noNborsForSite = new int [N];
		fracDeadNbors = new double [N];
		int maxNbors = findNoCircleNbors();	
		nborList = new int [N][maxNbors];
		findCircleNbors("Periodic", maxNbors);
		makeNborLists(maxNbors, "Periodic", true);
		int ct=0;
		for (int i = 0; i < N; i ++){
			if(aliveLattice[i]) ct+=1;
		}
		double [] ap = new double [ct];
		int liveSite = 0;
		for (int i = 0; i < N; i ++){
			if(aliveLattice[i]){
				double phi_i = 1.0 - fracDeadNbors[i];
				ap[liveSite] = 1-phi_i*(1-alpha);
//				alpha_iHist.accum(ap[liveSite]);
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
	
	void makeNborLists(int maxNbors, String boundaryConditions, boolean deadDissipation){
		int noNbors = findCircleNbors(R);
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
	
	void findCircleNbors(String boundaryConditions, int maxNbors){
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
