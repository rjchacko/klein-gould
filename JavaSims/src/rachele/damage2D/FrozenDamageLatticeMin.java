package rachele.damage2D;

import rachele.damage2D.AlphaDistribution;
import rachele.damage2D.multidx.Damage;
import rachele.util.FileUtil;
import rachele.util.MathTools;
import scikit.jobs.params.Parameters;

public class FrozenDamageLatticeMin extends AbstractFrozenDamage{

		public boolean [] aliveLattice;  // No of failures at each site: once a site is dead -> = -1;
		public int [] failedSitesList;
		public String boundaryConditions;
		public double aveGamma;
		public double [] gamma;
		int [] nearestNbors = new int [4]; 
		public int noLiveSites; 
		int noDeadSites;
		int nextSiteToFail;
		int maxNbors;
		int [] noNborsForSite;  // No of nbors at each site (for when different for each site.)
		public int [][]nborList;  // List of neighbors for each site: nborList[site = 0..N-1][nbor index = 0..maxNbors]  This is gonna be really big
		double nborDistanceSum;
		double [] fracDeadNbors;
		String damage;
		String infoFileName;
		boolean deadDissipation;
		public Damage dl;
		public AlphaDistribution ad;
		
		public FrozenDamageLatticeMin(Parameters params, String name) {
			infoFileName = name;
			FileUtil.initFile(infoFileName, params);
			FileUtil.printlnToFile(infoFileName, "# Runlog for OFC_DamageLattice");
			setParameters(params);
			allocate();
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
			noNbors = findNoCircleNbors(R);
			random.setSeed(params.iget("Random Seed"));	
			tStress = 1.0;
			rStress = params.fget("Residual Stress");
			maxResNoise = params.fget("Res. Max Noise");
			damage = params.sget("Type of Damage");
			boundaryConditions = params.sget("Boundary Conditions");
		}
		
		void allocate(){
			stress = new double [N];
			noNborsForSite = new int [N];
			fracDeadNbors = new double [N];
			alpha = new double [N];
			failedSitesList = new int [N];
			for(int i = 0; i<N; i++)
				failedSitesList[i] = -1;
		}
		
		public void initLattice(Parameters params){
			
			p("Initializing lattice");
			for (int i = 0; i < N; i++){
				stress[i] = random.nextFloat()*(tStress-rStress)+rStress;	
			}
			
			maxNbors = findCircleNbors(R);
			nborList = new int [N][maxNbors];
			p("Finding circle nbors");
			findCircleNbors(R, maxNbors);	
			p("Setting damage");
			dl = new Damage(pow, R, params.iget("Random Seed"), infoFileName);	
			boolean [] dAlive = Damage.setDamage(damage, params.iget("Dead Parameter"), params.fget("Init Percent Dead"), params.iget("Number Dead"));
			p("Damage Set");
			aliveLattice = new boolean [N];
			int ct = 0;
			for (int i =0; i < N; i++){
				aliveLattice[i]=dAlive[i];
				if(aliveLattice[i]) ct += 1;
			}
			noLiveSites = ct;	
			noDeadSites = N-ct;
			System.out.println("no dead = " + noDeadSites);

		for(int i = 0; i < N; i++) {
				if(aliveLattice[i]==false) stress[i] = 0;
			}
			
			p("Making neighbor lists...");
			makeNborLists();
			
			p("Setting gamma array");
			ad = new AlphaDistribution(pow, R, params.iget("Random Seed"), infoFileName);
			alpha = AlphaDistribution.setAlphaArray(params.sget("Alpha Distribution"), params.fget("Dissipation Param"));
			
			p("Calculating gamma and variance");
			double [] ap = calcGamma();
			aveGamma = ap[0];
			System.out.println("Gamma = " + ap[0] + ", var = " + ap[1]);// +  " Connection Param (dead Excluded) = " + ap[2] + " var = " + ap[3]+ " Conection Param = " + ap[4] + " var = " + ap[5]);
			FileUtil.printlnToFile(infoFileName, "gamma = " , ap[0], " var = ", ap[1]);
			double percentSitesDamaged = (double)noDeadSites/(double)N;
			FileUtil.printlnToFile(infoFileName, "# No of sites damaged = " , noDeadSites);
			FileUtil.printlnToFile(infoFileName, "# Percent of sites damaged = ", percentSitesDamaged);
			System.out.println("Percent of sites damaged = "+ percentSitesDamaged);
			params.set("Percent Damage", percentSitesDamaged);

			cg_time = 0;
			plateUpdates = 0;

			p("Lattice init complete.");
		}

		/**
		 * Calculate effective alpha (alpha' or ap) for each lattice site:
		 * phi_i*(1-alpha) = 1-alpha'_i
		 */
		double [] calcGamma(){
			gamma = new double [N];
			double [] gammaParam = new double [noLiveSites];
			int liveSite = 0;
			for (int i = 0; i < N; i ++){
				gamma[i] = 1.0;
				if(aliveLattice[i]){
					double phi_i = 1.0 - fracDeadNbors[i];
					gammaParam[liveSite] = 1-phi_i*(1-alpha[i]);
					gamma[i] = gammaParam[liveSite];
					liveSite += 1;
				}
			}
			
			double ave = MathTools.mean(gammaParam);
			double var = MathTools.variance(gammaParam);
			double [] ret = new double [2];
			ret[0] = ave;
			ret[1] = var;
			return ret;
		}
			
		
		/**
		 * Calculates connection parameter for a given range, maxrange
		 * Warning- this method reassigns neighbors.  Should not be called during a typical run,
		 * only by separate apps for calculating connection params.
		 */

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
		
		int findNoCircleNbors(int range){
			int count = 0;
			nborDistanceSum = 0;
			 for(int dy = -range; dy <= range; dy++){
				 for(int dx = -range; dx <= range; dx++){
					 double distance = Math.sqrt(dx*dx + dy*dy);
					 if (distance <= range){
							 count += 1;
					 }
				 }
			 }
			 count = count -1;
			 return count;
		}
		
		double [] findNborDistances(int range, int maxNo){
			int count = 0;
			int nearestNborIndex = 0;
			double [] ret = new double [maxNo];
			 for(int dy = -range; dy <= range; dy++){
				 for(int dx = -range; dx <= range; dx++){
					 double distance = Math.sqrt(dx*dx + dy*dy);
					 if (distance <= range){
						 if(distance == 1.0){
							 nearestNbors[nearestNborIndex]=count;
							 System.out.println("nearest neighbor = " + nearestNbors[nearestNborIndex] + " distance = " + distance);
							 nearestNborIndex+=1;
						 }
						 if(distance != 0){
							 ret[count] = distance/(double)R;
							 count += 1;
						 }else{
							 System.out.println("dist = 0");
						 }
					 }
				 }
			 }
			 return ret;
		}
		
		void findCircleNbors(int range, int maxNo){
			
			if(boundaryConditions == "Periodic"){
				int nborIndex = 0;
				for (int s = 0; s < N; s++){
					nborIndex = 0;
					int x = s%L;
					int y = s/L;
					for(int dy = -range; dy <= range; dy++){
						for(int dx = -range; dx <= range; dx++){
							double distance = Math.sqrt(dx*dx + dy*dy);
							if (distance <= range){
								if(distance !=0){
									int xx = (x+dx+L)%L;
									int yy = (y+dy+L)%L;
									int nborSite = yy*L+xx;
									nborList[s][nborIndex] = nborSite;
									nborIndex += 1;
								}
							}
						}
					}
				}
			}else if(boundaryConditions == "Open"){
				int nborIndex = 0;
				for (int s = 0; s < N; s++){
					nborIndex = 0;
					int x = s%L;
					int y = s/L;
					for(int dy = -range; dy <= range; dy++){
						for(int dx = -range; dx <= range; dx++){
							double distance = Math.sqrt(dx*dx + dy*dy);
							if (distance <= range){
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
		}
		
		public void ofcStep(){
			avSize=0;
			bringToFailure();
			failedSitesList[avSize]=epicenterSite;
			avSize = 1;
			dt=1;
			fail(epicenterSite);
			int nextSiteToFail = checkFail();
			while(nextSiteToFail >= 0){
				failedSitesList[avSize]=nextSiteToFail;
				fail(nextSiteToFail);
				nextSiteToFail = checkFail();
				avSize += 1;
				if (avSize%N==0){
					System.out.println("Error:  Av Size = " + avSize);
					nextSiteToFail = -1;
					FileUtil.printlnToFile(infoFileName, "Av Size = " , avSize);
				}
			}
			cg_time += dt;
			plateUpdates += 1;
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

		
		void p(String s){
			System.out.println(s);
		}

	}

