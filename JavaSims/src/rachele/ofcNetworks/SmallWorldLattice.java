package rachele.ofcNetworks;

import scikit.jobs.params.Parameters;

public class SmallWorldLattice extends AbstractOFC_Lattice{

	
	double prob;
	int randomOrder [];
	
	/**
	 * Constructor for Small World Lattice Class
	 */
	public SmallWorldLattice(Parameters params){
		
		setParams(params);
		initStress();
		setNearNbors();
		rewire();
	}
	
	void setParams(Parameters params){
		L = params.iget("L");
		N = L*L;
		prob = params.fget("Rewire Probability");
		System.out.println("Rewire Probability = " + prob);
		stress = new double [N];
		alpha = params.fget("alpha");
		nbor = new int [N][4];			// 4 neighbors for nearest neighbor lattice
		noNbors = 4;
		sitesToFail = new int [N];
		stressAdd = new double [N];
		randomOrder = new int [N];
	}
	
	void rewire(){
		randomizeOrder();
		for (int site1 = 0; site1 < N; site1++){
			double r = random.nextDouble();
			if (r < prob){												//test for rewire
				int site2 = -1;	
				int nbor1Index = -1;
				int nbor2Index = -1;
				boolean notReady = false;
				do{
					nbor1Index = (int)(random.nextDouble()*noNbors); 
					int nbor1 = nbor[site1][nbor1Index];
					
					site2 = (int)(random.nextDouble()*N);  //choose a random site
					if (site1 == site2){						//Make sure you haven't picked the 
						notReady = true;
					} else{
						if(connected(site1, site2)){
							notReady = true;
						}else{
							nbor2Index = (int)(random.nextDouble()*noNbors);			
							int nbor2 = nbor[site2][nbor2Index];
							if(nbor1 == nbor2){
								notReady = true;
							}else{	
								if(connected(nbor1, nbor2)){
									notReady = true;
								}else{
									// Rewire!!!
									// nbor[nbor1][site1Index] = site1
									// nbor[nbor2][site2Index] = site2
									int site1Index = noNbors;
									int site2Index = noNbors;
									for(int i = 0; i < noNbors; i++){
										if(nbor[nbor1][i]==site1) site1Index = i;
										if(nbor[nbor2][i]==site2) site2Index = i;
									}
									nbor[site1][nbor1Index] = site2;
									nbor[nbor1][site1Index] = nbor2;
									nbor[site2][nbor2Index] = site1;
									nbor[nbor2][site2Index] = nbor1;
									System.out.println("Rewired sites " + site1 + " and " + site2);
								}			
							}
						}
					}
				}while(notReady);
			}
		}
	}
	

	
	boolean connected(int site1, int site2){
		boolean connected = false;
		for(int j = 0; j < noNbors; j ++){//If site i is already connected to new random site, choose another
			if (nbor[site1][j]==site2)
				connected = true;
		}
		return connected;
	}

	void randomizeOrder(){
		boolean available [] = new boolean [N];
		for(int i = 0; i < N; i ++){
			available[i] = true;			//set all positions as initially available
			randomOrder[i] = -1;			//no sites have yet been set
		}
		for (int i = 0; i < N; i++){
			while(randomOrder[i]==-1){
				int s = (int)(random.nextDouble()*N);
				if (available[s]) randomOrder[i] = s; 
			}
		}
	}

}
