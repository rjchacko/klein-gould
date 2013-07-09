package rachele.ofcNetworks;

import scikit.jobs.params.Parameters;

public class SmallWorldLatticeFixedNoRewire extends AbstractOFC_Lattice{

	public double percentRewired;
	
	double prob;
	int randomOrder [];

	
	/**
	 * Constructor for Small World Lattice Class
	 */
	public SmallWorldLatticeFixedNoRewire(Parameters params){
		
		setParams(params);
		initStress();
		setNearNbors();
		rewireSetNo();
	}
	
	void setParams(Parameters params){
		random.setSeed(params.iget("Random Seed"));
		L = params.iget("L");
		N = L*L;
		prob = params.fget("Rewire Probability");
		System.out.println("Rewire Probability = " + prob);
		stress = new double [N];
		alpha = params.fget("alpha");
		nbor = new int [N][4];			// 4 neighbors for nearest neighbor lattice
		noNbors = 4;
		sitesToFail = new int [N+1];
		stressAdd = new double [N];
		randomOrder = new int [N];
	}
	void rewireSetNo(){
		randomizeOrder();
		int noToRewire = (int)(prob*N);
		System.out.println("No of nodes to rewire = " + noToRewire);
		for(int k = 0; k < noToRewire; k++){
			int site1 = randomOrder[k];


			int site2 = -1;	
			int nbor1Index = -1;
			int nbor2Index = -1;
			
			//Select a random link from ith node to rewire.
			int nbor1 = -1;
			while(nbor1==-1){
				nbor1Index = (int)(random.nextDouble()*noNbors); 
				nbor1 = nbor[site1][nbor1Index];					//Makes sure you didn't pick a dead site?
			}
			
			//choose a random second site
			site2 = (int)(random.nextDouble()*N);  
			while(site1 == site2){						
				site2 = (int)(random.nextDouble()*N); 	//Make sure you didn't choose the original site itself.
			}
			
			//Choose a random neighbor of that site
			int nbor2 = -1;
			while(nbor2 == -1){
				nbor2Index = (int)(random.nextDouble()*noNbors);			
				nbor2 = nbor[site2][nbor2Index];
			}
			// Rewire!!!
			// nbor[nbor1][site1Index] = site1
			// nbor[nbor2][site2Index] = site2
			int site1Index = noNbors;
			int site2Index = noNbors;
			for(int i = 0; i < noNbors; i++){
				//Possible problem if these nbors are boundary sites
				if(nbor[nbor1][i]==site1) site1Index = i;
				if(nbor[nbor2][i]==site2) site2Index = i;
			}
			nbor[site1][nbor1Index] = site2;
			nbor[site2][nbor2Index] = site1;
			nbor[nbor2][site2Index] = nbor1;
			nbor[nbor1][site1Index] = nbor2;
			//System.out.println("Rewired sites " + site1 + " and " + site2);
			//System.out.println("Not Ready = " + notReady);
		}
	}
	
//	void rewire(){
//		randomizeOrder();
//		int noRewired = 0;
//		for (int k = 0; k < N; k++){
//			int site1 = randomOrder[k];
//			double r = random.nextDouble();
//			if (r < prob){												//test for rewire
//				int site2 = -1;	
//				int nbor1Index = -1;
//				int nbor2Index = -1;
//				boolean notReady = false;
//				do{
//					int nbor1 = -1;
//					while(nbor1==-1){
//						nbor1Index = (int)(random.nextDouble()*noNbors); 
//						nbor1 = nbor[site1][nbor1Index];
//					}
//					site2 = (int)(random.nextDouble()*N);  //choose a random site
//					if (site1 == site2){						//Make sure you haven't picked the original site
//						notReady = true;
//					} else{
//						if(connected(site1, site2)){
//							notReady = true;
//						}else{
//							int nbor2 = -1;
//							while(nbor2 == -1){
//								nbor2Index = (int)(random.nextDouble()*noNbors);			
//								nbor2 = nbor[site2][nbor2Index];
//							}
//							if(nbor1 == nbor2){
//								notReady = true;
//							}else{	
//								if(connected(nbor1, nbor2)){
//									notReady = true;
//								}else{
//									// Rewire!!!
//									// nbor[nbor1][site1Index] = site1
//									// nbor[nbor2][site2Index] = site2
//									int site1Index = noNbors;
//									int site2Index = noNbors;
//									for(int i = 0; i < noNbors; i++){
//										//Possible problem if these nbors are boundary sites
//										if(nbor[nbor1][i]==site1) site1Index = i;
//										if(nbor[nbor2][i]==site2) site2Index = i;
//									}
//									nbor[site1][nbor1Index] = site2;
//									nbor[site2][nbor2Index] = site1;
//									nbor[nbor2][site2Index] = nbor1;
//									nbor[nbor1][site1Index] = nbor2;
//									//System.out.println("Rewired sites " + site1 + " and " + site2);
//									//System.out.println("Not Ready = " + notReady);
//									noRewired ++;
//									notReady = false;
//								}			
//							}
//						}
//					}
//				}while(notReady);
//			}
//		}
//		percentRewired = (double)noRewired/(double)N;
//		System.out.println("Percent Rewired = " + percentRewired);
//	}
//	

	
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
