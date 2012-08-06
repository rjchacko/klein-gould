package rachele.ofcNetworks;

import scikit.jobs.params.Parameters;

public class LongRangeLattice extends AbstractOFC_Lattice{
	
	int R; 		// Stress transfer range
		
	/**
	 * Constructor for Long Range Lattice Class
	 */
	public LongRangeLattice(Parameters params){
		setParams(params);
		initStress();
		setNbors();
	}
	
	void setParams(Parameters params){
		L = params.iget("L");
		N = L*L;
		stress = new double [N];
		alpha = params.fget("alpha");
		R = params.iget("R");
		noNbors = findNoCircleNbors();
		nbor = new int [N][noNbors];
		sitesToFail = new int [N+1];
		stressAdd = new double [N];
	}
	
	/**
	 * Boundary conditions for now are OPEN only!
	 */
	void setNbors(){
		int nborIndex = 0;
		for (int s = 0; s < N; s++){
			nborIndex = 0;
			int x = s%L;
			int y = s/L;
			for(int dy = -R; dy <= R; dy++){
				for(int dx = -R; dx <= R; dx++){
					double distance = Math.sqrt(dx*dx + dy*dy);
					if (distance <= R){
						int xx = (x+dx);						//define the location of the test site
						int yy = (y+dy);
						//System.out.println("site " + s);
						if(dx==0 & dy==0){						//exclude the site itself
							//System.out.println("no nbor set for site " + s);
						}else if(xx>=0 & xx<L & yy>=0 & yy<L){	//test site is on the lattice
							int nborSite = yy*L+xx;				//get index of test site
							nbor[s][nborIndex] = nborSite;		//set test site as next neighbor in nbor array
							//System.out.println("nbor index = " + s + " " + nborIndex  + " nbor = " + nborList[s][nborIndex]  + " nbor = " + dx + " " + dy);
							nborIndex += 1;
						}else{									//test site is not on lattice
							nbor[s][nborIndex] = -1;			//set next nbor as dead
							//System.out.println("nbor index = " + s + " "  + nborIndex  + " nbor = "   + nborList[s][nborIndex]  + " nbor = " + dx + " " + dy);
							nborIndex += 1;
						}

					}
				}
			}
		}
	}

	/**
	 * Returns the max number of circle neighbors.
	 */
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
		 count = count -1; // Subtract off the site itself.
		 return count;
	}
}
