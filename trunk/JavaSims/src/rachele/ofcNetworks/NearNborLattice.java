package rachele.ofcNetworks;

import scikit.jobs.params.Parameters;

public class NearNborLattice extends AbstractOFC_Lattice{

	/**
	 * Constructor for Near Nbor Lattice Class
	 */
	public NearNborLattice(Parameters params){
		
		setParams(params);
		initStress();
		setNbors();
	}
	
	void setParams(Parameters params){
		L = params.iget("L");
		N = L*L;
		stress = new double [N];
		alpha = params.fget("alpha");
		nbor = new int [N][4];			// 4 neighbors for nearest neighbor lattice
		noNbors = 4;
		sitesToFail = new int [N];
	}
	
	/**
	 * Boundary conditions for now are OPEN only!
	 */
	void setNbors(){
		for(int i = 0; i < N; i++){
			int x = i % L;
			int y = i / L;
			
			// left neighbor
			if (x == 0){
				nbor[i][0] = -1;
			}else{
				nbor[i][0] = i - 1;
			}
			// bottom neighbor
			if (y == 0){
				nbor[i][1] = -1;
			}else{
				nbor[i][1] = i - L;
			}
			// right neighbor
			if (x == L-1){
				nbor[i][2] = -1;
			}else{
				nbor[i][2] = i + 1;
			}
			// top neighbor
			if (y == L-1){
				nbor[i][3] = -1;
			}else{
				nbor[i][3] = i + L;
			}
		}
	}


	

	

	

}
