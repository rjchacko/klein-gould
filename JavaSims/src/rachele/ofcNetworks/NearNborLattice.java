package rachele.ofcNetworks;

import scikit.jobs.params.Parameters;

public class NearNborLattice extends AbstractOFC_Lattice{

	/**
	 * Constructor for Near Nbor Lattice Class
	 */
	public NearNborLattice(Parameters params){
		
		setParams(params);
		initStress();
		setNearNbors();
	}
	
	void setParams(Parameters params){
		L = params.iget("L");
		N = L*L;
		stress = new double [N];
		alpha = params.fget("alpha");
		nbor = new int [N][4];			// 4 neighbors for nearest neighbor lattice
		noNbors = 4;
		sitesToFail = new int [N+1];
		stressAdd = new double [N];
	}

}
