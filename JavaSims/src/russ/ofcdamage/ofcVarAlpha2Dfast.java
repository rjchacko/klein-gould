package russ.ofcdamage;

import scikit.jobs.params.Parameters;

public class ofcVarAlpha2Dfast extends ofc2Dfast {

	private double a0[], da[];
	private boolean an;
	
	public ofcVarAlpha2Dfast(Parameters params) {
		
		super(params);
		ofcVarAlpha2DfastConstructor(params);
		return;
	}

	public void ofcVarAlpha2DfastConstructor(Parameters params) {
		
		a0 = new double[N];
		da = new double[N];
		// fill a0[] using some criterion
		
		return;
	}

	// rewrite all methods which use a0 or da and replace these
	// variables with a0[] and da[]

	
	public void evolve(int mct, boolean takedata){
		// evolve the ofc model starting with a stable configuration
		// until the next stable configuration is reached.
		//
		// mct is the monte carlo time step or the number of forced failures
		//
		// takedata specifies whether or not to record data
		
		int a,b, tmpfail, tmpnb;
		double release;
		
		// force failure in the zero velocity limit
		forceZeroVel(mct, takedata, true);
		GR = 1; // the seed site
		
		// discharge site and repeat until lattice is stable
		while(newindex > index){
			a     = index;
			b     = newindex;
			index = newindex;
			for (int jj = a ; jj < b ; jj++){
				tmpfail = fs[jj];
				release = (1-nextAlpha(tmpfail))*(stress[tmpfail]-sr[tmpfail])/qN;
				for(int kk = 0 ; kk < qN ; kk++){
					tmpnb = getNbr(tmpfail,kk);
					if(tmpnb == -1 || failed[tmpnb]) continue; // -1 is returned if neighbor is self or is off lattice for open BC
					stress[tmpnb] += release;
					if(stress[tmpnb] > sf[tmpnb]){
						fs[newindex++] = tmpnb;	
						failSite(tmpnb);
					}
				}
				resetSite(tmpfail);
			}
		}
		return;
	}

	protected double nextAlpha(int id){
		
		return an ? a0[id] + 2*da[id]*(getRand().nextDouble()-0.5) : a0[id];
	}
	
	public double getAlpha0(int id){
		
		return a0[id];
	}

}
