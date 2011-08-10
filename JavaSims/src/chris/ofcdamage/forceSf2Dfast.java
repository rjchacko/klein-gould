package chris.ofcdamage;

import chris.util.MathUtil;
import scikit.jobs.params.Parameters;

public class forceSf2Dfast extends damage2Dfast{

	public forceSf2Dfast(Parameters params) {
		
		super(params);
		contructor_forceSf2Dfast(params);
		return;
	}

	public void contructor_forceSf2Dfast(Parameters params){
		
		return;
	}
	
	public void evolveCO(int mct, boolean takedata){
		// evolve the ofc model starting with a stable configuration
		// until the next stable configuration is reached.
		//
		// mct is the monte carlo time step or the number of forced failures
		//
		// takedata specifies whether or not to record data
		
		double release;
		int a, b, tmpfail, tmpnb;
		index = 0;
		newindex = index;
		
		// force failure
		forceFailure(mct, takedata);
		
		// discharge site and repeat until lattice is stable
		while(newindex > index){
			a     = index;
			b     = newindex;
			index = newindex;
			for (int jj = a ; jj < b ; jj++){
				tmpfail = fs[jj];
				release = (1-nextAlpha())*(stress[tmpfail]-sr[tmpfail])/qN;
				for(int kk = 0 ; kk < qN ; kk++){
					tmpnb = getNbr(fs[jj],kk);
					if(tmpnb == -1 || failed[tmpnb]) continue; // -1 is returned if neighbor is self or is off lattice for open BC
					stress[tmpnb] += release;
					if(stress[tmpnb] > sf[tmpnb]){
						fs[newindex++] = tmpnb;	
						failed[tmpnb] = true;
					}
				}
				GR += newindex - index;
				resetSite(tmpfail);
			}
		}
		return;
	}
	

	protected void forceFailure(int mct, boolean takedata){
		// force failure in the zero velocity limit

		double dsigma, tmpbar;
		int jjmax, ndt;
		index = 0;
		ndt   = 0;
		newindex = index;
		
		// force failure in the zero velocity limit
		jjmax = 0;
		if(takedata){
			tmpbar = 0;
			Omega  = 0;
			for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 1
				// find next site to fail
				if( (sf[jj]-stress[jj]) < (sf[jjmax]-stress[jjmax])) jjmax = jj;
				
				// calculate metric (PART 1)
				sbar[jj] += stress[jj];
				tmpbar   += sbar[jj];
				ndt      += MathUtil.bool2bin(ftt[jj]);
			}
			dsigma = sf[jjmax]-stress[jjmax];
			tmpbar = tmpbar / N;
			Omega  = 0;
			for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 2
				// reduce failure threshold to fail site
				sf[jj] -= dsigma;
				
				//calculate metric (PART 2)
				Omega += (sbar[jj] - tmpbar)*(sbar[jj] - tmpbar);
				ftt[jj] = false;
			}
			//calculate metric (PART 3)
			Omega = Omega/((double)(mct)*(double)(mct)*(double)(N));

			// save and/or write data
			if(mct%dlength == 0 && mct > 0){
				writeData(mct);
			}
			saveData(mct, false, dsigma, ndt);
		}
		else{
			for (int jj = 0 ; jj < N ; jj++){
				// find next site to fail
				if( (sf[jj]-stress[jj]) < (sf[jjmax]-stress[jjmax])) jjmax = jj;
			}
			// add stress to fail site
			dsigma = sf[jjmax]-stress[jjmax];
			for (int jj = 0 ; jj < N ; jj++){
				stress[jj] += dsigma;
			}
		}
				
		fs[newindex++] = jjmax;
		failed[jjmax]  = true;
		GR             = 1;
		
		return;
	}
	
	
	protected void forceZeroVel(int mct, boolean takedata, boolean eqmode){
		// force failure in the zero velocity limit

		double dsigma, tmpbar;
		int jjmax, ndt;
		index = 0;
		ndt   = 0;
		newindex = index;
	
		jjmax = 0;
		if(takedata){
			tmpbar = 0;
			Omega  = 0;
			for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 1
				hstrs.accum(stress[jj]);
				// find next site to fail
				if( ((sf[jj]-stress[jj]) < (sf[jjmax]-stress[jjmax])) && !failed[jj] ) jjmax = jj;
				// calculate metric (PART 1)
				sbar[jj] += MathUtil.bool2bin(!failed[jj])*stress[jj];
				tmpbar   += MathUtil.bool2bin(!failed[jj])*sbar[jj];
				// sum the sites from last time
				ndt += MathUtil.bool2bin(ftt[jj]);
			}
			dsigma = sf[jjmax]-stress[jjmax];
			tmpbar = tmpbar / (N-Ndead);
			//tmpbar = tmpbar / N;
			Omega  = 0;
			for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 2
				// add stress to fail site
				sf[jj] -= MathUtil.bool2bin(!failed[jj])*dsigma;
				//calculate metric (PART 2)
				Omega += MathUtil.bool2bin(!failed[jj])*(sbar[jj] - tmpbar)*(sbar[jj] - tmpbar);
				// reset ftt
				ftt[jj] = false;
			}
			//calculate metric (PART 3)
			//Omega = Omega/((double)(mct)*(double)(mct)*(double)(N-Ndead));
			Omega = Omega/((double)(mct)*(double)(mct)*(double)(N));

			
			// save and/or write data
			if(mct%dlength == 0 && mct > 0){
				writeData(mct);
			}
			saveData(mct, eqmode, dsigma, -77); //-77 to test if this is  called
		}
		else{
			for (int jj = 0 ; jj < N ; jj++){
				// find next site to fail
				if( ((sf[jj]-stress[jj]) < (sf[jjmax]-stress[jjmax])) && !failed[jj] ) jjmax = jj;
			}
			// add stress to fail site
			dsigma = sf[jjmax]-stress[jjmax];
			for (int jj = 0 ; jj < N ; jj++){
				sf[jj] -= MathUtil.bool2bin(!failed[jj])*dsigma;
			}
		}
		
		fs[newindex++] = jjmax;
		failSite(jjmax,mct);	
		return;
	}
	
	protected double nextSf(int site){

		return sf[site];
	}
	
	public double sfNow(){
		
		return sf[0];
	}
	
	
}
