package chris.ofcdamage;

import scikit.jobs.params.Parameters;

public class damage2Dfast extends ofc2Dfast{
	
	private int Ndead, liveNbs[], Nl0, Lives[];
	
	
	public damage2Dfast(Parameters params){
		
		super(params); //call to ofc2Dfast's constructor
		contructor_damage2Dfast(params);
	
		return;
	}
	
	public void contructor_damage2Dfast(Parameters params){
		
		int dN;
		
		Nl0 = params.iget("Number of Lives");
		dN  = params.iget("Full width of NL");
		
		Ndead   = 0;
		liveNbs = new int[N];
		Lives   = new int[N];
		
		if(Nl0 - dN > 0){	// cannot have negative lives
			dN = Nl0;
			params.set("Full width of NL", dN);
		}
		
		for (int jj = 0 ; jj < N ; jj++){
			liveNbs[jj] = qN;
			Lives[jj]   = (dN > 0) ? Nl0 + getRand().nextInt(2*dN)-dN : Nl0;
		}
		
		if(dN == Nl0){
			for (int jj = 0 ; jj < N ; jj++){
				if(Lives[jj] == 0){
					Ndead++;
					failed[jj] = true;
					for (int kk = 0 ; kk < qN ; kk++){
						liveNbs[getNbr(jj,kk)]--;
					}
				}
			}
		}
		
		return;
	}
	
	
	public int evolveD(int mct, boolean takedata){
		// evolve the ofc model *IN DAMAGE MODE* starting with a stable 
		//configurationu ntil the next stable configuration is reached.
		//
		// mct is the monte carlo time step or the number of forced failures
		//
		// takedata specifies whether or not to record data
		
		
		double dsigma, release, tmpbar;
		int jjmax, index, newindex, a, b, tmpfail, tmpnb;
		boolean lastlife = false;
		index = 0;
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
			}
			dsigma = sf[jjmax]-stress[jjmax];
			tmpbar = tmpbar / N;
			Omega  = 0;
			for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 2
				// add stress to fail site
				stress[jj] += dsigma;
				
				//calculate metric (PART 2)
				Omega += (sbar[jj] - tmpbar)*(sbar[jj] - tmpbar);
			}
			//calculate metric (PART 3)
			Omega = Omega/((double)(mct)*(double)(mct)*(double)(N));

			// save and/or write data
			if(mct%dlength == 0 && mct > 0){
				writeData(mct);
			}
			saveData(mct, false);
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
				
		// discharge site and repeat until lattice is stable
		while(newindex > index){
			a     = index;
			b     = newindex;
			index = newindex;
			for (int jj = a ; jj < b ; jj++){
				tmpfail = fs[jj];
				if(liveNbs[tmpfail] == 0){
					// NOT QUITE -- FIX ME!
					if(Lives[tmpfail] == 1){
						killSite(tmpfail);
					}
					else{
						resetSite(tmpfail);
					}
					continue;
				}
				release = (1-nextAlpha())*(stress[tmpfail]-sr[tmpfail])/liveNbs[tmpfail];
				if(Lives[tmpfail] == 1) lastlife = true;
				for(int kk = 0 ; kk < qN ; kk++){
					tmpnb = getNbr(fs[jj],kk);
					if(tmpnb == -1 || failed[tmpnb]) continue; // -1 is returned if neighbor is self or is off lattice for open BC
					stress[tmpnb] += release;
					if(stress[tmpnb] > sf[tmpnb]){
						fs[newindex++] = tmpnb;	
						failed[tmpnb] = true;
					}
					if(lastlife) liveNbs[tmpnb]--;
				}
				GR += newindex - index;
				if(lastlife){
					lastlife = false;
					killSite(tmpfail);
				}
				else{
					resetSite(tmpfail);
				}
			}
		}
		return Ndead;
	}
	
	private void killSite(int site){
		
		Ndead++;
		sr[site]     = 0;
		sf[site]     = 0;
		stress[site] = 0;
		failed[site] = true;
		return;
	}
	
	public int getN(){
		
		return N;
	}
	
	public double[] getStress(){
		
		return stress;
	}
	
	public int[] getDorA(){
		
		int[] ret = new int[N];
		
		for(int jj = 0 ; jj < N ; jj++){
			ret[jj] = (Lives[jj] > 0) ? 1 : 0;
		}
		return ret;
	}
	
}
