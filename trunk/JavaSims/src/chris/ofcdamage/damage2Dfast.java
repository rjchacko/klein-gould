package chris.ofcdamage;

import chris.util.MathUtil;
import chris.util.PrintUtil;
import scikit.jobs.params.Parameters;

public class damage2Dfast extends ofc2Dfast{
	
	private int liveNbs[], Nl0, Lives[];
	
	
	public damage2Dfast(Parameters params){
		
		super(params); //call to ofc2Dfast's constructor
		contructor_damage2Dfast(params);
	
		return;
	}
	
	public void contructor_damage2Dfast(Parameters params){
		
		int dN;
		
		Nl0 = params.iget("Number of Lives");
		dN  = params.iget("NL width");
		fs  = new int[3*N];
		
		Ndead   = 0;
		liveNbs = new int[N];
		Lives   = new int[N];
		
		if(Nl0 - dN < 0){	// cannot have negative lives
			dN = Nl0;
			params.set("Full width of NL", dN);
		}
		
		for (int jj = 0 ; jj < N ; jj++){
			liveNbs[jj] = qN;
			Lives[jj]   = (dN > 0) ? Nl0 + getRand().nextInt(2*dN+1)-dN : Nl0; // rand in [0 , 2*dN]
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
		
		double release;
		int a, b, tmpfail, tmpnb;
		boolean lastlife = false;
		
		// force failure
		forceFailure(mct, takedata);
		GR = 1; // the seed site
		
		// discharge forced site(s) and repeat until lattice is stable
		while(newindex > index){
			a     = index;
			b     = newindex;
			index = newindex;
			for (int jj = a ; jj < b ; jj++){
				tmpfail = fs[jj];
				if(liveNbs[tmpfail] == 0){
					if(Lives[tmpfail] == 1){
						killSite(tmpfail);
					}
					else{
						resetSiteD(tmpfail);
					}
					continue;
				}
				release = (1-nextAlpha())*(stress[tmpfail]-sr[tmpfail])/liveNbs[tmpfail];
				if(Lives[tmpfail] == 1) lastlife = true;
				for(int kk = 0 ; kk < qN ; kk++){
					tmpnb = getNbr(fs[jj],kk);
					if(tmpnb == -1 || failed[tmpnb]) continue; // -1 is returned if neighbor is self or is off lattice for open BC
					if(lastlife) liveNbs[tmpnb]--;
					stress[tmpnb] += release;
					if((stress[tmpnb] > sf[tmpnb])){
						fs[newindex++] = tmpnb;	
						failSite(tmpnb);
					}
				}
				if(lastlife){
					lastlife = false;
					killSite(tmpfail);
				}
				else{
					resetSiteD(tmpfail);
				}
			}
		}
		return Ndead;
	}
	
	protected void forceFailure(int mct, boolean takedata){
		// force failure in the zero velocity limit

		double dsigma, tmpbar, dsmin;
		int jjmax;
		index = 0;
		newindex = index;
		dsmin    = 1e5;
		
		// force failure in the zero velocity limit
		jjmax = 0;
		if(takedata){
			tmpbar = 0;
			Omega  = 0;
			for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 1
				// find next site to fail (NB, must be alive)
				if( (sf[jj]-stress[jj]) < dsmin && !failed[jj]){
					jjmax = jj;
					dsmin = sf[jj]-stress[jj];
				}
				// calculate metric (PART 1)
				sbar[jj] += MathUtil.bool2bin(!failed[jj])*stress[jj];
				tmpbar   += MathUtil.bool2bin(!failed[jj])*sbar[jj];
			}
			dsigma = sf[jjmax]-stress[jjmax];
			tmpbar = tmpbar / (N-Ndead);
			Omega  = 0;
			for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 2
				// add stress to fail site
				stress[jj] += MathUtil.bool2bin(!failed[jj])*dsigma;
				//calculate metric (PART 2)
				Omega += MathUtil.bool2bin(!failed[jj])*(sbar[jj] - tmpbar)*(sbar[jj] - tmpbar);
			}
			//calculate metric (PART 3)
			Omega = Omega/((double)(mct)*(double)(mct)*(double)(N-Ndead));

			// save and/or write data
			if(mct%dlength == 0 && mct > 0){
				writeData(mct);
			}
			saveData(mct, false);
		}
		else{
			for (int jj = 0 ; jj < N ; jj++){
				// find next site to fail
				if( (sf[jj]-stress[jj]) < dsmin && !failed[jj]){
					jjmax = jj;
					dsmin = sf[jj]-stress[jj];
				}
			}
			// add stress to fail site
			dsigma = sf[jjmax]-stress[jjmax];
			for (int jj = 0 ; jj < N ; jj++){
				stress[jj] += MathUtil.bool2bin(!failed[jj])*dsigma;
			}
		}
				
		fs[newindex++] = jjmax;
		failSite(jjmax);
		
		return;
	}
	
	protected void failSite(int index){
		
		GR++;
		failed[index]  = true;
		return;
	}
	
	
	private void killSite(int site){
		
		Ndead++;
		sr[site]     = 0;
		sf[site]     = 1;
		stress[site] = 0;
		Lives[site]  = 0;
		failed[site] = true;
		return;
	}
	
	protected void resetSiteD(int site){
		
		super.resetSite(site);
		Lives[site]--;
		return;
	}
	
	
	public int getN(){
		
		return N;
	}
	
	public double[] getStress(){
		
		return stress;
	}
	
	public double getStress(int site){
		
		return (site < N) ? stress[site] : -1;
	}
	
	public int[] getDorA(){
		
		int[] ret = new int[N];
		
		for(int jj = 0 ; jj < N ; jj++){
			ret[jj] = (Lives[jj] > 0) ? 1 : 0;
		}
		return ret;
	}
	
	protected void saveData(int mct, boolean EQ){
		
		data[0][mct%dlength] = 1./Omega;
		data[1][mct%dlength] = GR;
		data[2][mct%dlength] = MathUtil.bool2bin(EQ);
		data[3][mct%dlength] = 1. - (double)(Ndead)/(double)(N); 
		return;
	}
	
	public double getData(int mct, int dindex){
		
		if(dindex >= ofc2Dfast.dcat) return -77;
		
		return data[dindex][mct%dlength];
	}
	
	public void printFitParams(String fout, double slope, double offset, double width){
		
		PrintUtil.printlnToFile(fout,"slope = ", slope);
		PrintUtil.printlnToFile(fout,"intercept = ", offset);
		PrintUtil.printlnToFile(fout,"width = ", width);
		return;
	}
	
	public int getLN(int site){
		
		return (site < N) ? liveNbs[site] : -1;
	}
	
	public boolean isAlive(int site){
		
		return (!failed[site]);
	}
	
}
