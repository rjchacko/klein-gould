package chris.ofcdamage;

import scikit.jobs.params.Parameters;
import chris.util.PrintUtil;

public class damage2DfastCG extends ofc2Dfast{
	
	private int liveNbs[], Nl0, Lives[];
	
	
	public damage2DfastCG(Parameters params){
		
		super(params); //call to ofc2Dfast's constructor
		contructor_damage2DfastCG(params);
	
		return;
	}
	
	public void contructor_damage2DfastCG(Parameters params){
		
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
		forceZeroVel(mct, takedata, false);
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
	
	public int[] getDorA(){
		
		int[] ret = new int[N];
		
		for(int jj = 0 ; jj < N ; jj++){
			ret[jj] = (Lives[jj] > 0) ? 1 : 0;
		}
		return ret;
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
	
}
