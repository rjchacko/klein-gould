package chris.ofcdamage;

import scikit.jobs.Job;
import scikit.jobs.params.Parameters;
import chris.util.MathUtil;

public class heal2Dfast extends ofc2Dfast{
	
	private int ht, dht, heal[][], htstop[], tmax, Ndead;
	private double phi;
	private boolean htN;

	public heal2Dfast(Parameters params) {
	
		super(params);
		constructor_heal2Dfast(params);
		return;
	}

	public void constructor_heal2Dfast(Parameters params){
		
		ht     = params.iget("Heal Time");
		dht    = params.iget("HT Width");
		tmax   = ht + dht;
		heal   = new int[tmax][N]; 
		htstop = new int[tmax];
		phi    = 1.;
		htN    = (dht > 0);
		Ndead  = 0;
		
		for (int jj = 0 ; jj < tmax ; jj++){
			htstop[jj] = 0;
		}
		
		return;
	}
	
	public double evolveH(int mct, boolean takedata){
		// evolve the ofc model *IN DAMAGE MODE* starting with a stable 
		//configuration until the next stable configuration is reached.
		//
		// mct is the monte carlo time step or the number of forced failures
		//
		// takedata specifies whether or not to record data
		double release;
		int a, b, tmpfail, LN, tmpnb;
		
		//heal appropriate sites
		healSites(mct);

		
		//force failure
		forceZeroVel(mct, takedata, false);
		
		// discharge forced site(s) and repeat until lattice is stable
		while(newindex > index){
			a = index;
			b = newindex;
			index = newindex;
			for (int jj = a ; jj < b ; jj++){
				tmpfail = fs[jj];
				LN = 0;
				for (int kk = 0 ; kk < qN ; kk++){
					tmpnb = getNbr(tmpfail,kk);
					if(failed[tmpnb] || tmpnb == -1) continue;
					LN++;
				}
				if(LN == 0) continue;
				release = (1-nextAlpha())*(stress[tmpfail]-sr[tmpfail])/LN;
				for (int kk = 0 ; kk < qN ; kk++){
					tmpnb = getNbr(tmpfail,kk);
					if(failed[tmpnb] || tmpnb == -1) continue;
					stress[tmpnb] += release;
					if((stress[tmpnb] > sf[tmpnb])){
						fs[newindex++] = tmpnb;	
						failSite(tmpnb,mct);
					}
				}
//				resetSite(tmpfail);
//				// this is where you would, in principle, reset the site
			}
			Job.animate();  // TAKE PICTURE
		}
		
		phi = 1.-(double)(Ndead)/(double)(N);
		return phi;
	}
	
	
	protected void failSite(int index, int mct){
		int tmp;
		
		failSite(index);
		Ndead++;
		tmp = (mct + nextHT())%tmax;
		heal[tmp][htstop[tmp]++] = index;
		return;
	}
	
	private int nextHT(){
		
		return htN ? ht + getRand().nextInt(2*dht+1)-dht : ht; // rand in [0 , 2*dht]
	}
	
	private void healSites(int mct){
		
		for( int jj = 0 ; jj < htstop[mct%tmax] ; jj++){
			resetSite(heal[mct%tmax][jj]);
			Ndead--;
		}
		htstop[mct%tmax] = 0;
		return;
	}
	
	protected void saveData(int mct, boolean EQ){
		
		data[0][mct%dlength] = 1/Omega;
		data[1][mct%dlength] = GR;
		data[2][mct%dlength] = MathUtil.bool2bin(EQ);
		data[3][mct%dlength] = phi; 
		return;
	}
	
	public double[] getStress(){
		
		return stress;
	}
	
	public int[] getDorA(){
		
		int[] ret = new int[N];
		
		for(int jj = 0 ; jj < N ; jj++){
			ret[jj] = 1 - MathUtil.bool2bin(failed[jj]);
		}
		return ret;
	}
	
	public int getGR(){
		
		return GR;
	}
	
}
