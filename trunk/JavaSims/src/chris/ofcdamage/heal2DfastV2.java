package chris.ofcdamage;

import scikit.jobs.params.Parameters;

public class heal2DfastV2 extends damage2Dfast {

	protected int ht, dht, htm, eoa[], heal[][]; 	//eoa marks the end of the array 				
	protected boolean htN;							//(jj < eoa[time] in heal[time][jj])
	
	public heal2DfastV2(Parameters params) {

		super(params);
		head2DfastV2_constructor(params);
	}
	
	public void head2DfastV2_constructor(Parameters params){
		
		ht   = params.iget("Heal Time");
		dht  = params.iget("HT Width");
		htN  = (dht > 0);
		htm  = ht + dht;
		heal = new int[htm][N];
		eoa  = new int[htm];
		
		return;
	}

	public int evolveH(int mct, boolean takedata){
		// evolve the ofc model *IN DAMAGE MODE* starting with a stable 
		//configuration until the next stable configuration is reached.
		//
		// mct is the monte carlo time step or the number of forced failures
		//
		// takedata specifies whether or not to record data
	
		healSites(mct);
		return evolveD(mct, takedata);
	}
	
	protected void killSite(int index, int mct){

		killSite(index);
		int tmp = (mct + nextHT())%htm;
		heal[tmp][eoa[tmp]++] = index;
		return;
	}
	
	private int nextHT(){
		
		return htN ? ht + getRand().nextInt(2*dht+1)-dht : ht; // rand in [0 , 2*dht]
	}
	
	private void healSites(int mct){
		int tmp;
		for( int jj = 0 ; jj < eoa[mct%htm] ; jj++){
			int st = heal[mct%htm][jj];
			resetSite(st);
			Ndead--;
			Lives[st] = nextNumLives();
			for (int kk = 0 ; kk < qN ; kk++){
				tmp = getNbr(st,kk);
				if(tmp == -1) continue;
				liveNbs[tmp]++;
			}
		}
		eoa[mct%htm] = 0;
		return;
	}
	
	
	
}

