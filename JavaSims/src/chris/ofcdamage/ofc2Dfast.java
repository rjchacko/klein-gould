package chris.ofcdamage;

import java.util.Random;
import chris.util.LatticeNeighbors;
import scikit.jobs.params.Parameters;

public class ofc2Dfast {

	@SuppressWarnings("unused")
	private double sr0, sf0, a0, dsr, dsf, da, sr[], sf[], stress[], sbar[], data[][];
	private int L, N, R, nbArray[], nbSeed, fs[], GR, qN;
	private boolean srn, sfn, an;
	@SuppressWarnings("unused")
	private String outdir, bcs;
	private Random rand;
	private LatticeNeighbors nbs;
	private static int dlength = 100000;
	private static int dcat = 4;
	
	
	public ofc2Dfast(Parameters params){
		
		constructor(params);
		return;
	}
	
	public void constructor(Parameters params){
		// constructor for ofc2Dfast object
		
		int seed;
		String shape;
		
		// read in parameters of model
		outdir = params.sget("Data Directory");
		seed   = params.iget("Random Seed");
		shape  = params.sget("Interaction Shape");
		R      = params.iget("Interaction Radius (R)");
		L      = params.iget("Lattice Size");
		bcs    = params.sget("Boundary Condtions");
		sf0     = params.fget("Failure Stress (\u03C3_f)");
		sr0    = params.fget("Residual Stress (\u03C3_r)");
		dsr    = params.fget("\u03C3_r width");
		a0     = params.fget("Dissipation (\u03B1)");
		da     = params.fget("\u03B1 width");
		
		// initialize variables
		rand = new Random(seed);
		
		srn = (dsr > 0);
		sfn = (dsr > 0);
		an  = (da > 0);
		
		N = L*L;
		
		sr      = new double[N];
		sf      = new double[N];
		stress  = new double[N];
		sbar    = new double[N];
		nbArray = new int[N];
		fs      = new int[3*N];
		data    = new double[dcat][dlength];
		
		for(int jj = 0 ; jj < N ; jj++){
			sr[jj]     = srn ? sr0+2*dsr*(rand.nextDouble()-0.5) : sr0;
			sf[jj]     = sfn ? sf0+2*dsf*(rand.nextDouble()-0.5) : sf0;
			stress[jj] = sr[jj] + (sf[jj]-sr[jj])*rand.nextDouble();
		}
				
		// set up the lattice neighbors array
		nbArray = setupNBS(shape);
		qN      = nbArray.length;
		
		return;
	}
	
	private int[] setupNBS(String shape){

		if(shape.equals("All Sites")){
			nbs = new LatticeNeighbors((int) L,(int) L,0,N/2,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.All);
			return nbs.get(nbSeed);	
		}

		if (bcs.equals("Bordered")){
			nbSeed = (int)((1 + L)*(L/2));
			if(shape.equals("Circle")){
				nbs = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Circle);
			}
			else if(shape.equals("Square")){
				nbs = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Square);
			}
			else if(shape.equals("Diamond")){
				nbs = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Diamond);
			}
		}
		else{
			nbSeed = 0;
			if(shape.equals("Circle")){
				nbs = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Circle);
			}
			else if(shape.equals("Square")){
				nbs = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Square);
			}
			else if(shape.equals("Diamond")){
				nbs = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Diamond);
			}
		}

		return nbs.get(nbSeed);	
	}
	
	public void evolve(int mct, boolean takedata){
		// evolve the ofc model starting with a stable configuration
		// until the next stable configuration is reached.
		//
		// mct is the monte carlo time step or the number of forced failures
		//
		// takedata specifies whether or not to record data
		
		double dsigma, release;
		int jjmax, index, newindex, a, b, tmpfail, tmpnb;
		index = 0;
		newindex = index;
		
		// force failure in the zero velocity limit
		jjmax = 0;
		for (int jj = 0 ; jj < N ; jj++){
			if( (sf[jj]-stress[jj]) < (sf[jjmax]-stress[jjmax])) jjmax = jj;
		}
		dsigma = sf[jjmax]-stress[jjmax];
		for (int jj = 0 ; jj < N ; jj++){
			stress[jj] += dsigma;
		}
		
		fs[newindex++] = jjmax;
		GR             = 1;
		
		while(newindex > index){
			a     = index;
			b     = newindex;
			index = newindex;
			for (int jj = a ; jj < b ; jj++){
				tmpfail = fs[jj];
				release = (1-nextAlpha())*(stress[tmpfail]-sr[tmpfail])/qN;
				for(int kk = 0 ; kk < qN ; kk++){
					tmpnb = nbs.getJ(fs[jj],nbSeed,nbArray,kk);
					stress[tmpnb] += release;
					if(stress[tmpnb] > sf[tmpnb]) fs[newindex++] = tmpnb;
				}
				GR += newindex - index;
				resetSite(tmpfail);
			}
		}
		return;
	}
	
	private double nextAlpha(){
		
		return an ? a0 + 2*da*(rand.nextDouble()-0.5) : a0;
	}
	
	private void resetSite(int site){
		
		sr[site]     = srn ? sr0+2*dsr*(rand.nextDouble()-0.5) : sr0;
		sf[site]     = sfn ? sf0+2*dsf*(rand.nextDouble()-0.5) : sf0;
		stress[site] = sr[site];
		return;
	}

}
