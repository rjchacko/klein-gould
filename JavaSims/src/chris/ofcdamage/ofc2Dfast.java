package chris.ofcdamage;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;

import scikit.jobs.params.Parameters;
import chris.util.LatticeNeighbors;

public class ofc2Dfast {

	private double sr0, sf0, a0, dsr, dsf, da, sr[], sf[], stress[], sbar[], data[][], Omega;
	private int L, N, R, nbArray[], nbSeed, fs[], GR, qN;
	private boolean srn, sfn, an, failed[];
	private String outdir, bcs, bname;
	private Random rand;
	private LatticeNeighbors nbs;
	public static int dlength = 150000;
	private static int dcat = 2;
	
	
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
		bname  = params.sget("Data File");
		seed   = params.iget("Random Seed");
		shape  = params.sget("Interaction Shape");
		R      = params.iget("Interaction Radius (R)");
		L      = params.iget("Lattice Size");
		bcs    = params.sget("Boundary Condtions");
		sf0    = params.fget("Failure Stress (\u03C3_f)");
		dsf    = params.fget("\u03C3_f width");
		sr0    = params.fget("Residual Stress (\u03C3_r)");
		dsr    = params.fget("\u03C3_r width");
		a0     = params.fget("Dissipation (\u03B1)");
		da     = params.fget("\u03B1 width");
		
		// initialize variables
		rand = new Random(seed);
		
		srn = (dsr > 0);
		sfn = (dsf > 0);
		an  = (da > 0);
		
		N = L*L;
		
		sr      = new double[N];
		sf      = new double[N];
		stress  = new double[N];
		sbar    = new double[N];
		nbArray = new int[N];
		fs      = new int[3*N];
		data    = new double[dcat][dlength];
		failed  = new boolean[N];
		
		for(int jj = 0 ; jj < N ; jj++){
			sr[jj]     = srn ? sr0+2*dsr*(rand.nextDouble()-0.5) : sr0;
			sf[jj]     = sfn ? sf0+2*dsf*(rand.nextDouble()-0.5) : sf0;
			stress[jj] = sr[jj] + (sf[jj]-sr[jj])*rand.nextDouble();
			failed[jj] = false;
		}

		// set up the lattice neighbors array
		if(shape.equals("Open")) shape = "Bordered";
		nbArray = setupNBS(shape);
		qN      = nbArray.length;
		
		// set up output file
		configOF();
		
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
		
		double dsigma, release, tmpbar;
		int jjmax, index, newindex, a, b, tmpfail, tmpnb;
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
			saveData(mct);
			
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
				release = (1-nextAlpha())*(stress[tmpfail]-sr[tmpfail])/qN;
				for(int kk = 0 ; kk < qN ; kk++){
					tmpnb = nbs.getJ(fs[jj],nbSeed,nbArray,kk);
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
	
	private double nextAlpha(){
		
		return an ? a0 + 2*da*(rand.nextDouble()-0.5) : a0;
	}
	
	private void resetSite(int site){
		
		sr[site]     = srn ? sr0+2*dsr*(rand.nextDouble()-0.5) : sr0;
		sf[site]     = sfn ? sf0+2*dsf*(rand.nextDouble()-0.5) : sf0;
		stress[site] = sr[site];
		failed[site] = false;
		return;
	}

	public void writeData(int mct){
		
		int ub;
		int offset = (int)((mct-1)/dlength);
		offset = offset*dlength;

		ub = mct%dlength;
		if(ub==0) ub = dlength;
		
		try{
			File file = new File(outdir+File.separator+bname+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int jj = 0 ; jj < ub ; jj++){
				pw.print(jj+offset);
				pw.print("\t");
				for (int kk = 0 ; kk < dcat ; kk++){
					pw.print(data[kk][jj]);
					pw.print("\t");
				}
				pw.println();
			}			
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		return;
	}
	
	private void saveData(int mct){
		
		data[0][mct%dlength] = Omega;
		data[1][mct%dlength] = GR;
		return;
	}
	
	private void configOF(){
		
		try{
			File file = new File(outdir+File.separator+bname+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print("Time");
			pw.print("\t");
			pw.print("Metric");
			pw.println();
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		
		return;
	}
	
	
}
