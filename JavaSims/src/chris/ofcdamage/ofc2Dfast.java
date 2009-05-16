package chris.ofcdamage;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Random;

import scikit.jobs.params.Parameters;
import chris.util.LatticeNeighbors;
import chris.util.MathUtil;
import chris.util.PrintUtil;

public class ofc2Dfast {

	private double sr0, sf0, a0, dsr, dsf, da;
	protected double Omega, sr[], sf[], stress[],sbar[], data[][];
	private int L, R, nbArray[], nbSeed;
	protected int N, qN, fs[], GR, index, newindex, Ndead;
	private boolean srn, sfn, an;
	protected boolean failed[];
	private String outdir, bcs, bname;
	private Random rand;
	private LatticeNeighbors nbs;
	public static int dlength = 150000;
	public static int dcat = 4;
	private static DecimalFormat fmtI = new DecimalFormat("0000");
	
	public ofc2Dfast(Parameters params){
		
		constructor_ofc2Dfast(params);
		return;
	}
	
	public void constructor_ofc2Dfast(Parameters params){
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
		srn  = (dsr > 0);
		sfn  = (dsf > 0);
		an   = (da > 0);
		
		N     = L*L;
		Ndead = 0;
		
		sr      = new double[N];
		sf      = new double[N];
		stress  = new double[N];
		sbar    = new double[N];
		nbArray = new int[N];
		fs      = new int[2*N];
		data    = new double[dcat][dlength];
		failed  = new boolean[N];
		
		for(int jj = 0 ; jj < N ; jj++){
			sr[jj]     = srn ? sr0+2*dsr*(getRand().nextDouble()-0.5) : sr0;
			sf[jj]     = sfn ? sf0+2*dsf*(getRand().nextDouble()-0.5) : sf0;
			stress[jj] = sr[jj] + (sf[jj]-sr[jj])*getRand().nextDouble();
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
		
		int a,b, tmpfail, tmpnb;
		double release;
		GR = 0;
	
		// force failure in the zero velocity limit
		forceZeroVel(mct, takedata);
		
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
						failSite(tmpnb);
					}
				}
				resetSite(tmpfail);
			}
		}
		return;
	}
	
	protected void forceZeroVel(int mct, boolean takedata){
		// force failure in the zero velocity limit

		double dsigma, tmpbar;
		int jjmax;
		index = 0;
		newindex = index;
		
		jjmax = 0;
		if(takedata){
			tmpbar = 0;
			Omega  = 0;
			for (int jj = 0 ; jj < N ; jj++){ //use this loop to calculate the metric PART 1
				// find next site to fail
				if( ((sf[jj]-stress[jj]) < (sf[jjmax]-stress[jjmax])) && !failed[jj] ) jjmax = jj;
				
				// calculate metric (PART 1)
				sbar[jj] += MathUtil.bool2bin(!failed[jj])*stress[jj];
				tmpbar   += sbar[jj];
			}
			dsigma = sf[jjmax]-stress[jjmax];
			tmpbar = tmpbar / N;
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
			saveData(mct, true);
			fs[newindex++] = jjmax;
			failSite(jjmax,mct);
			
		}
		else{
			for (int jj = 0 ; jj < N ; jj++){
				// find next site to fail
				if( ((sf[jj]-stress[jj]) < (sf[jjmax]-stress[jjmax])) && !failed[jj] ) jjmax = jj;
			}
			// add stress to fail site
			dsigma = sf[jjmax]-stress[jjmax];
			for (int jj = 0 ; jj < N ; jj++){
				stress[jj] += MathUtil.bool2bin(!failed[jj])*dsigma;
			}
			fs[newindex++] = jjmax;
			failSite(jjmax);
		}

				
		return;
	}

	protected void failSite(int index, int mct){
		
		failSite(index);
		return;
	}
	
	protected void failSite(int index){
		
		GR++;
		failed[index]  = true;
		return;
	}
	
	protected double nextAlpha(){
		
		return an ? a0 + 2*da*(getRand().nextDouble()-0.5) : a0;
	}
	
	protected void resetSite(int site){

		stress[site] = sr[site];
		sr[site]     = nextSr(site);
		sf[site]     = nextSf(site);
		failed[site] = false;
		return;
	}
	
	protected double nextSr(int site){
	
		return srn ? sr0+2*dsr*(getRand().nextDouble()-0.5) : sr0;
	}
	
	protected double nextSf(int site){

		return sfn ? sf0+2*dsf*(getRand().nextDouble()-0.5) : sf0;
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
	
	protected void saveData(int mct, boolean EQ){
		
		data[0][mct%dlength] = 1/Omega;
		data[1][mct%dlength] = GR;
		data[2][mct%dlength] = MathUtil.bool2bin(EQ);
		data[3][mct%dlength] = 0; 
		return;
	}
	
	private void configOF(){
		
		try{
			File file = new File(outdir+File.separator+bname+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print("Time");
			pw.print("\t");
			pw.print("Inverse Metric");
			pw.print("\t");
			pw.print("Shower Size");
			pw.print("\t");
			pw.print("EQ Mode");
			pw.print("\t");
			pw.print("Phi");
			pw.println();
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		
		return;
	}
	
	public void PrintParams(String fout, Parameters prms){
		
		PrintUtil.printlnToFile(fout,prms.toString());
		PrintUtil.printlnToFile(fout,"Neighbors = " + fmtI.format(qN));
		return;
	}
	
	public String getOutdir(){
		
		return outdir;
	}
	
	public String getBname(){
		
		return bname;
	}
	
	public void setBname(String bn){
		
		bname = bn;
		return;
	}

	protected Random getRand() {
	
		return rand;
	}
	
	protected int getNbr(int ffs, int nbindex){
		
		return nbs.getJ(ffs,nbSeed,nbArray,nbindex);
	}
	
	public void clearData(){
		
		for(int jj = 0 ; jj < N ; jj++){
			for (int kk = 0 ; kk < dcat ; kk++){
				data[kk][jj] = 0;
			}
			sbar[jj] = 0;
		}
		return;
	}
}
