package chris.ofcdamage;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import scikit.jobs.params.Parameters;
import chris.util.LatticeNeighbors;
import chris.util.MathUtil;
import chris.util.PrintUtil;
import chris.util.Random;
import chris.util.SortUtil;

public class ofc2Dfast{

	private double sr0, sf0, a0, dsr, dsf, da;
	protected double Omega, sr[], sf[], stress[], sbar[], data[][];
	private int L, R, nbArray[], nbSeed;
	protected int N, qN, fs[], GR, index, newindex, Ndead;
	private boolean srn, sfn, an;
	protected boolean failed[], ftt[];
	private String outdir, bcs, bname;
	private Random rand;
	private LatticeNeighbors nbs;
	public static int dlength = 150000;
	public static int dcat = 6;
	private static DecimalFormat fmtI = new DecimalFormat("0000");
	
	public ofc2Dfast(Parameters params){
		
		constructor_ofc2Dfast(params);
		return;
	}

	public ofc2Dfast(Parameters params, double[] srIN, double[] sfIN, double sIN[], double sbIN[], double dIN[][], int grIN, boolean lsIN[]){
		
		constructor_ofc2Dfast(params);
		copy_sys(srIN, sfIN, sIN, sbIN, dIN, grIN, lsIN);
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
		fs      = new int[2*N];
		data    = new double[dcat][dlength];
		failed  = new boolean[N];
		ftt     = new boolean[N];

		
		for(int jj = 0 ; jj < N ; jj++){
			sr[jj]     = srn ? sr0+2*dsr*(getRand().nextDouble()-0.5) : sr0;
			sf[jj]     = sfn ? sf0+2*dsf*(getRand().nextDouble()-0.5) : sf0;
			stress[jj] = sr[jj] + (sf[jj]-sr[jj])*getRand().nextDouble();
			failed[jj] = false;
			ftt[jj]    = false;
		}

		// set up the lattice neighbors array
		if(shape.equals("Open")) shape = "Bordered";
		nbArray = setupNBS(shape);
		if(nbArray == null){
			qN = N;
		}
		else{
			qN = nbArray.length;
		}
		
		// set up output file
		configOF();
		
		return;
	}
	
	private void copy_sys(double[] srIN, double[] sfIN, double sIN[], double sbIN[], double dIN[][], int grIN, boolean lsIN[]){
		
		for(int jj = 0 ; jj < N ; jj++){
			sr[jj]     = srIN[jj];
			sf[jj]     = sfIN[jj];
			stress[jj] = sIN[jj];
			sbar[jj]   = sbIN[jj];
			ftt[jj]    = lsIN[jj];

		}
		for(int jj = 0 ; jj < dcat ; jj++){
			for (int kk = 0 ; kk < dlength ; kk++){
				data[jj][kk]   = dIN[jj][kk];
			}
		}
		GR = grIN;
		
		return;
	}
	
	public void setRandToClone(Random r2clone){
		
		rand = r2clone.clone();
	}
	
	public void freezeIn(Parameters params, int[] liveNbs, int[] Lives){
		
		double p = params.fget("\u03D5");
		if(p < 1){
			double rseq[] = new double[N];
			for (int jj = 0 ; jj < N ; jj++){
				rseq[jj]    = getRand().nextDouble();
				liveNbs[jj] = qN;
				Lives[jj]   = 1;
			}

			int[] rs = SortUtil.S2LindexSort(rseq);
			Ndead = (int)(p*N);
			for (int jj = 0 ; jj < Ndead ; jj++){
				sr[rs[jj]]     = getSr0();
				sf[rs[jj]]     = getSf0();
				stress[rs[jj]] = 0;
				Lives[rs[jj]]  = 0;
				failed[rs[jj]] = true;
				for (int kk = 0 ; kk < qN ; kk++) liveNbs[getNbr(rs[jj],kk)]--;
			}
		}
		else{
			double rseq[] = new double[N];
			for (int jj = 0 ; jj < N ; jj++){
				rseq[jj]    = getRand().nextDouble();
				liveNbs[jj] = qN;
				Lives[jj]   = 1;
			}

			int[] rs = SortUtil.S2LindexSort(rseq);
			Ndead = (int)(p);
			for (int jj = 0 ; jj < Ndead ; jj++){
				sr[rs[jj]]     = getSr0();
				sf[rs[jj]]     = getSf0();
				stress[rs[jj]] = 0;
				Lives[rs[jj]]  = 0;
				failed[rs[jj]] = true;
				for (int kk = 0 ; kk < qN ; kk++) liveNbs[getNbr(rs[jj],kk)]--;
			}
		}
		return;
	}
	
	public void reRandomize(int newseed, int[] Lives){
		rand = new Random(newseed);
		for(int jj = 0 ; jj < N ; jj++){
			if(Lives[jj] == 0) continue; // KEEP THE DEAD SITES THE SAME
			sr[jj]     = srn ? sr0+2*dsr*(getRand().nextDouble()-0.5) : sr0;
			sf[jj]     = sfn ? sf0+2*dsf*(getRand().nextDouble()-0.5) : sf0;
			stress[jj] = sr[jj] + (sf[jj]-sr[jj])*getRand().nextDouble();
		}
		return;
	}
	
	private int[] setupNBS(String shape){

		if(shape.equals("All Sites")){
			nbs = new LatticeNeighbors((int) L,(int) L,0,N/2,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.All);
			return null;
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
//			Job.animate();
//			// FOR DEBUGGING ONLY
		}
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
				stress[jj] += MathUtil.bool2bin(!failed[jj])*dsigma;
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
			saveData(mct, eqmode, dsigma, ndt);
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
		}
		
		fs[newindex++] = jjmax;
		failSite(jjmax,mct);		
		return;
	}

	protected void failSite(int index, int mct){
		
		failSite(index);
		return;
	}
	
	protected void failSite(int index){
		
		GR++;
		ftt[index]     = true;
		failed[index]  = true;
		return;
	}
	
	protected double nextAlpha(){
		
		return an ? a0 + 2*da*(getRand().nextDouble()-0.5) : a0;
	}
	
	protected void resetSite(int site){

		sr[site]     = nextSr(site);
		sf[site]     = nextSf(site);
		stress[site] = sr[site];
		failed[site] = false;
		return;
	}
	
	protected double nextSr(int site){
	
		return srn ? sr0+2*dsr*(getRand().nextDouble()-0.5) : sr0;
	}
	
	protected double nextSf(int site){

		return sfn ? sf0+2*dsf*(getRand().nextDouble()-0.5) : sf0;
	}
	
	public double getSW(){
		
		return (sf0+dsf-sr0-dsr);
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
	
	protected void saveData(int mct, boolean EQ, double ds, int grp){
		
		data[0][mct%dlength] = 1/Omega;
		data[1][mct%dlength] = GR;
		data[2][mct%dlength] = ds;
		data[3][mct%dlength] = MathUtil.bool2bin(EQ);
		data[4][mct%dlength] = 1.-(double)(Ndead)/(double)(N); 
		data[5][mct%dlength] = grp;
		return;
	}

	protected void configOF(){
		
		try{
			File file = new File(outdir+File.separator+bname+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print("Time");
			pw.print("\t");
			pw.print("Inverse Metric");
			pw.print("\t");
			pw.print("Shower Size");
			pw.print("\t");
			pw.print("d(sigma)");
			pw.print("\t");
			pw.print("EQ Mode");
			pw.print("\t");
			pw.print("Phi");
			pw.print("\t");
			pw.print("Unique Shower Size");
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
	
	public void PrintParams(String fout, Parameters prms, Object obj){
		
		PrintUtil.printlnToFile(fout,prms.toString());
		PrintUtil.printlnToFile(fout,"Neighbors = " + fmtI.format(qN));
		PrintUtil.printlnToFile(fout, obj.getClass().getName());
		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        Date date = new Date();
        PrintUtil.printlnToFile(fout, dateFormat.format(date));		
		
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
		configOF();
		return;
	}
	
	public void setBname(String bn, boolean bl){
		if(bl){
			setBname(bn);
			return;
		}
		bname = bn;
		return;
	}

	public Random getRand() {
	
		return rand;
	}
	
	protected int getNbr(int ffs, int nbindex){
		
		return nbs.getJ(ffs,nbSeed,nbArray,nbindex);
	}
	
	public void clearData(){
		
		if (N > dlength){
			for (int jj = 0 ; jj < N ; jj++){
				if(jj < dlength){
					for (int kk = 0 ; kk < dcat ; kk++){
						data[kk][jj] = 0;
					}
				}
				sbar[jj] = 0;
			}	
		}
		else{
			for (int jj = 0 ; jj < dlength ; jj++){
				for (int kk = 0 ; kk < dcat ; kk++){
					data[kk][jj] = 0;
				}
				if(jj < N) sbar[jj] = 0;
			}	
		}

		return;
	}
	
	public int getL(){
		
		return L;
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
	
	public boolean isAlive(int site){
		
		return (!failed[site]);
	}
	
	public double getData(int mct, int dindex){
		
		if(dindex >= ofc2Dfast.dcat) return -77;
		
		return data[dindex][mct%dlength];
	}
	
	public double getSr0(){
		
		return sr0;
	}
	
	public double getSf0(){
		
		return sf0;
	}
	
	public double getAlpha0(){
		
		return a0;
	}
	
	public int getGR(){
		
		return GR;
	}
		
	public double getSr(int st){
		
		return sr[st];
	}
	
	public double[] getSr(){
		
		return sr;
	}
	
	public double getSf(int st){
		
		return sf[st];
	}
	
	public double[] getSf(){
		
		return sf;
	}
	
	public boolean[] getLastShower(){
		
		return ftt;
	}
	
	public double[] getSbar(){
		
		return sbar;
	}
	
	public double[][] getData(){
		
		return data;
	}
	
}