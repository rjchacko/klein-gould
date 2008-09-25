package chris.ofcdamage;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Random;
import javax.imageio.ImageIO;
import scikit.graphics.dim2.Grid;
import scikit.jobs.params.Parameters;
import chris.util.CopyArray;
import chris.util.DirUtil;
import chris.util.LatticeNeighbors;
import chris.util.PrintUtil;

public class damage {	
	
	// Model Parameters
	private double  R, L, N, Sr0, Sr[], SrW, Sf0, Sf[], SfW, alpha0, alphaW,
					stress[], Scum[], t1, t2, dt2, Nrichter,
					Omega, sbarT, dsRMS, sT, nDS;
	private boolean SrN, SfN, alphaN, deadSite[], ks[], equil, clocks, randSi;
	private Random  rand;
	private String  shape, bc;
	private int     alive[], nbs0, nbs[], nlMin, nlMax, numLives[];
	private LatticeNeighbors neighbors;
	
	// Data Collection Parameters
	private String outdir, datafile1, picdir;
	
	// I/O Parameters
	private DecimalFormat fmtI = new DecimalFormat("00000000");
	
	
	
	public damage(Parameters params){
	
		SKdamageConstr(params);
		return;
	}
	
	public void SKdamageConstr(Parameters params){
	
		int seed;
		
		outdir = params.sget("Data Directory");
		seed   = params.iget("Random Seed");
		shape  = params.sget("Interaction Shape");
		R      = params.fget("Interaction Radius (R)");
		L      = params.iget("Lattice Size");
		nlMin  = params.iget("Min Lives");	
		nlMax  = params.iget("Max Lives");	
		bc     = params.sget("Boundary Condtions");
		Sf0    = params.fget("Failure Stress (\u03C3_f)");
		SfW    = params.fget("\u03C3_f width");
		Sr0    = params.fget("Residual Stress (\u03C3_r)");
		SrW    = params.fget("\u03C3_r width");
		alpha0 = params.fget("Dissipation (\u03B1)");
		alphaW = params.fget("\u03B1 width");
		randSi = params.sget("Intitial Stess").equals("Random");

		
		rand = new Random(seed);

		t1     = 0;
		t2     = 0;
		N      = L*L;
		SrN    = (SrW > 0);
		SfN    = (SfW > 0);
		alphaN = (alphaW > 0);
		sbarT  = 0;
		nDS    = 0;
		equil  = true;
		clocks = false;
		
		picdir    = outdir +"/GridPics/";
		datafile1 = outdir + File.separator + "StressData.txt"; 
		DirUtil.MkDir(picdir);
	
		InitializeArrays();
		
		return;
	}
	
	private void InitializeArrays(){
		
		Sr             = new double[intN()];
		Sf             = new double[intN()];
		stress         = new double[intN()];
		Scum           = new double[intN()];
		alive          = new int[intN()];
		numLives       = new int[intN()];
		deadSite       = new boolean[intN()];
		ks             = new boolean[intN()];
		
		for (int jj = 0 ; jj < intN() ; jj++){
			
			alive[jj]          = 1;
			deadSite[jj]       = false;
			ks[jj]             = false;
			numLives[jj]       = nlMin + rand.nextInt(nlMax - nlMin + 1);
			
			if(randSi){
				stress[jj] = Sr0 + (Sf0 - Sr0)*rand.nextDouble();
			}
			else{
				stress[jj] = Sr0 + SrW;
			}
			if(SrN){
				Sr[jj] = Sr0 + SrW*(0.5 - rand.nextDouble());
			}
			else{
				Sr[jj] = Sr0;
			}
			if(SfN){
				Sf[jj] = Sf0 + SfW*(0.5 - rand.nextDouble());
			}
			else{
				Sf[jj] = Sf0;
			}
	
		}

		
		nbs = setupNBS();
		
		return;
	}
	
	private int[] setupNBS(){
		
		if(shape.equals("All Sites")){
			neighbors = new LatticeNeighbors((int) L,(int) L,0,N/2,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.All);
			return neighbors.get(nbs0);	
		}
		
		if (bc.equals("Bordered")){
			nbs0 = (int)((1 + L)*(L/2));
			if(shape.equals("Circle")){
				neighbors = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Circle);
			}
			else if(shape.equals("Square")){
				neighbors = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Square);
			}
			else if(shape.equals("Diamond")){
				neighbors = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Diamond);
			}
		}
		else{
			nbs0 = 0;
			if(shape.equals("Circle")){
				neighbors = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Circle);
			}
			else if(shape.equals("Square")){
				neighbors = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Square);
			}
			else if(shape.equals("Diamond")){
				neighbors = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Diamond);
			}
		}

		return neighbors.get(nbs0);	
	}
	
	
	public void equilibrate(){
		
		Nrichter = 0;
		
		int[] seedsites = forceFailure();
		
		Nrichter++; // the forced failure
		
		while((seedsites = distributeStress(seedsites)) != null);
		
		return;
	}
	
	public boolean avalanche(){

		Nrichter = 0;
		
		int[] seedsites = forceFailure();
		if (sysFail()) return false;
		
		Nrichter++; // the forced failure
		
		while((seedsites = distributeStress(seedsites)) != null);
		
		ergMetric();
		
		return true;
	}
	
	private int[] forceFailure(){

		double stressMax;
		int imax = 0;
		
		dt2   = 0;
		dsRMS = 0;

		for (int jj = 0 ; jj < N ; jj++){
			if((Sf[jj] - stress[jj]) < (Sf[imax] - stress[imax])) imax = jj;
		}

		stressMax = stress[imax];		
		for (int jj = 0; jj<N; jj++){
			stress[jj] += Sf[imax]-stressMax;
		}
		
		if(clocks){	
			dt2 = Sf[imax] - stressMax;
			t1++;
			t2 += dt2;
		}

		return CopyArray.copyArray(imax,1);
	}
	
	
	protected boolean sysFail(){
		
		return (nDS == N);
	}
	
	private int[] distributeStress(int[] sites){
		
		int Nalive, dsnbs[], newlykilled[], counter;
		double release;
		
		newlykilled = new int[intN()];
		counter     = 0;
		
		if (nDS == N) return null;
		
		failSites(sites);
		
		for (int jj = 0 ; jj < sites.length ; jj++){
		
			int st = sites[jj]; 
			Nalive = 0;
			dsnbs = neighbors.get(st, nbs0, nbs);
		
			for (int kk = 0 ; kk < dsnbs.length ; kk++){
				Nalive += alive[dsnbs[kk]];
			}
			
			if(Nalive == 0) continue;
			release = (stress[st] - Sr[st])/Nalive;
			
			for (int kk = 0 ; kk < dsnbs.length ; kk++){

				stress[dsnbs[kk]] += alive[dsnbs[kk]]*release*(1-nextalpha());
				if( !(ks[dsnbs[kk]]) && (alive[dsnbs[kk]]*stress[dsnbs[kk]] > Sf[dsnbs[kk]]) ){
					ks[dsnbs[kk]] = true;
					newlykilled[counter++] = dsnbs[kk];
					dsRMS += (stress[dsnbs[kk]] - Sf[dsnbs[kk]])*(stress[dsnbs[kk]] - Sf[dsnbs[kk]]);
				}
			}
		}
		
		Nrichter += counter;
		resetSites(sites);
		
		return CopyArray.copyArray(newlykilled,counter);
	}
	
	private void resetSites(int[] sites){
		
		for (int jj = 0 ; jj < sites.length ; jj++){
			int st = sites[jj];
			
			if(equil){
				alive[st]  = 1;
				Sf[st]     = nextSf(st);
				Sr[st]     = nextSr(st);
				ks[st]     = false;
				stress[st] = Sr[st];
			
			}
			else{
				if(numLives[st] == 1){
					killSite(st);
				}
				else{
					alive[st]  = 1;
					Sf[st]     = nextSf(st);
					Sr[st]     = nextSr(st);
					ks[st]     = false;
					stress[st] = Sr[st];
					numLives[st]--;
				}
			}
		}
		
		return;
	}
	
	private void failSites(int[] sites){
		
		for (int jj = 0 ; jj < sites.length ; jj++){
			alive[sites[jj]] = 0;
		}
		
		return;
	}
	
	private void killSite(int site){
				
		nDS++;
		numLives[site]--;
		alive[site]       = 0;
		ks[site]          = true;
		deadSite[site]    = true;
		stress[site]      = 0;
		
		add2cluster(site);
	
		return;
	}
	
	protected void add2cluster(int site){
		
		return;
	}
	
	private double nextalpha(){

		return(alphaN) ? alpha0 + alphaW*(0.5 - rand.nextDouble()) : alpha0;
	}
	
	private double nextSf(int site){
		
		return (SfN) ? Sf0 + SfW*(0.5 - rand.nextDouble()) : Sf0;
	}
	
	private double nextSr(int site){
			
		return (SrN) ? Sr0 + SrW*(0.5 - rand.nextDouble()) : Sr0;
	}
	
	public int intN(){
		return (int) N;
	}
	
	public String getOutdir(){
		
		return outdir;
	}
	
	public void writeDataHeaders(){
		
		PrintUtil.printlnToFile(datafile1,"t1","t2","N_Sts/avlnch","<s-s_f>_rms","<s>(t)","<s(t)>","N_dead", "Omega");
		return;
	}
	
	public void takeData(){
		
		PrintUtil.printlnToFile(datafile1,t1, t2, Nrichter, Math.sqrt(dsRMS/(Nrichter - 1)), sbarT, sT, nDS, Omega);
		return;
	}
	
	public static void printParams(String fout, Parameters prms){
		
		PrintUtil.printlnToFile(fout,prms.toString());
		return;
	}
	
	public void takePicture(Grid grid, boolean display){
		
		if(!display) return;

		String SaveAs = picdir + File.separator + grid.getTitle()+fmtI.format(t1)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
		return;
	}
	
	public void takePicture(Grid grid, boolean display, String fn){
		
		if(!display) return;

		try {
			ImageIO.write(grid.getImage(), "png", new File(fn));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + fn);
		}
		
		return;
	}
	
	public int getL(){
		
		return (int) L;
	}
	
	public double[] getStress(){
		
		return stress;
	}
	
	public double getTime(int id){

		if(id == 1){
			return t1;
		}
		else if (id == 2){
			return t2;
		}
		else{
			return -77;
		}
		
	}
	
	public double getAvlnchSize(){
		
		return Nrichter;
	}
	
	private void ergMetric(){
		
		double Sbar = 0;
		int count = 0;
		
		if(nDS < N){
			
			for (int jj = 0 ; jj < N ; jj++){
				Scum[jj] += alive[jj]*stress[jj]*dt2;
				count += alive[jj];
				Sbar += alive[jj]*Scum[jj];
			}
			
			Sbar = Sbar/count;

			Omega = 0;
			for (int jj = 0 ; jj < N ; jj++){
				Omega += alive[jj]*(Scum[jj] - Sbar)*(Scum[jj] - Sbar);
//				PrintUtil.printlnToFile("/Users/cserino/Desktop/debug.txt",jj,alive[jj],Scum[jj],Sbar);
			}
	
			Omega = Omega/(t2*t2*count);
			sbarT = Sbar / t2;
			sT    = Scum[(int)(N/3)] / t2; 
			
		}
			
		return;
	}
	
	public double getSbar(){
		
		return sbarT;
	}
	
	public double getNdead(){
		
		return nDS;
	}
	
	public void setEquil(boolean bool){
		
		equil = bool;
		return;
	}
	
	public int maxLives(){
		
		return nlMax;
	}
	
	public int[] getLives(){
		
		return numLives;
	}
	
	public void runClocks(boolean bool){
		clocks = bool;
	}
	
	public void printStress(String fn){
		PrintUtil.printArrayToFile(fn, stress, (int) L, (int) L);
		return;
	}
	
	public void printSr(String fn){
		PrintUtil.printArrayToFile(fn, Sr, (int) L, 1);
		return;
	}
	
	public String getBC(){
		
		return bc;
	}
	
	public void setDataFile(int str){
		
		datafile1 = outdir + File.separator + "StressData_" + str + ".txt";
	}
	
}
