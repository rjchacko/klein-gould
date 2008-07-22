package chris.sigmaFmodel;

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

public class SKdamageModel {	
	
	// Model Parameters
	@SuppressWarnings("unused")
	private double  R, L, N, Sr0, Sr[], SrW, dSr, dSrW, Sf0, Sf[], SfW, alpha0, alphaW,
					stress[], Scum[][], t1, t2, t3, Nrichter, numFailures[], globalFailures[];
	@SuppressWarnings("unused")
	private boolean SrN, dSrN, SfN, alphaN, deadSite[], ks[];
	private Random  rand;
	private String  shape, bc;
	private int     alive[], seedsites[], nbs0, nbs[], liveSites[];
	private LatticeNeighbors neighbors;
	
	// Data Collection Parmeters
	@SuppressWarnings("unused")
	private String outdir, datafile1;
	
	// I/O Parameters
	private String picdir;
	private DecimalFormat fmtI = new DecimalFormat("00000000");
	
	
	
	public SKdamageModel(Parameters params){
	
		SKdamageConstr(params);
		return;
	}
	
	public void SKdamageConstr(Parameters params){
	
		int seed;
		
		outdir = params.sget("Data Directory");
		seed   = params.iget("Seed");
		shape  = params.sget("Interaction Shape");
		R      = params.fget("Interaction Radius (R)");
		L      = params.iget("Lattice Size (L)");
		bc     = params.sget("Boundary Condtions");
		Sr0    = params.fget("Residual Stress (\u03C3_r)");
		SrW    = params.fget("\u03C3_r width");
		dSr    = params.fget("d(\u03C3_r)");
		dSrW   = params.fget("d(\u03C3_r) width");
		alpha0  = params.fget("Dissipation (\u03B1)");
		alphaW = params.fget("\u03B1 width");

		rand = new Random(seed);

		t1     = 0;
		t2     = 0;
		t3     = 0;
		N      = L*L;
		SrN    = (SrW > 0);
		SfN    = (SfW > 0);
		alphaN = (alphaW > 0);

		picdir    = outdir +"/GridPics/";
		datafile1 = outdir + File.separator + "Data.txt"; 
		DirUtil.MkDir(picdir);
	
		InitializeArrays();
		
		return;
	}
	
	private void InitializeArrays(){
		
		Sr             = new double[Nint()];
		Sf             = new double[Nint()];
		stress         = new double[Nint()];
		Scum           = new double[3][Nint()];
		alive          = new int[Nint()];
		numFailures    = new double[Nint()];
		deadSite       = new boolean[Nint()];
		ks             = new boolean[Nint()];
		globalFailures = new double[Nint()];
		liveSites	   = new int[Nint()];
		
		for (int jj = 0 ; jj < Nint() ; jj++){
			
			liveSites[jj]      = jj;
			alive[jj]          = 1;
			numFailures[jj]    = 0;
			deadSite[jj]       = false;
			ks[jj]             = false;
			globalFailures[jj] = 0;
			
			for (int kk = 0; kk < 3 ; kk++){
				Scum[kk][jj] = 0;
			}
			
			if(true){
				stress[jj] = Sr0 + (Sf0 - Sr0)*rand.nextDouble();
			}
			else{
				
			}
			if(SrN){
				Sr[jj] = Sr0 + SrW*rand.nextDouble();
			}
			else{
				Sr[jj] = Sr0;
			}
			if(SfN){
				Sf[jj] = Sf0 + SfW*rand.nextDouble();
			}
			else{
				Sf[jj] = Sf0;
			}
		}
		globalFailures[0] = N;
		
		nbs = setupNBS();
		
		return;
	}
	
	private int[] setupNBS(){
		
		if(shape.equals("All Sites")) return null;
		
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
	
	public boolean avalanche(){

		Nrichter = 0;
		
		seedsites = forceFailure();
		if (seedsites == null) return false;
		
		Nrichter++; // the forced failure

		while((seedsites = distributeStress(seedsites)) != null){
			System.out.println("stuck");
		}
		
		return true;
	}
	
	private int[] forceFailure(){
		
		int trialsite;
		int length;
		double dt3;
		
		if ((length = liveSites.length) == 0) return null;
		
	
		while(true){
			trialsite = liveSites[rand.nextInt(length)];
			dt3 = Sf[trialsite] - (Sr[trialsite] + (Sf[trialsite] - Sr[trialsite])*rand.nextDouble());
			t3 += dt3;
			Sf[trialsite] -= dt3;
			t2++;
			if(stress[trialsite] > Sf[trialsite]) break;	
		}
		t1++;
		
		return CopyArray.copyArray(trialsite,1);
	}
	
	private int[] distributeStress(int[] sites){
		
		int Nalive, dsnbs[], newlykilled[], counter;
		double release;
		
		newlykilled = new int[Nint()];
		counter     = 0;
		
		if (liveSites.length == 0) return null;
		
		failSites(sites);
		
		for (int jj = 1 ; jj < sites.length ; jj++){
		
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
				if( !(ks[dsnbs[kk]]) && (alive[dsnbs[kk]]*stress[dsnbs[kk]] > Sr[dsnbs[kk]]) ){
					ks[dsnbs[kk]] = true;
					newlykilled[counter++] = dsnbs[kk];
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
			globalFailures[(int) numFailures[st]++]--;	
			globalFailures[(int) numFailures[st]]++;
			if(Sr[st] == 0){
				killSite(st);
			}
			else{
				alive[st] = 1;
				Sf[st]    = nextSf(st);
				Sr[st]    = nextSr(st);
				ks[st]    = false;
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
		
		alive[site]       = 0;
		ks[site]          = true;
		deadSite[site]    = true;
		numFailures[site] = -1;
		stress[site]      = 0;
		
		int[] foo = new int[liveSites.length - 1];
		for (int jj = 0 ; jj < site ; jj++){
			foo[jj] = liveSites[jj];
		}
		for (int jj = site+1 ; jj < Nint() ; jj++){
			foo[jj] = liveSites[jj];
		}
		liveSites = foo;
		
		return;
	}
	
	private double nextalpha(){
		
		return(alphaN) ? alpha0 + alphaW*rand.nextDouble() : alpha0;
	}
	
	private double nextSf(int site){
		
		return Sf0;
	}
	
	private double nextSr(int site){
		
		return Sr0;
	}
	
	public int Nint(){
		return (int) N;
	}
	
	public String getOutdir(){
		
		return outdir;
	}
	
	public void writeDataHeaders(){
	
		return;
	}
	
	public void takeData(){
		
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
	
	public int getL(){
		
		return (int) L;
	}
	
	public double[] getFailures(){
		
		return numFailures;
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
		else if (id == 3){
			return t3;
		}
		else{
			return -77;
		}
		
	}
	
	
}
