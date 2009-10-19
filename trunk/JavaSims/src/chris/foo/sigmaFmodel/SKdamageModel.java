package chris.foo.sigmaFmodel;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Random;

import javax.imageio.ImageIO;

import scikit.graphics.dim2.Grid;
import scikit.jobs.params.Parameters;
import chris.util.CopyUtil;
import chris.util.DirUtil;
import chris.util.LatticeNeighbors;
import chris.util.PrintUtil;

public class SKdamageModel {	
	
	// Model Parameters
	private double  R, L, N, Sr0, Sr[], SrW, dSr, dSrW, Sf0, Sf[], SfW, alpha0, alphaW,
					stress[], Scum[][], t1, t2, t3, dt2, dt3, Nrichter, numFailures[],
					globalFailures[], Omega1, Omega2, Omega3, sbarT, nfbar, dsRMS, sT,
					nDS;

	private boolean SrN, dSrN, SfN, alphaN, deadSite[], ks[];
	private Random  rand;
	private String  shape, bc;
	private int     alive[], seedsites[], nbs0, nbs[], liveSites[];
	private LatticeNeighbors neighbors;
	
	// Data Collection Parameters
	private String outdir, datafile1, datafile2, picdir;
	
	// I/O Parameters
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
		Sf0    = params.fget("Failure Stress (\u03C3_f)");
		SfW    = params.fget("\u03C3_f width");
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
		sbarT  = 0;
		nDS    = 0;
		
		picdir    = outdir +"/GridPics/";
		datafile1 = outdir + File.separator + "StressData.txt"; 
		datafile2 = outdir + File.separator + "MetricData.txt"; 
		DirUtil.MkDir(picdir);
	
		InitializeArrays();
		
		return;
	}
	
	private void InitializeArrays(){
		
		Sr             = new double[intN()];
		Sf             = new double[intN()];
		stress         = new double[intN()];
		Scum           = new double[3][intN()];
		alive          = new int[intN()];
		numFailures    = new double[intN()];
		deadSite       = new boolean[intN()];
		ks             = new boolean[intN()];
		globalFailures = new double[intN()];
		liveSites	   = new int[intN()];
		
		for (int jj = 0 ; jj < intN() ; jj++){
			
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
//			else{
//				
//			}
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
		if (sysFail()) return false;
		
		Nrichter++; // the forced failure
		
		while((seedsites = distributeStress(seedsites)) != null);
		
		sbarT = ergMetric();
		
		return true;
	}
	
	private int[] forceFailure(){
		
		int trialsite;
		int length;
		double trialMove;
		
		dt2   = 0;
		dt3   = 0;
		dsRMS = 0;

		if ((length = liveSites.length) == 0) return null;
		
	
		while(true){
			trialsite = liveSites[rand.nextInt(length)];
			trialMove = Sf[trialsite] - (Sr[trialsite] + (Sf[trialsite] - Sr[trialsite])*rand.nextDouble());
			dt3 += trialMove;
			Sf[trialsite] -= trialMove;
			dt2++;
			if(stress[trialsite] > Sf[trialsite]) break;	
		}
		t1++;
		t2 += dt2;
		t3+= dt3;
		
		return CopyUtil.copyArray(trialsite,1);
	}
	
	private boolean sysFail(){
		
		return (seedsites == null);
	}
	
	private int[] distributeStress(int[] sites){
		
		int Nalive, dsnbs[], newlykilled[], counter;
		double release;
		
		newlykilled = new int[intN()];
		counter     = 0;
		
		if (liveSites.length == 0) return null;
		
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
		
		return CopyUtil.copyArray(newlykilled,counter);
	}
	
	private void resetSites(int[] sites){
		
		for (int jj = 0 ; jj < sites.length ; jj++){
			int st = sites[jj];
			globalFailures[(int) numFailures[st]++]--;	
			globalFailures[(int) numFailures[st]]++;
			if(Sr[st] <= 0){
				killSite(st);
			}
			else{
				alive[st]  = 1;
				Sf[st]     = nextSf(st); // ASK / THINK ABOUT ORDER
				Sr[st]     = nextSr(st);
				ks[st]     = false;
				stress[st] = Sr[st];
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
		alive[site]       = 0;
		ks[site]          = true;
		deadSite[site]    = true;
		numFailures[site] = -1;
		stress[site]      = 0;
		
		int Nlength = liveSites.length - 1;
		
		int[] foo   = new int[Nlength];
		int counter = 0;
		
		for (int jj = 0 ; jj < (Nlength + 1) ; jj++){
			if(liveSites[jj] == site) continue;
			foo[counter++] = liveSites[jj];
		}

		liveSites = foo;		
		return;
	}
	
	private double nextalpha(){

		return(alphaN) ? alpha0 + alphaW*(0.5 - rand.nextDouble()) : alpha0;
	}
	
	private double nextSf(int site){
		
		return (Sr[site] + (Sf0 - Sr[site])*rand.nextDouble());
	}
	
	private double nextSr(int site){
		
		double val = Sr[site] - dSr;
		
		if(dSrN)     val += dSrW*(0.5 - rand.nextDouble());
		if(val < 0) val = 0;	
			
		return val;
	}
	
	public int intN(){
		return (int) N;
	}
	
	public String getOutdir(){
		
		return outdir;
	}
	
	public void writeDataHeaders(){
		
		PrintUtil.printlnToFile(datafile1,"t1","t2","t3","N_Sts/avlnch","<s-s_f>_rms","<s>(t)","<s(t)>","<failures>");
		PrintUtil.printlnToFile(datafile2,"t1","t2","t3", "Omega1","Omega2","Omega3","N_dead");
		return;
	}
	
	public void takeData(){
		
		PrintUtil.printlnToFile(datafile1,t1, t2, t3, Nrichter, Math.sqrt(dsRMS/(Nrichter - 1)), sbarT, sT, nfbar);
		PrintUtil.printlnToFile(datafile2,t1, t2, t3, Omega1, Omega2, Omega3, nDS);
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
	
	public double getAvlnchSize(){
		
		return Nrichter;
	}
	
	private double ergMetric(){
		
		double[] aTime = new double[]{t1, t2, t3};
		double[] adt = new double[]{1, dt2, dt3};
		double Sbar[] = new double[]{0, 0, 0};
		double[] Omegas = new double[]{0, 0, 0};
		
		nfbar = 0;
		
		int aliveSites = liveSites.length;
		
		if(aliveSites == 0) return -1;
		
		for (int jj = 0 ; jj < intN() ; jj++){
			for (int kk = 0 ; kk < 3 ; kk++){
				Scum[kk][jj] += alive[jj]*stress[jj]*adt[kk];
				Sbar[kk] += Scum[kk][jj];
			}
			nfbar += jj*globalFailures[jj];
		}
		nfbar = nfbar / N;
		
		for (int kk = 0 ; kk < 3 ; kk++){
			Sbar[kk] = (Sbar[kk] / aliveSites);
		}
		
		for (int jj = 0 ; jj < intN() ; jj++){
			for (int kk = 0 ; kk < 3 ; kk++){
				Omegas[kk] += (Scum[kk][jj] - Sbar[kk])*(Scum[kk][jj] - Sbar[kk]);
			}
		}
		
		for (int kk = 0 ; kk < 3 ; kk++){
			Omegas[kk] = Omegas[kk] / (aTime[kk]*aTime[kk]*aliveSites);
		}
		
		Omega1 = Omegas[0];
		Omega2 = Omegas[1];
		Omega3 = Omegas[2];
		
		sT = Scum[0][(int)(N/3)]/aTime[0];
		return (Sbar[0]/aTime[0]);
	}
	
	public double getSbar(){
		
		return sbarT;
	}
	
	public double getNdead(){
		
		return nDS;
	}
	
}
