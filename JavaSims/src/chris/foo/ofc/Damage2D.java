package chris.foo.ofc;

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


public class Damage2D {

	// OFC Parameters
	
	private int alive[], Nrichter;
	protected int seeds[];
	private double L, N, alpha0, alphaW, Sf0, Sr0;
	protected double Sf[];
	private double 	Sr[], SfW, SrW, stress[], R, t1, t2,
					dt, pt, Scum[], Lcum[], Omega, OmegaL, tMAX, isites;
	
	private String shape, bc;
	private boolean alphaN, SfN, SrN, killed[];
	
	protected Random rand;
	private LatticeNeighbors neighbors;
	
	private int nbs0, nbs[];

	// Damage Parameters

	private int Nlives0, LivesLeft[];
	private double NlivesW;
	private boolean NlivesN, crack, EQmode;
	private String LifeStyle;
	
	// I/O Parameters
	
	@SuppressWarnings("unused")
	private String outdir, outfile1, outfile2, PicDir;
	private boolean showGrid;
	private DecimalFormat fmtT = new DecimalFormat("000000");
	private DecimalFormat fmtH = new DecimalFormat("000");
	
	// Metric Study Parameters
	//private String outfileOmega;
	
	public Damage2D(Parameters params){
		Constructor(params);
	}

	public void Constructor(Parameters params){
		int seed;
		String animate, Mode;
		
		outdir    = params.sget("Data Directory");
		Mode      = params.sget("Mode");
		seed      = params.iget("Random Seed",0);
		shape     = params.sget("Interaction Shape");
		R         = params.fget("Interaction Radius (R)");
		L         = params.fget("Lattice Size");
		bc        = params.sget("Boundary Condtions");
		pt        = params.iget("Equilibrate Time");
		Nlives0   = params.iget("Number of Lives");
		LifeStyle = params.sget("Life Distribution");
		NlivesW   = params.fget("LD Width");
		Sf0       = params.fget("Failure Stress (\u03C3_f)");
		SfW       = params.fget("\u03C3_f width");
		Sr0       = params.fget("Residual Stress (\u03C3_r)");
		SrW       = params.fget("\u03C3_r width");
		alpha0    = params.fget("Dissipation (\u03B1)");
		alphaW    = params.fget("\u03B1 Width");
		animate   = params.sget("Animation");

		N = L*L;
		
		rand  = new Random(seed);
		
		outfile1 = outdir + File.separator+"DamageData.txt";
		outfile2 = outdir + File.separator+"DamageLives.txt";
		PicDir   = outdir + "/Pics/";
		DirUtil.MkDir(PicDir);
		
		EQmode   = (Mode.equals("Earthquake"));
		SfN      = (!(SfW == 0));
		SrN      = (!(SrW == 0));
		alphaN   = (!(alphaW == 0));
		NlivesN  = (!(EQmode) && !(LifeStyle.equals("Constant")));
		showGrid = (animate.equals("On"));

		return;
		
	}
	
	public void Initialize(){
		
		// Lattice is intact
		crack = false;
		
		// Initialize the arrays
		alive  = new int[(int)(2*N)];
		Sf     = new double[(int) N];
		Sr     = new double[(int) N];
		stress = new double[(int) N];
		Scum   = new double[(int) N];
		killed = new boolean[(int) N];
		Lcum   = new double[(int) N];	
		// Set up Lattice Neighbors
		nbs = SetUpNbs();
	
		// store initial values in arrays
		InitArrays();

		// Set up N_lives Left
		LivesLeft = SetupNlives(EQmode);
		
		// Set up First Failure
		if(pt <= 0){
			ForceFailure();
			pt = 0;
		}
		else{
			EQmode = true;
			ForceFailure();
		}
		
		// start clocks
		t1  = 0;
		t2  = 0;
		pt  = -pt;
		
	}
	
	private int[] SetUpNbs(){
		
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
	
	protected void InitArrays(){
				
		for (int jj = 0 ; jj < N ; jj++){
			
			stress[jj] = Sr0+(Sf0-Sr0)*rand.nextDouble();
			Scum[jj] = 0;
			Lcum[jj] = 0;
			
			alive[jj+(int)N] = 1;
			if(!NlivesN){ // no noise in the number of lives
				alive[jj] = Nlives0;
			}
			else{
				alive[jj] = (int) Math.round(Nlives0 + NlivesW*(rand.nextDouble()-0.5));
			}
			
			if(!SfN){	// no noise in the critical stress
				Sf[jj] = Sf0;
			}
			else{
				Sf[jj] = Sf0 + SfW*(rand.nextDouble()-0.5);
			}
			
			if(!SrN){	// no noise in the residual stress
				Sr[jj] = Sr0;
			}
			else{
				Sr[jj] = Sr0 + SrW*(rand.nextDouble()-0.5);
			}
			
		}
				
		return;
	}
	
	protected int[] SetupNlives(boolean md){
		
		
		if (!(EQmode)){
			
			int ret[];
			int max = alive[0];

			for (int jj = 0 ; jj < N ; jj++){
				if(alive[jj] > max ) max = alive[jj];
			}

			ret = new int[max + 1];

			for (int jj = 0 ; jj < max + 1 ; jj++){
				ret[jj] = 0;
			}

			for (int jj = 0 ; jj < N ; jj++){
				ret[alive[jj]]++;
			}
			return ret;
		}
		else{
			tMAX = Nlives0;
			int[] ret = {0};
			return ret;
		}
		
		
	}
	
	public boolean Equilibrate(){

		pt++;
		
		Nrichter = seeds.length; // the forced failure
		
		int[] avlnch = DistributeStress(seeds);
		ResetPlate(seeds);
		CascadeSites(avlnch);
	
		while(avlnch != null){
			int[] copy = avlnch;
			avlnch = DistributeStress(avlnch);
			ResetPlate(copy);
			CascadeSites(avlnch);
			//Job.animate();
		}
				 
		ForceFailure();
				
		return (pt < 0);
	}
	
	public boolean Avalanche(){
					
		t1++;
		t2+=dt;
		
		Nrichter = seeds.length; // the forced failure(s)
		isites = Nrichter;
		
		int[] avlnch = DistributeStress(seeds);
		ResetPlate(seeds);
		CascadeSites(avlnch);
	
		while(avlnch != null){
			int[] copy = avlnch;
			avlnch = DistributeStress(avlnch);
			ResetPlate(copy);
			CascadeSites(avlnch);
			//Job.animate();
		}
		
		Omega = StressMetric();
		
		ForceFailure();
		if(seeds == null) return false;
				
		return (EQmode) ? (t1 < tMAX) : (!crack);
	}
	
	
	protected void ForceFailure(){

		double stressMax;
		int imax = 0;

		for (int jj = 0 ; jj < N ; jj++){
			if((Sf[jj] - stress[jj]) < (Sf[imax] - stress[imax])) imax = jj;
			killed[jj] = false;
		}

		stressMax = stress[imax];		
		for (int jj = 0; jj<N; jj++){
			stress[jj] += Sf[imax]-stressMax;
		}
		dt = Sf[imax] - stressMax;
		
		seeds = CopyArray.copyArray(imax,1);
		
		ManageLives(seeds);

		return;
	}
		
	protected void ManageLives(int[] sites){

		for (int jj = 0 ; jj < sites.length ; jj++){
		
			int imax = sites[jj];

			// do only if equilibration is done and running in damage mode
			if(!(EQmode) && !(pt < 0)){	

				if(alive[imax] <= 0){
					System.out.println("All Sites Failed!");
					crack=true;
					return;
				}

				LivesLeft[alive[imax]]--;
				LivesLeft[alive[imax]-1]++;
				alive[imax]--;
			}	
			// do independent of time and mode

			alive[imax+(int)N]=0;
		
		}
		
		return;
	}
	
	private int[] DistributeStress(int[] ds){
		
		int length = ds.length;
		int toKillindex = 0;
		int[] toKill = new int[(int) N];
		int[] dsnbs;
		int Nalive = 0;
		double release;
		
		if(shape.equals("All Sites")){
			
			Nalive = (int)N - length - LivesLeft[0];
			
			if (Nalive > 0){
				for (int jj = 0 ; jj < (length - 1) ; jj++){	// loop over all dead sites except the last
					// one (the last time we record dead sites)

					if(alive[ds[jj]] > 0){
						release = (stress[ds[jj]] - Sr[ds[jj]])/Nalive;
					}
					else{
						release = stress[ds[jj]]/Nalive;
					}

					if(alphaN){	// alpha is a random variable
						for (int kk = 0 ; kk < ds[jj] ; kk++){
							stress[kk] += alive[kk+(int)N]*release*(1 - alphaW*rand.nextGaussian() - alpha0);
						}
						for (int kk = ds[jj] + 1 ; kk < N ; kk++){
							stress[kk] += alive[kk+(int)N]*release*(1 - alphaW*rand.nextGaussian() - alpha0);
						}
					}
					else{	// alpha is fixed
						for (int kk = 0 ; kk < ds[jj] ; kk++){
							stress[kk] += alive[kk+(int)N]*release*(1 - alpha0);
						}
						for (int kk = ds[jj] + 1 ; kk < N ; kk++){
							stress[kk] += alive[kk+(int)N]*release*(1- alpha0);

						}
					}
				}
				
				int jj = length - 1;	// now, the last site in the list
				
				if(alive[ds[jj]] > 0){
					release = (stress[ds[jj]] - Sr[ds[jj]])/Nalive;
				}
				else{
					release = stress[ds[jj]]/Nalive;
				}

				if(alphaN){	// alpha is a random variable
					for (int kk = 0 ; kk < ds[jj] ; kk++){
						stress[kk] += alive[kk+(int)N]*release*(1 - alphaW*rand.nextGaussian() - alpha0);
						if (alive[kk+(int)N]*stress[kk] > Sf[kk]) toKill[toKillindex++] = kk;						
					}
					for (int kk = ds[jj] + 1 ; kk < N ; kk++){
						stress[kk] += alive[kk+(int)N]*release*(1 - alphaW*rand.nextGaussian() - alpha0);
						if (alive[kk+(int)N]*stress[kk] > Sf[kk]) toKill[toKillindex++] = kk;
					}
				}
				else{	// alpha is fixed
					for (int kk = 0 ; kk < ds[jj] ; kk++){
						stress[kk] += alive[kk+(int)N]*release*(1 - alpha0);
						if (alive[kk+(int)N]*stress[kk] > Sf[kk]) toKill[toKillindex++] = kk;
					}
					for (int kk = ds[jj] + 1 ; kk < N ; kk++){
						stress[kk] += alive[kk+(int)N]*release*(1- alpha0);
						if (alive[kk+(int)N]*stress[kk] > Sf[kk]) toKill[toKillindex++] = kk;
					}
				}
			}
		}
		else{	// loop over a subset of the lattice sites
			
			for (int jj = 0 ; jj < length ; jj++){

				Nalive = 0;
				
				dsnbs = neighbors.get(ds[jj],nbs0, nbs);
				//dsnbs = neighbors.get(ds[jj]);

				
				for (int kk = 0 ; kk < dsnbs.length ; kk++){
					Nalive += alive[dsnbs[kk]+(int)N];
				}
				
				if(Nalive == 0) continue;
				
				if(alive[ds[jj]] > 0){
					release = (stress[ds[jj]] - Sr[ds[jj]])/Nalive;
				}
				else{
					release = stress[ds[jj]]/Nalive;
				}
						
				if(alphaN){	// alpha is a random variable
					
					for (int kk = 0 ; kk < dsnbs.length ; kk++){
						stress[dsnbs[kk]] += alive[dsnbs[kk]+(int)N]*release*(1 - alphaW*rand.nextGaussian() - alpha0);	
						if ((alive[dsnbs[kk]+(int)N]*stress[dsnbs[kk]] > Sf[dsnbs[kk]]) && !(killed[dsnbs[kk]])){
							toKill[toKillindex++] = dsnbs[kk];
							killed[dsnbs[kk]] = true;
						}
					}
				}
				else{	// alpha is fixed
					for (int kk = 0 ; kk < dsnbs.length ; kk++){
						stress[dsnbs[kk]] += alive[dsnbs[kk]+(int)N]*release*(1 - alpha0);
						if ((alive[dsnbs[kk]+(int)N]*stress[dsnbs[kk]] > Sf[dsnbs[kk]]) && !(killed[dsnbs[kk]])){
							toKill[toKillindex++] = dsnbs[kk];
							killed[dsnbs[kk]] = true;
						}
					}
				}
			}
			
			
		}

		return CopyArray.copyArray(toKill, toKillindex);
	}
	
	protected void ResetPlate(int[] ds){
		
		int length = ds.length;
		
		if(EQmode){
			if(SrN){	// annealed noise on Sr
				for (int jj = 0 ; jj < length; jj++){
					stress[ds[jj]] = Sr[ds[jj]];
					Sr[ds[jj]] = nextSr(ds[jj]);
					alive[ds[jj]+(int)N] = 1;
				}
			}
			else{
				for (int jj = 0 ; jj < length; jj++){
					stress[ds[jj]] = Sr[ds[jj]];
					alive[ds[jj]+(int)N] = 1;
	
				}
			}
		}
		else{	// damage mode
			if(SrN){	// annealed noise on Sr
				for (int jj = 0 ; jj < length; jj++){
					if(alive[ds[jj]] > 0){
						stress[ds[jj]] = Sr[ds[jj]];
						Sr[ds[jj]] = nextSr(ds[jj]);
						alive[ds[jj]+(int)N] = 1;
					}
					else{
						stress[ds[jj]] = -2;
					}
				}
			}
			else{
				for (int jj = 0 ; jj < length; jj++){
					if(alive[ds[jj]] > 0){
						stress[ds[jj]] = Sr[ds[jj]];
						alive[ds[jj]+(int)N] = 1;
					}
					else{
						stress[ds[jj]] = -2;
					}
				}
			}
		}
		
		
		return;
	}
	
	private void CascadeSites(int[] ds){
		
		if(ds == null) return;
		
		adjustNF(ds);
		
		int length = ds.length;
		
		if(EQmode){
			for (int jj = 0 ; jj < length ; jj++){
				alive[ds[jj]+(int)N] = 0;
				killed[ds[jj]] = false;
			}
		}
		else{
			for (int jj = 0 ; jj < length ; jj++){
				alive[ds[jj]+(int)N] = 0;
				LivesLeft[alive[ds[jj]]]--;
				if(alive[ds[jj]] > 0) LivesLeft[alive[ds[jj]]-1]++;
				alive[ds[jj]]--;
				killed[ds[jj]] = false;
			}
		}
		
		Nrichter+=length;
		return;
		
	}
	
	protected double nextSr(int site){
		
		return (Sr0 + SrW*(rand.nextDouble()-0.5));
	}
	
	private double StressMetric(){
		
		if(N > LivesLeft[0]){
			double Sbar = 0;

			for (int jj = 0 ; jj < N ; jj++){
				Scum[jj] += alive[jj+(int)N]*stress[jj]*dt;
			}

			for (int jj = 0 ; jj < N ; jj++){
				Sbar += Scum[jj];
			}

			Sbar = Sbar/((int)N - LivesLeft[0]); // all the sites save the dead ones

			Omega = 0;
			for (int jj = 0 ; jj < N ; jj++){
				Omega += (Scum[jj] - Sbar)*(Scum[jj] - Sbar);
			}

			Omega = Omega/(t2*t2*((int)N - LivesLeft[0]));
			
			
			//PrintUtil.printlnToFile(outfile2,t2,Sbar/t2,Scum[1]/t2,Scum[(int)(getN()/2)]/t2);	// FOR TESTING ONLY!!!!!!
			
			/*
			 * 
			 * NB	For ForceSf this LivesLeft[0] is probably not a denominator 
			 * 		since sites that fail once are, in principle, dead. Maybe run in 
			 * 		damage mode with Nlives = 1?
			 * 
			 */
			
			return Omega;
			
		}
		
		return -1;
	}
	
	public void LifeMetric(){
		
		OmegaL = 0;
		double Lbar = 0;
		
		if(N > LivesLeft[0]){

			for (int jj = 0 ; jj < N ; jj++){
				Lcum[jj] += alive[jj]*dt;
				Lbar += Lcum[jj];
			}

	

			Lbar = Lbar/((int)N - LivesLeft[0]); // all the sites save the dead ones


			for (int jj = 0 ; jj < N ; jj++){
				OmegaL += (Lcum[jj] - Lbar)*(Lcum[jj] - Lbar);
			}

			OmegaL = OmegaL/(t2*t2*((int)N - LivesLeft[0]));
			
			
			//PrintUtil.printlnToFile(outfile2,t2,Sbar/t2,Scum[1]/t2,Scum[(int)(getN()/2)]/t2);	// FOR TESTING ONLY!!!!!!
			
			/*
			 * 
			 * NB	For ForceSf this LivesLeft[0] is probably not a denominator 
			 * 		since sites that fail once are, in principle, dead. Maybe run in 
			 * 		damage mode with Nlives = 1?
			 * 
			 */
		}
		
			return;
		
		
	}
	
	public void WriteDataHeaders(){
		
		PrintUtil.printlnToFile(outfile1,"Sweeps","t_vel","N_sts/avl","Omega","Omega_L","Init Sites","N_dead");
		//PrintUtil.printlnToFile(outfile2,"Sweeps","t_vel","N_lvs=0","N_lvs=1","N_lvs=N_max/2","N_lvs=N_max-1","N_lvs=N_max");
		;
		
	}
	
	public void TakeDate(){
		PrintUtil.printlnToFile(outfile1,t1,t2,Nrichter,Omega,OmegaL,isites, LivesLeft[0]);
		//PrintUtil.printlnToFile(outfile2,t1,t2,LivesLeft[0],LivesLeft[1],LivesLeft[(int)(Nlives0/2)],LivesLeft[Nlives0-1],LivesLeft[Nlives0]);
	}
	
	public static void PrintParams(String fout, Parameters prms){
		
		PrintUtil.printlnToFile(fout,prms.toString());
		return;
	}
	
	public void TakePicture(Grid grid){
		
		if(showGrid){
		
			String SaveAs = PicDir + File.separator + grid.getTitle()+fmtT.format(t1)+"-"+fmtH.format(Nrichter)+".png";

			try {
				ImageIO.write(grid.getImage(), "png", new File(SaveAs));
			} catch (IOException e) {
				System.err.println("Error in Writing File" + SaveAs);
			}
		
		}
		
		return;
	}
	
	public double[] getStress(){
		
		return stress;
	}
	
	public int maxLives(){
		
		return LivesLeft.length - 1;
	}
	
	public int[] getLives(){
		return CopyArray.copyArray(alive,(int)N);
	}
	
	public double getTime(int t){
		
		if(t==1){
			return t1;
		}
		else if(t==2){
			return t2;
		}
		else if(t==-1){
			return pt;
		}
		else{
			return Integer.MIN_VALUE;
		}
		
	}
	
	public boolean draw(){
		
		return showGrid; 
	}
	
	public String getOutdir(){
		
		return outdir;
	}
	
	public double getSmax(){
		
		return Sf0 + 5*SfW;
	}
	
	public int getN(){
		
		return (int) N;
	}
	
	public int getL(){
		
		return (int) L;
	}
	
	public int getAS(){
		
		return Nrichter;
	}
	
	public double getAveLL(){

		double ret = 0;
		
		for(int jj = 0 ; jj < LivesLeft.length ; jj++){
			ret += jj*LivesLeft[jj];
		}
		
		return ret/getN();
	}
	
	public void resetMode(Parameters params){
		
		EQmode = (params.sget("Mode").equals("Earthquake"));		
		
		return;
	}
	
	public boolean getEQmode(){
		
		return EQmode;
	}
	
	public double getSr0(){
		
		return Sr0;
	}
	
	public double getSrW(){
		
		return SrW;
	}
	
	public double getSf0(){
		
		return Sf0;
	}
	
	public double getStress(int site){
		
		return (site < N) ? stress[site] : Integer.MIN_VALUE; 
		
	}
	
	public double getSr(int site){
		
		return (site < N) ? Sr[site] : Integer.MIN_VALUE;
	}
	
	protected void setT1T2(){
		
		t2 = t1;
		dt = 1;
	}
	
	protected void setStress(int site, double str){
		
		if(site < N){
			stress[site] = str;
		}
		
		return;
	}
	
	
	public void setSf(int site, double val){
		
		if(site < N){
			Sf[site] = val;
		}
		
		return;
	}
	
	protected void adjustNF(int[] sites){
		return;
	}
	
	protected void killSite(int site){
		
		if(site < N) alive[site] = 0;
		return;
	}
	
	public void setOutfile(String str){
		
		outfile1 = outdir + File.separator + str + ".txt";
		
		return;
	}
}
