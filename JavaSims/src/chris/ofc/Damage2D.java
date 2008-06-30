package chris.ofc;

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



/*
 * 
 * 
 * Make EQ an option and thus merge, for example
 * Avalanche and AvalancheEQ, or what is more important
 * maybe, ResetPlate so one change changes everything?
 * 
 * 
 */


public class Damage2D {

	// OFC Parameters
	
	private int imax, alive[], Nrichter;
	private double L, N, alpha0, alphaW, Sf0, Sr0, Sf[], Sr[], SfW, SrW, stress[], R, t1,
				   t2, dt, pt, pt2, Scum[], Scum2[], Omega, Omega2, tMAX, Omega3;
	/*
	 * 
	 * Add a SrQnch[]
	 * 
	 */
	
	private String shape, bc;
	private boolean alphaN, SfN, SrN, killed[];
	
	private Random rand;
	private LatticeNeighbors neighbors;
	
	private int nbs0, nbs[];

	// Damage Parameters

	private int Nlives0, LivesLeft[];
	private double NlivesW;
	private boolean NlivesN, crack;
	private String LifeStyle, Mode;
	
	// I/O Parameters
	
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
		String animate;
		
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
		//outfileOmega = outdir + File.separator+"Omega.txt";
		PicDir   = outdir + "/Pics/";
		DirUtil.MkDir(PicDir);
		
		if(SfW == 0){
			SfN = false;
		}
		else{
			SfN = true;
		}
		
		if(SrW == 0){
			SrN = false;
		}
		else{
			SrN = true;
		}
		
		if(alphaW == 0){
			alphaN = false;
		}
		else{
			alphaN = true;
		}
		if(Mode.equals("Damage")){
			if(LifeStyle.equals("Constant")){
				NlivesN = false;
			}
			else if(LifeStyle.equals("Gaussian")){
				NlivesN = true;
			}
			else{
				NlivesN = true;
			}
		}
		
		if(animate.equals("On")){
			showGrid = true;
		}
		else{
			showGrid = false;
		}

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
		Scum2  = new double[(int) N];
		killed = new boolean[(int) N];
				
		// Set up Lattice Neighbors
		nbs = SetUpNbs();
	
		// store initial values in arrays
		InitArrays();

		// Set up N_lives Left
		LivesLeft = SetupNlives(Mode);
		
		// Set up First Failure
		if(pt <= 0){
			ForceFailure();
			pt = 0;
		}
		else{
			ForceFailureEQ();
		}
		
		// start clocks
		t1  = 0;
		t2  = 0;
		pt2 = 0;
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
	
	private void InitArrays(){
				
		for (int jj = 0 ; jj < N ; jj++){
			
			stress[jj] = Sr0+(Sf0-Sr0)*rand.nextDouble();
			Scum[jj] = 0;
			
			alive[jj+(int)N] = 1;
			if(!NlivesN){ // no noise in the number of lives
				alive[jj] = Nlives0;
			}
			else{
				alive[jj] = (int) Math.round(Nlives0 + NlivesW*rand.nextGaussian());
			}
			
			if(!SfN){	// no noise in the critical stress
				Sf[jj] = Sf0;
			}
			else{
				Sf[jj] = Sf0 + SfW*rand.nextGaussian();
			}
			
			if(!SrN){	// no noise in the residual stress
				Sr[jj] = Sr0;
			}
			else{
				// Sr[jj] = Sr0 + SrW*rand.nextGaussian();
				Sr[jj] = Sr0 + SrW*(rand.nextDouble()-0.5);
			}
			
		}
				
		return;
	}
	
	private int[] SetupNlives(String md){
		
		
		if (md.equals("Damage")){
			
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
	
	public boolean Avalanche(){
			
			
		//System.out.println(rand.nextDouble());
		
		t1++;
		t2+=dt;
		
		Nrichter = 1; // the forced failure
		
		int[] orig = {imax};
		
		int[] avlnch = DistributeStress(orig);
		ResetPlate(orig);
		CascadeSites(avlnch);
	
		while(avlnch.length > 0){
			int[] copy = avlnch;
			avlnch = DistributeStress(avlnch);
			ResetPlate(copy);
			CascadeSites(avlnch);

		}
		
		Omega = StressMetric();
		
		ForceFailure();
				
		return (!crack);
	}
	
	private void ForceFailure(){

		double stressMax;
		imax = 0;
		
		/*
		 * replace this with
		 * minimizing the difference
		 * between the stress on a 
		 * site and that site's 
		 * failure stress 
		 * 
		 */		
		
		for (int jj = 0 ; jj < N ; jj++){
			if(stress[jj] > stress[imax]){ 
				imax = jj;
				killed[jj] = false;
			}
		}

		

		
		if(alive[imax] <= 0){
			System.out.println("All Sites Failed!");
			crack=true;
			return;
		}
				
		LivesLeft[alive[imax]]--;
		LivesLeft[alive[imax]-1]++;
		
		alive[imax]--;
		alive[imax+(int)N]=0;


		stressMax = stress[imax];		
		for (int jj = 0; jj<N; jj++){
			stress[jj] += Sf[imax]-stressMax;
		}
		dt = Sf[imax] - stressMax;

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
	
	private void ResetPlate(int[] ds){
		
		int length = ds.length;
		
		
		if(SrN){	
			for (int jj = 0 ; jj < length; jj++){
				if(alive[ds[jj]] > 0){
					stress[ds[jj]] = Sr[ds[jj]];
					Sr[ds[jj]] = Sr0 + SrW*(rand.nextDouble()-0.5); 
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
		
		
		
		return;
	}
	
	private void CascadeSites(int[] ds){
		
		if(ds == null) return;
		
		int length = ds.length;
		
		for (int jj = 0 ; jj < length ; jj++){
			alive[ds[jj]+(int)N] = 0;
			LivesLeft[alive[ds[jj]]]--;
			if(alive[ds[jj]] > 0) LivesLeft[alive[ds[jj]]-1]++;
			alive[ds[jj]]--;
			killed[ds[jj]] = false;
		}
		
		Nrichter+=length;
	
	}
	
	public boolean Equilibrate(){
		
		
		/*
		 * replace this with
		 * minimizing the difference
		 * between the stress on a 
		 * site and that site's 
		 * failure stress 
		 * 
		 */
		
		pt++;
		pt2 += dt; 
		
		Nrichter = 1; // the forced failure
		
		int[] orig = {imax};
		
		int[] avlnch = DistributeStress(orig);
		ResetPlateEQ(orig);
		CascadeSitesEQ(avlnch);
	
		while(avlnch.length > 0){
			int[] copy = avlnch;
			avlnch = DistributeStress(avlnch);
			ResetPlateEQ(copy);
			CascadeSitesEQ(avlnch);

		}
		
		Omega = StressMetric();
		
		ForceFailureEQ();
				
		return (pt < 0);
	}
	
	public boolean AvalancheEQ(){
		
		
		/*
		 * replace this with
		 * minimizing the difference
		 * between the stress on a 
		 * site and that site's 
		 * failure stress 
		 * 
		 */

		t1++;
		t2+=dt;
		pt2+=dt;
				 
		Nrichter = 1; // the forced failure
		
		int[] orig = {imax};
		
		int[] avlnch = DistributeStress(orig);
		ResetPlateEQ(orig);
		CascadeSitesEQ(avlnch);
	
		while(avlnch.length > 0){
			int[] copy = avlnch;
			avlnch = DistributeStress(avlnch);
			ResetPlateEQ(copy);
			CascadeSitesEQ(avlnch);

		}
		
		Omega  = StressMetric();
		Omega2 = StressMetricV2();
		
		ForceFailureEQ();
				
		return (t1 < tMAX);
	}
	
	private void ForceFailureEQ(){
		double stressMax;
		imax = 0;
		
		
		for (int jj = 0 ; jj < N ; jj++){
			if(stress[jj] > stress[imax]){ 
				imax = jj;
				killed[jj] = false;
			}
		}

		alive[imax+(int)N]=0;

		stressMax = stress[imax];		
		for (int jj = 0; jj<N; jj++){
			stress[jj] += Sf[imax]-stressMax;
		}
		
		dt = Sf[imax] - stressMax;

		return;
	}
	
	
	
	private void ResetPlateEQ(int[] ds){
		int length = ds.length;
		
		
		if(SrN){	
			for (int jj = 0 ; jj < length; jj++){
				stress[ds[jj]] = Sr[ds[jj]];
				Sr[ds[jj]] = Sr0 + SrW*(rand.nextDouble()-0.5); 
				alive[ds[jj]+(int)N] = 1;
			}
		}
		else{
			for (int jj = 0 ; jj < length; jj++){
				stress[ds[jj]] = Sr[ds[jj]];
				alive[ds[jj]+(int)N] = 1;
			}
		}
		
		
		
		
		return;
	}
	
	private void CascadeSitesEQ(int[] ds){
		if(ds == null) return;
		
		int length = ds.length;
		
		for (int jj = 0 ; jj < length ; jj++){
			alive[ds[jj]+(int)N] = 0;
			killed[ds[jj]] = false;
		}
		
		Nrichter+=length;
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

			
			Omega3 = Omega/(t2*t2*((int)N - LivesLeft[0]));
			
			Omega = Omega/(pt2*pt2*((int)N - LivesLeft[0]));
			
			//Sbar = Sbar/t2;
			
			//PrintUtil.printlnToFile(outfileOmega,t2,stress[(1+getL())*getL()/2],Scum[(1+getL())*getL()/2]/t2,Sbar,Omega);
			
			return Omega;
			
		}
		
		return -1;
	}
	
	private double StressMetricV2(){
		
		if(N > LivesLeft[0]){
			double Sbar = 0;

			for (int jj = 0 ; jj < N ; jj++){
				Scum2[jj] += alive[jj+(int)N]*stress[jj]*dt;
			}

			for (int jj = 0 ; jj < N ; jj++){
				Sbar += Scum2[jj];
			}

			Sbar = Sbar/((int)N - LivesLeft[0]); // all the sites save the dead ones

			Omega2 = 0;
			for (int jj = 0 ; jj < N ; jj++){
				Omega2 += (Scum2[jj] - Sbar)*(Scum2[jj] - Sbar);
			}

			Omega2 = Omega2/(t2*t2*((int)N - LivesLeft[0]));
			
			//Sbar = Sbar/t2;
			
			//PrintUtil.printlnToFile(outfileOmega,t2,stress[(1+getL())*getL()/2],Scum[(1+getL())*getL()/2]/t2,Sbar,Omega);
			
			return Omega2;
			
		}
		
		return -1;
	}
	
	public void WriteDataHeaders(){
		
		PrintUtil.printlnToFile(outfile1,"Sweeps","t_vel","pt","N_sts/avl","Omega","Omega2","Omega3");
		PrintUtil.printlnToFile(outfile2,"Sweeps","t_vel","N_lvs=0","N_lvs=1",". . . ","N_lvs=N_max");
		//PrintUtil.printlnToFile(outfileOmega,"t_vel","sigma_cntr","sigma(t)_cntr","sigma(t)--bar","Omega");
		
	}
	
	public void TakeDate(){
		PrintUtil.printlnToFile(outfile1,t1,t2,pt2,Nrichter,Omega,Omega2,Omega3);
		//PrintUtil.print2TimeAndVectorToFile(outfile2,t1,t2,LivesLeft);
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
	
}
