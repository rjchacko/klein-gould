package chris.ofc;

import java.io.File;
import java.text.DecimalFormat;
import java.util.Random;
import scikit.graphics.dim2.Grid;
import scikit.jobs.params.Parameters;
import chris.util.DirUtil;
import chris.util.LatticeNeighbors;


public class Damage2D {

	// OFC Parameters
	
	private int imax, alive[], Nrichter;
	private double L, N, alpha0, alphaW, Sc0, Sr0, Sc[], Sr[], ScW, SrW, stress[], R, t1, t2, dt, pt, Scum[];
	private String shape, bc;
	private Boolean alphaN, ScN, SrN;
	
	private Random rand;
	private LatticeNeighbors neighbors;
	
	private int nbs0, nbs[];

	// Damage Parameters

	private int Nlives0, LivesLeft[];
	private double NlivesW;
	private Boolean NlivesN, crack;
	private String LifeStyle;
	
	// I/O Parameters
	
	private String outdir, outfile1, outfile2, PicDir;
	private Boolean showGrid;
	private DecimalFormat fmtT = new DecimalFormat("000000");
	private DecimalFormat fmtH = new DecimalFormat("000");

	public Damage2D(Parameters params){
		Constructor(params);
	}

	public void Constructor(Parameters params){
		int seed;
		String animate;
		
		outdir    = params.sget("Data Directory");
		seed      = params.iget("Random Seed",0);
		shape     = params.sget("Interaction Shape");
		R         = params.fget("Interaction Radius (R)");
		L         = params.fget("Lattice Size");
		bc        = params.sget("Boundary Condtions");
		pt        = params.iget("Equilibrate Time");
		Nlives0   = params.iget("Number of Lives");
		LifeStyle = params.sget("Life Distribution");
		NlivesW   = params.fget("LD Width");
		Sc0       = params.fget("Critical Stress (\u03C3_c)");
		ScW       = params.fget("\u03C3_c width");
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
		
		if(ScW == 0){
			ScN = false;
		}
		else{
			ScN = true;
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
		
		if(LifeStyle.equals("Constant")){
			NlivesN = false;
		}
		else if(LifeStyle.equals("Gaussian")){
			NlivesN = true;
		}
		else{
			NlivesN = true;
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
		Sc     = new double[(int) N];
		Sr     = new double[(int) N];
		stress = new double[(int) N];
		Scum   = new double[(int) N];
				
		// Set up Lattice Neighbors
		nbs = SetUpNbs();
	
		// store initial values in arrays
		InitArrays();

		// Set up N_lives Left
		LivesLeft = SetupNlives();
		
		// Set up First Failure
		ForceFailure();
		
		// start clocks
		t1 = 0;
		t2 = 0;
		pt = (-1)*pt;
		
	}
	
	private int[] SetUpNbs(){
		
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
			
			stress[jj] = Sc0*rand.nextDouble();
			Scum[jj] = 0;
			
			alive[jj+(int)N] = 1;
			if(!NlivesN){ // no noise in the number of lives
				alive[jj] = Nlives0;
			}
			else{
				alive[jj] = (int) Math.round(Nlives0 + NlivesW*rand.nextGaussian());
			}
			
			if(!ScN){	// no noise in the critical stress
				Sc[jj] = Sc0;
			}
			else{
				Sc[jj] = Sc0 + ScW*rand.nextGaussian();
			}
			
			if(!SrN){	// no noise in the residual stress
				Sr[jj] = Sr0;
			}
			else{
				Sr[jj] = Sr0 + SrW*rand.nextGaussian();
			}
			
		}
		return;
	}
	
	private int[] SetupNlives(){
		
		int max = alive[0];
		int ret[];
		
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
	
	public boolean Equilibrate(){
		
		pt++;
		Nrichter = 1; // the forced failure
		
		int[] orig = {imax};
		
		int[] avlnch = DistributeStress(orig);
		ResetPlate(orig);
		CascadeSites(avlnch);
	
		while(avlnch != null){
			int[] copy = avlnch;
			avlnch = DistributeStress(avlnch);
			ResetPlate(copy);
			CascadeSites(avlnch);
		}
		
		ForceFailure();
		
		return (pt < 0);
	}
	
	public boolean Avalanche(){
		
		t1++;
		t2+=dt;
		
		Nrichter = 1; // the forced failure
		
		int[] orig = {imax};
		
		int[] avlnch = DistributeStress(orig);
		ResetPlate(orig);
		CascadeSites(avlnch);
	
		while(avlnch != null){
			int[] copy = avlnch;
			avlnch = DistributeStress(avlnch);
			ResetPlate(copy);
			CascadeSites(avlnch);
		}
		
		ForceFailure();
		
		StressMetric(); //DELETE THIS
		
		return (!crack);
	}
	
	private void ForceFailure(){

		double stressMax;
		imax = 0;
		
		for (int jj = 0 ; jj < N ; jj++){
			if(stress[jj] > stress[imax]) imax = jj;
		}

		if(alive[imax] == 0){
			System.out.println("All Sites Failed!");
			crack=true;
			return;
		}
				
		LivesLeft[alive[imax]]--;
		LivesLeft[alive[imax]-1]++;
		
		alive[imax]--;
		alive[imax+(int)N]=0;
	
		stressMax = stress[imax];
		for (int i = 0; i<N; i++){
			stress[i]+=Sc[imax]-stressMax;
			Scum[i]=stress[i];
		}
		dt = Sc[imax] - stressMax;
		
		return;
	}
	
	private int[] DistributeStress(int[] ds){
		
		int length = ds.length;
		int toKillindex = 0;
		int[] toKill = new int[(int) N];
		int[] dsnbs;
		int Nalive = 0;
		double release;
		
		for (int jj = 0 ; jj < (length - 1) ; jj++){ // loop over all dead sites
			
			if (shape.equals("All Sites")){
				Nalive = (int)N - length - LivesLeft[0];
				
				if (Nalive > 0){

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
			}
			else{	// loop over a subset of the lattice sites
				
				dsnbs = neighbors.get(ds[jj],nbs0, nbs);
				
				for (int kk = 0 ; kk < dsnbs.length ; kk++){
					Nalive += alive[dsnbs[kk]+(int)N];
				}
				
				if(alive[ds[jj]] > 0){
					release = (stress[ds[jj]] - Sr[ds[jj]])/Nalive;
				}
				else{
					release = stress[ds[jj]]/Nalive;
				}

				if(alphaN){	// alpha is a random variable
					for (int kk = 0 ; kk < dsnbs.length ; kk++){
						stress[dsnbs[kk]] += alive[dsnbs[kk]+(int)N]*release*(1 - alphaW*rand.nextGaussian() - alpha0);
					}
				}
				else{	// alpha is fixed
					for (int kk = 0 ; kk < dsnbs.length ; kk++){
						stress[dsnbs[kk]] += alive[dsnbs[kk]+(int)N]*release*(1 - alpha0);
					}
				}
			}
		}
		
		int jj = length - 1;

		if (shape.equals("All Sites")){
			Nalive = (int)N - length - LivesLeft[0];
			
			if (Nalive > 0){

				if(alive[ds[jj]] > 0){
					release = (stress[ds[jj]] - Sr[ds[jj]])/Nalive;
				}
				else{
					release = stress[ds[jj]]/Nalive;
				}

				if(alphaN){	// alpha is a random variable
					for (int kk = 0 ; kk < ds[jj] ; kk++){
						stress[kk] += alive[kk+(int)N]*release*(1 - alphaW*rand.nextGaussian() - alpha0);
						if (stress[kk] > Sc[kk]) toKill[toKillindex++] = kk;
					}
					for (int kk = ds[jj] + 1 ; kk < N ; kk++){
						stress[kk] += alive[kk+(int)N]*release*(1 - alphaW*rand.nextGaussian() - alpha0);
						if (stress[kk] > Sc[kk]) toKill[toKillindex++] = kk;
					}
				}
				else{	// alpha is fixed
					for (int kk = 0 ; kk < ds[jj] ; kk++){
						stress[kk] += alive[kk+(int)N]*release*(1 - alpha0);
						if (stress[kk] > Sc[kk]) toKill[toKillindex++] = kk;
					}
					for (int kk = ds[jj] + 1 ; kk < N ; kk++){
						stress[kk] += alive[kk+(int)N]*release*(1- alpha0);
						if (stress[kk] > Sc[kk]) toKill[toKillindex++] = kk;
					}
				}
			}
		}
		else{	// loop over a subset of the lattice sites
			
			dsnbs = neighbors.get(ds[jj],nbs0, nbs);
			
			for (int kk = 0 ; kk < dsnbs.length ; kk++){
				Nalive += alive[dsnbs[kk]+(int)N];
			}
			
			if(alive[ds[jj]] > 0){
				release = (stress[ds[jj]] - Sr[ds[jj]])/Nalive;
			}
			else{
				release = stress[ds[jj]]/Nalive;
			}

			if(alphaN){	// alpha is a random variable
				for (int kk = 0 ; kk < dsnbs.length ; kk++){
					stress[dsnbs[kk]] += alive[dsnbs[kk]+(int)N]*release*(1 - alphaW*rand.nextGaussian() - alpha0);
					if (stress[dsnbs[kk]] > Sc[dsnbs[kk]]) toKill[toKillindex++] = dsnbs[kk];
				}
			}
			else{	// alpha is fixed
				for (int kk = 0 ; kk < dsnbs.length ; kk++){
					stress[dsnbs[kk]] += alive[dsnbs[kk]+(int)N]*release*(1 - alpha0);
					if (stress[dsnbs[kk]] > Sc[dsnbs[kk]]) toKill[toKillindex++] = dsnbs[kk];
				}
			}
		}
		
		return CopyArray(toKillindex, toKill);
	}
	
	private void ResetPlate(int[] ds){
		
		int length = ds.length;
		
		for (int jj = 0 ; jj < length; jj++){
			if(alive[ds[jj]] > 0){
				stress[ds[jj]] = Sr[ds[jj]];
				alive[ds[jj]+(int)N] = 1;
			}
			else{
				stress[ds[jj]] = -2;
			}
		}
		return;
	}
	
	private void CascadeSites(int[] ds){
		
		int length = ds.length;
		
		for (int jj = 0 ; jj < length ; jj++){
			alive[ds[jj]+(int)N] = 0;
			LivesLeft[alive[ds[jj]]]--;
			if(alive[ds[jj]] > 0) LivesLeft[alive[ds[jj]]-1]++;
			alive[ds[jj]]--;
		}
		
		Nrichter+=length;
	
	}
	
	private double StressMetric(){
		
		return 0;
	}
	
	private int[] CopyArray(int length, int[] array){
		
		if (length > 0){
			int[] ret = new int[length];

			for (int jj = 0 ; jj < length ; jj++){
				ret[jj] = array[jj];
			}

			return ret;
		}
		else{
			return null;
		}
	}
	
	public void WriteDataHeaders(){
		if (showGrid){
			System.out.println(outfile1);
			System.out.println(outfile2);
			System.out.println(fmtT);
			System.out.println(fmtH);
		}
	}
	
	public void TakeDate(){
		
	}
	
	public static void PrintParams(Parameters params){
		
	}
	
	public void TakePicture(Grid grid){
		
	}
	
}
