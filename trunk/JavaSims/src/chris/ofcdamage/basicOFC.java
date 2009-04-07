package chris.ofcdamage;

import java.util.Random;

import scikit.jobs.params.Parameters;
import chris.util.CopyUtil;
import chris.util.LatticeNeighbors;
import chris.util.PrintUtil;

public class basicOFC {	
	
	// Model Parameters
	private double  R, L, N, Sr0, Sr[], SrW, Sf0, alpha0, stress[], Scum[], t, dt, Omega;
	private boolean SrN, clocks, ks[];
	private Random  rand;
	private String  shape, bc;
	private int     alive[], nbs0, nbs[];
	private LatticeNeighbors neighbors;
	
	
	public basicOFC(Parameters params){
	
		SKdamageConstr(params);
		return;
	}
	
	public void SKdamageConstr(Parameters params){
	
		int seed;
		
		seed   = params.iget("Random Seed");
		shape  = params.sget("Interaction Shape");
		R      = params.fget("Interaction Radius (R)");
		L      = params.iget("Lattice Size");
		bc     = params.sget("Boundary Condtions");
		Sf0    = params.fget("Failure Stress (\u03C3_f)");
		Sr0    = params.fget("Residual Stress (\u03C3_r)");
		SrW    = params.fget("\u03C3_r width");
		alpha0 = params.fget("Dissipation (\u03B1)");
		
		rand = new Random(seed);

		t      = 0;
		N      = L*L;
		SrN    = (SrW > 0);
		clocks = false;
		
	
		InitializeArrays();
		
		return;
	}
	
	private void InitializeArrays(){
		
		Sr             = new double[intN()];
		stress         = new double[intN()];
		Scum           = new double[intN()];
		alive          = new int[intN()];
		ks             = new boolean[intN()];
		
		for (int jj = 0 ; jj < intN() ; jj++){
			alive[jj]  = 1;
			stress[jj] = Sr0 + (Sf0 - Sr0)*rand.nextDouble();
			if(SrN){
				Sr[jj] = Sr0 + SrW*(0.5 - rand.nextDouble());
			}
			else{
				Sr[jj] = Sr0;
			}	
			ks[jj] = false;
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
		
		int[] seedsites = forceFailure();

		while((seedsites = distributeStress(seedsites)) != null);
		
		return;
	}
	
	public boolean avalanche(){

		int[] seedsites = forceFailure();

		while((seedsites = distributeStress(seedsites)) != null);
		
		ergMetric();
		
		return true;
	}
	
	private int[] forceFailure(){

		double stressMax;
		int imax = 0;
		dt   = 0;
		
		for (int jj = 0 ; jj < N ; jj++){
			if((Sf0 - stress[jj]) < (Sf0 - stress[imax])) imax = jj;
		}
		
		stressMax = stress[imax];		
		for (int jj = 0; jj<N; jj++){
			stress[jj] += Sf0-stressMax;
		}
		
		if(clocks){	
			dt = Sf0 - stressMax;
			t += dt;
		}

		return CopyUtil.copyArray(imax,1);
	}
	

	private int[] distributeStress(int[] sites){
		
		int Nalive, dsnbs[], newlykilled[], counter;
		double release;
		
		newlykilled = new int[intN()];
		counter     = 0;
		
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
				stress[dsnbs[kk]] += alive[dsnbs[kk]]*release*(1-alpha0);
				if(!(ks[dsnbs[kk]]) &&  alive[dsnbs[kk]]*stress[dsnbs[kk]] > Sf0){
					ks[dsnbs[kk]] = true;
					newlykilled[counter++] = dsnbs[kk];
				}
			}
			
		}
		
		resetSites(sites);
		return CopyUtil.copyArray(newlykilled,counter);
	}
	
	private void resetSites(int[] sites){
		
		for (int jj = 0 ; jj < sites.length ; jj++){
			alive[sites[jj]]  = 1;
			Sr[sites[jj]]     = nextSr(sites[jj]);
			stress[sites[jj]] = Sr[sites[jj]];
			ks[sites[jj]]     = false;
		}
		
		return;
	}
	
	private void failSites(int[] sites){
		
		for (int jj = 0 ; jj < sites.length ; jj++){
			alive[sites[jj]] = 0;
		}
		
		
		
		
		return;
	}
	
	private double nextSr(int site){
			
		return (SrN) ? Sr0 + SrW*(0.5 - rand.nextDouble()) : Sr0;
	}
	
	public int intN(){
		return (int) N;
	}

	public static void printParams(String fout, Parameters prms){
		
		PrintUtil.printlnToFile(fout,prms.toString());
		return;
	}
	
	
	public int getL(){
		
		return (int) L;
	}
	
	public double[] getStress(){
		
		return stress;
	}

	
	private void ergMetric(){
		
		if(!(clocks)) return;
		
		double Sbar = 0;

		for (int jj = 0 ; jj < N ; jj++){
			Scum[jj] += alive[jj]*stress[jj]*dt;
			Sbar += alive[jj]*Scum[jj];
		}
		
		Sbar = Sbar/N;
		Omega = 0;
		
		for (int jj = 0 ; jj < N ; jj++){
			Omega += alive[jj]*(Scum[jj] - Sbar)*(Scum[jj] - Sbar);
		}

		Omega = Omega/(t*t*N);

		return;
	}

	public void runClocks(boolean bool){
		clocks = bool;
	}
	
	public double getSimTime(){

		return t;
	}
	
	public double getMetric(){

		return Omega;
	}
	
}
