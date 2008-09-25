package chris.TFB;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Random;

import javax.imageio.ImageIO;

import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.jobs.params.Parameters;
import chris.util.PrintUtil;

public class DTFBmodel {

	private double 	stress, L, N, Stot, beta, kappa, Sf[], Sf0, SfW, t, sCum[], 
					Omega, energy, D0, D[], Dw, Nalive, elastic;
	private int 	state[];
	private String 	dir, of;
	private boolean draw;
	public Random 	rand;
	
	private Accumulator gplotdata = new Accumulator(1/N);
	
	private DecimalFormat fmt = new DecimalFormat("0000000");
	
	public DTFBmodel(Parameters params){
		
		tfbConstructor(params);
		
		return;
	}
	
	public void tfbConstructor(Parameters params){
		
		stress = params.fget("Stress / site");
		L      = params.iget("Lattice Size");
		N      = L*L;
		Stot   = N*stress;
		beta   = 1/(params.fget("Temperature"));
		kappa  = params.fget("Kappa");
		D0     = params.fget("D");
		Dw     = params.fget("D width");
		rand   = new Random(params.iget("Random Seed"));
		Sf     = new double[(int) N];
		state  = new int[(int) N];
		sCum   = new double[(int) N];
		Sf0    = params.fget("Failure Stress");
		SfW    = params.fget("\u03C3_f width");
		Nalive = (int) N;
		draw   = (params.sget("Animation").equals("On"));
		dir    = params.sget("Data Directory");
		t      = 0;
		of     = dir + File.separator + "TFBdata.txt";		
		
		initArrays();
	}
	
	protected void initArrays(){
		
		energy  = 0;
		elastic = 0.5*Stot*Stot/kappa;
		for (int jj = 0 ; jj < N ; jj++){
			Sf[jj]    = Sf0 + SfW*(0.5 - rand.nextDouble());
			state[jj] = 1;
			sCum[jj]  = 0;
			D[jj]     = D0 * Dw*(rand.nextDouble() - 0.5);
			energy    += elastic - D[jj];
		}
				
		return;
	}
	
	public boolean nextLattice(){
		
		t++;
		int deltaN = 0;
		
		int testsite = rand.nextInt((int) N);
		if(nextState(testsite)){
			state[testsite] = 1-state[testsite];	// state = 0 --> 1 // state = 1 --> 0 
			deltaN += 2*state[testsite] - 1;  // state = 0 --> -1 // state = 1 --> +1
		}	
		
		Nalive += deltaN;
		ergMetric();
		
		return (Nalive > 0);
	}
	
	private boolean nextState(int st){
			
		double dE = (2*state[st] - 1)*(elastic - D[st]);

		if (dE > 0) return true;
		if (rand.nextDouble() < Math.exp(-beta*dE)) return true;
		return false;
		
	}
	
//	private double pFail(double phiold){
//		
//		double phinew = (N*phiold-1)/N;
//		
//		double dg = gEnergy(phinew) - gEnergy(phiold);
//		
//		return Math.exp(-beta*dg);
//	}
//	
//	private double pHeal(double phiold){
//		
//		double phinew = (N*phiold+1)/N;
//		
//		double dg = gEnergy(phinew) - gEnergy(phiold);
//		
//		return Math.exp(-beta*dg);
//	}
	
//	private double pFail(double phiold){
//		
//		double phinew = (N*phiold-1)/N;
//		
//		if (phinew == 0) return 0;
//		
//		double dg = hamiltonian(phinew) - hamiltonian(phiold);
//		
//		return Math.exp(-beta*dg);
//	}
//	
//	private double pHeal(double phiold){
//		
//		if (phiold == 0) return 100;
//		
//		double phinew = (N*phiold+1)/N;
//		
//		double dg = hamiltonian(phinew) - hamiltonian(phiold);
//		
//		return Math.exp(-beta*dg);
//	}

	
	public static void printParams(String fout, Parameters prms){
		
		PrintUtil.printlnToFile(fout, prms.toString());
		return;
	}
	
	public void writeHeader(){
		
		PrintUtil.printlnToFile(of, "Time", "Phi", "Energy", "Omega");
		return;
	}
	
	public void takeData(){
		
		PrintUtil.printlnToFile(of, t, Nalive/N, energy, Omega);
		return;
	}
	
	public void takePicture(Grid grid){

		if(!draw) return;
		
		String SaveAs = dir + File.separator + grid.getTitle()+fmt.format(t)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}

		return;
	}
	
	private void ergMetric(){
		
		double sbar = 0;
		
		if(Nalive > 0){
			
			for (int jj = 0 ; jj < N ; jj++){
				sCum[jj] += state[jj];
				sbar += sCum[jj];
			}
			
			sbar = sbar/Nalive;		
			Omega = 0;
			
			for (int jj = 0 ; jj < N ; jj++){
				Omega += (sCum[jj] - sbar)*(sCum[jj] - sbar);
			}
			
			Omega = Omega/(t*t*Nalive);
		}
		
		return;
	}
	
	public double getTime(){
		
		return t;
	}
	
	public double getPhi(){
		
		return Nalive/N;
	}
	
	public int[] getState(){
		
		return state;
	}
	
	public int getL(){
		
		return (int) L;
	}
	
	public Accumulator getGplot(){
		
		return gplotdata;
	}
	
	public int getN(){
		
		return (int) N;
	}
	
	public String getdir(){
	
		return dir;
	}
	
	public double getST(){
		
		return Stot;
	}
	
	public double getK(){
		
		return kappa;
	}
	
}
