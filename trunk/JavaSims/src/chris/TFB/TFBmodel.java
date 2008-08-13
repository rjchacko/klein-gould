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

public class TFBmodel {

	private double 	stress, L, N, Stot, beta, kappa, Sf[], Sf0, SfW, t, sCum[], 
					Omega, energy, D, probF, probH;
	private int 	state[], Nalive;
	private String 	dir, of;
	private boolean draw;
	private Random 	rand;
	
	private Accumulator gplotdata = new Accumulator(1/N);
	
	private DecimalFormat fmt = new DecimalFormat("0000000");
	
	public TFBmodel(Parameters params){
		
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
		D      = params.fget("D");
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
	
	private void initArrays(){
		
		for (int jj = 0 ; jj < N ; jj++){
			Sf[jj]    = Sf0 + SfW*(0.5 - rand.nextDouble());
			state[jj] = 1;
			sCum[jj]  = 0;
		}
		
		plotGlandscape();
		
		return;
	}
	
	public boolean nextLattice(){
		
		t++;
		int deltaN = 0;
		
		probF = pFail(Nalive/N);
		probH = pHeal(Nalive/N);
		
//		for(int jj = 0 ; jj < N ; jj++){
//			if(nextState(jj)){
//				state[jj] = 1-state[jj];	// state = 0 --> 1 // state = 1 --> 0 
//				deltaN += 2*state[jj] - 1;  // state = 0 --> -1 // state = 1 --> +1
//			}
//		}
		
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
			
		if(state[st] == 1){
			// If site passes failure requirement, return TRUE   
			if(rand.nextDouble() < probF) return true;
			return false;
		}
		else{
			// If site passes healing requirement, return TRUE   
			if(rand.nextDouble() < probH) return true;
			return false;
		}
	}
	
	private double pFail(double phiold){
		
		double phinew = (N*phiold-1)/N;
		
		if(phinew==0) return Double.MAX_VALUE;
		
		if(phiold==1) return Math.exp(-beta*(-0.5*(Stot*Stot/kappa)*(1/phinew - 1/phiold) - D*(phinew - phiold) +
		(1/beta)*(phinew*Math.log(phinew)  + (1 - phinew)*Math.log(1 - phinew) )));
		
		double dg = -0.5*(Stot*Stot/kappa)*(1/phinew - 1/phiold) - D*(phinew - phiold) +
					(1/beta)*(phinew*Math.log(phinew) - phiold*Math.log(phiold) + 
					(1 - phinew)*Math.log(1 - phinew) - (1 - phiold)*Math.log(1 - phiold));
		
		return Math.exp(-beta*dg);
	}
	
	private double pHeal(double phiold){
		
		double phinew = (N*phiold+1)/N;
		
		if (phiold==0) return 0;
		
		if(phinew==1) return Math.exp(-beta*(-0.5*(Stot*Stot/kappa)*(1/phinew - 1/phiold) - D*(phinew - phiold) +
							(-1/beta)*(phiold*Math.log(phiold) + (1 - phiold)*Math.log(1 - phiold))));
		
		double dg = -0.5*(Stot*Stot/kappa)*(1/phinew - 1/phiold) - D*(phinew - phiold) +
					(1/beta)*(phinew*Math.log(phinew) - phiold*Math.log(phiold) + 
					(1 - phinew)*Math.log(1 - phinew) - (1 - phiold)*Math.log(1 - phiold));
		
		return Math.exp(-beta*dg);
	}
	
	public double gEnergy(double phi){
		
		if(phi==0) return (-1)*Double.MAX_VALUE;
		if(phi==1) return -0.5*(Stot*Stot/kappa);
		
		return -0.5*(Stot*Stot/kappa)*(1/phi) - D*phi +
				(1/beta)*(phi*Math.log(phi) + (1 - phi)*Math.log(1 - phi));
	}
	
	public double hamiltonian(double phi){
		
		return (0.5*Stot*Stot/(kappa*phi) - D*phi);
	}
	
	public static void printParams(String fout, Parameters prms){
		
		PrintUtil.printlnToFile(fout, prms.toString());
		return;
	}

	public void plotGlandscape(){
		
		double tmp[] = new double[(int)N + 1];
		
		for(int jj = 0 ; jj < N+1 ; jj++){
			gplotdata.accum(jj,gEnergy(jj/N));
			tmp[jj] = gEnergy(jj/N);
		}
		
		PrintUtil.printArrayToFile(dir+File.separator+"FEdata.txt",tmp,(int)N+1,1);
		
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
	
}
