package russ.ofcdamage.apps;

import java.io.File;
import java.text.DecimalFormat;

import russ.ofcdamage2.damage2Dfast;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.Parameters;
import chris.util.LatticeNeighbors;
import chris.util.MathUtil;
import chris.util.PrintUtil;
import chris.util.dummyParamUtil;

public class fastCorrDamageApp extends Simulation{
	
	private int eqt, Ndead, L;
	private double N;
	private double grd[][] = new double[10][201];
	private boolean mes[] = new boolean[11];
	private damage2Dfast model;
	private LatticeNeighbors nbs;
	private DecimalFormat ifmt = new DecimalFormat("###");
	
	
	public static void main(String[] args) {
		new Control(new fastCorrDamageApp(), "Damage Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
//		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)");
		params.add("Lattice Size");
		params.set("Interaction Radius (R)", (int) 20);
		params.set("Lattice Size", (int) 512);	
//		params.add("Boundary Condtions", new ChoiceValue("Periodic","Open"));
		params.add("Equil Time", 600000);
//		params.add("Sim Time", 500000);
//		params.add("Number of Lives",(int) 10);
//		params.add("NL width", (int) 8);
		params.add("Failure Stress (\u03C3_f)", 2.);
//		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.2);
		params.add("Dissipation (\u03B1)", 0.05);
//		params.add("\u03B1 width", 0.);
		params.add("Cycle");
		params.add("Mode");
		
		return;
	}
	
	public void run() {

		int ccl  = 0;
		int ndo;
		
		while(ccl < 1000){	
			
			for (int jj = 0 ; jj < 10 ; jj++){
				mes[jj] = true;
			}
			mes[10] = false;

			ccl++;	
			params.set("Random Seed", params.iget("Random Seed")+1);

			// Setup model
			params.set("Mode", "Intializing");
			Job.animate();
			
			L      = params.iget("Lattice Size");
			N      = L*L;
			Parameters dparams = dummyParamUtil.ofcParams(params);
			model  = new damage2Dfast(dparams);
			if(ccl == 1) model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".log",params,model);	
			eqt    = params.iget("Equil Time");
			Ndead  = 0;
			double[] dummy = new double[(int)(0.45*N)];
			int dummycount = 0;

			
			params.set("Cycle", ccl);
			params.set("Mode","Equilibrate");
			Job.animate();

			// Equilibrate the system

			for (int jj = 0 ; jj < eqt ; jj++){
				model.evolve(jj,false);
				if(jj%1000 == 0){
					params.set("Mode", (jj-eqt));
					Job.animate();
				}
			}

			// Simulate the model with damage
			
			ndo = 0;
			params.set("Mode","Damage");
			Job.animate();
			while(Ndead < N){
				Ndead = model.evolveD(1,false);
				if(Ndead - ndo > 1000){
					ndo = Ndead;
					params.set("Mode", Ndead);
					Job.animate();
				}
				// save lattice (ugh)
				// calculate g(r) and a phi
				if((double)(Ndead)/N > 0.475) continue;
				int phiIndex = getPhiIndex();
				if(mes[phiIndex]){ // check to make sure you aeren't 'double' counting / phi is in range
					mes[phiIndex] = false;
					calcGr(phiIndex);
					params.set("Mode", Ndead);
					Job.animate();
				}
				// save a list of phi for bug checking
				dummy[dummycount++] = Ndead/N;
			}


			params.set("Mode", "Next Cycle");
			params.set("Random Seed", params.iget("Random Seed") + 1);
			// print the data so far
			PrintUtil.printMxNarrayToFile(model.getOutdir()+File.separator+model.getBname()+"_"+ifmt.format(ccl)+".txt", grd, 201, 10);
			PrintUtil.printVectorToFile(model.getOutdir()+File.separator+model.getBname()+"phiStop"+"_"+ifmt.format(ccl)+".txt", dummy,dummycount);
			Job.animate();

		}
		params.set("Mode", "Done");
		Job.animate();
		return;
	}

	public void animate() {
		
		return;
	}

	public void clear() {

		return;
	}
	
	private int getPhiIndex(){
		
		double tmp = 1 - Ndead/N;
		if( .975 > tmp && .925 < tmp) return 0;
		else if( .925 > tmp && .875 < tmp) return 1;
		else if( .875 > tmp && .825 < tmp) return 2;
		else if( .825 > tmp && .775 < tmp) return 3;
		else if( .775 > tmp && .725 < tmp) return 4;
		else if( .725 > tmp && .675 < tmp) return 5;
		else if( .675 > tmp && .625 < tmp) return 6;
		else if( .625 > tmp && .575 < tmp) return 7;
		else if( .575 > tmp && .525 < tmp) return 8;
		else if( .525 > tmp && .475 < tmp) return 9;
		return 10;
	}
	
	public void calcGr(int idx){
		

		// first calculate g(r) from this lattice
		double gr[] = new double[200];
			// loop over the radii 1 <= r <= 200
			for (int rr = 1 ; rr < 201 ; rr++){
				params.set("Mode", "Calculating g(r)");
				Job.animate();
				// get neighbors at R = rr
				nbs = new LatticeNeighbors((int) L,(int) L,rr,rr+1,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Circle);
				// loop over every dead site on the lattice
				int sum   = 0;
				int count = 0;
				for (int jj = 0 ; jj < N ; jj++){
					if(model.isAlive(jj)) continue;
					int[] ln = nbs.get(jj);
					count = ln.length;
					for (int kk = 0 ; kk < count ; kk++){
						sum += MathUtil.bool2bin(!(model.isAlive(ln[kk])));
					}
				}
				gr[rr-1] = (double)(sum)/(Ndead*count);
			}
		// now add this measurement to previous ones	
		if(grd[idx][0] > 0){
			// do the average properly
			int M = (int)(grd[idx][0]);
			for(int jj = 1 ; jj < 201 ; jj++){
				grd[idx][jj] = (M*grd[idx][jj]+gr[jj-1])/(M+1);
			}	
		}
		else{
			// copy the data into the array
			for(int jj = 1 ; jj < 201 ; jj++){
				grd[idx][jj] = gr[jj-1];
			}
		}
		grd[idx][0]++;	// update the number of data sets
		
		return;
	}
	
	
	// end of class
}
