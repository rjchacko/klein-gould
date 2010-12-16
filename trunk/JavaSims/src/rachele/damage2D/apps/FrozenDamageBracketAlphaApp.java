package rachele.damage2D.apps;

import java.io.File;

import rachele.damage2D.multidx.FrozenDamageLattice;
import rachele.util.FileUtil;
//import scikit.dataset.Accumulator;
//import scikit.dataset.DatasetBuffer;
//import scikit.dataset.Histogram;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.util.DoubleArray;
//import scikit.util.Utilities;

public class FrozenDamageBracketAlphaApp extends Simulation{
	int dt;
	int pow;
	int maxSize;	
//	int cg_dt;

	String infoFile;
	
	double stressMetric;
	double alphaDissRate;
		
	FrozenDamageLattice ofc;
	Grid deadGrid = new Grid("Lattice");

	public static void main(String[] args) {
		new Control(new FrozenDamageBracketAlphaApp(), "OFC Damage Model Bracketing App");
	}

	public void load(Control c) {
		deadGrid = new Grid("Dead Blocks dx = " + 1);
		c.frame( deadGrid);

		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/contract2/testruns/"));
		String rd = "Random";
		String br = "Random Blocks";
		String ds = "Dead Strip";
		String pr = "Place Random Dead";
		String db = "Dead Block";
		String cs = "Cascade";
		String dr = "Dead Rectangle";
		String cr = "Cascade Random";
		String bl = "Dead Blocks";
		String pd = "Place Dead Blocks";
		params.add("Type of Damage", new ChoiceValue( cr, pd, bl, cs, br, ds, pr, rd, db, dr));
		params.add("Dead dissipation?", new ChoiceValue("Yes", "No") );
		params.add("Boundary Conditions", new ChoiceValue("Periodic", "Open"));
		params.addm("Random Seed", 1);
		params.addm("Size Power",8);
		params.addm("R", 16);
		params.addm("Init Percent Dead", 0.25);
		params.addm("Dead Parameter", 4);
		params.addm("Number Dead", 1024);
//		params.addm("Coarse Grained dt (PU)", 1);
//		params.addm("Equilibration Updates", 100000);
//		params.addm("Max PU",1000000);
//		params.addm("Data points per write", 100000);
		params.addm("Residual Stress", 0.625);
		params.addm("Res. Max Noise", 0.125);
		params.addm("Dissipation Param", 0.0);
		params.addm("Target Alpha", 0.2);
		params.addm("Tolerance", 0.01);
//		params.addm("Damage Tolerance", 1.0);
//		params.add("Av Size");
//		params.add("Plate Updates");
		params.add("Percent Damage");
		params.add("Alpha Prime");
		params.add("Iterations");
//		params.add("Error");
	}

	public void animate() {
		params.set("Alpha Prime", ofc.alphaPrime);

	}

	public void clear() {
	}

	
	public void run() {
		infoFile = params.sget("Data Dir") + File.separator + "info.txt";
//		String alphaHistFile = params.sget("Data Dir") + File.separator + "ah.txt";
//		ofc = new FrozenDamageLattice(params, infoFile);
//		FileUtil.initFile(alphaHistFile, params);
//		FileUtil.printHistToFile(alphaHistFile, ofc.alpha_iHist);
//		ofc.alpha_iHist.clear();		
//		pow = params.iget("Size Power");
		
		double target = params.fget("Target Alpha");
		double tol = params.fget("Tolerance");
		double upBound = target + tol;
		double lowBound = target - tol;
	
		boolean outOfBounds = true;
		int it = 0;
		boolean over;
		double ap=0;
		ofc = new FrozenDamageLattice(params, infoFile);
		ap=ofc.alphaPrime;
		double lastAp = ap;
		drawLattices();
		it += 1;
		double ipd = params.fget("Init Percent Dead");
		if(ap > lowBound && ap < upBound) outOfBounds = false;
		//if out of bounds, find initial high and low ipd.
//		double hBracket = 0.99999;
//		double lBracket = 0.00001;
		double step = 0.01;
		if (ap > target) over = true;
		else over =false;
		while(outOfBounds){

			if(over){
				ipd -= step;
				ofc = new FrozenDamageLattice(params, infoFile);
				ap=ofc.alphaPrime;
				if (ap < target){
					step /=2.0;
					over=false;
				}
			}else{
				ipd += step;
				ofc = new FrozenDamageLattice(params, infoFile);
				ap=ofc.alphaPrime;
				if(ap > target){
					step /=2.0;
					over = true;
				}
			}

			drawLattices();
			it += 1;
			params.set("Init Percent Dead", ipd);
			params.set("Iterations", it);
			if(ap > lowBound && ap < upBound) outOfBounds = false;	
			System.out.println("over = " + over);
			
			Job.animate();
			lastAp = ap;
		}
		System.out.println("DONE after " + it + " iteration.");		

		if(ap > lowBound && ap < upBound){
			System.out.println("DONE");
		}else if(ap > upBound){
			over = true;
		}else if(ap < lowBound){
			over = false;
		}
		
		while(true){
			
			Job.animate();
		}
	}
	
	void drawLattices(){

			int dx = 1;
			int Lp = ofc.L/dx;
			int Np = Lp*Lp;
			double [] deadSites = new double [Np];
			for (int i = 0; i < Np; i++){
				if(ofc.aliveLattice[i]) deadSites[i]=1;
				else deadSites[i]=0;
			}
			double percentAlive = DoubleArray.sum(deadSites)/(double)Np;
			deadGrid.setScale(0.0, 1.0);
			deadGrid.registerData(Lp, Lp, deadSites);
			FileUtil.printlnToFile(infoFile, "# Percent alive for dx=" + dx + " is " + percentAlive);
	}


}
