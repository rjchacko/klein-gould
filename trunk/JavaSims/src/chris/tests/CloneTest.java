package chris.tests;


import java.io.File;
import java.text.DecimalFormat;
import java.text.Format;

import scikit.graphics.ColorGradient;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.ofc2Dfast;
import chris.util.CloneUtil;
import chris.util.MathUtil;
import chris.util.PrintUtil;
import chris.util.dummyParamUtil;

public class CloneTest extends Simulation {
	
	private ofc2Dfast model, cmodel;
	private int eqt, simt, ct;
	private boolean draw;
	
	Grid gridO = new Grid ("Original");
	Grid gridC = new Grid ("Clone");
	Grid gridOO = new Grid ("Original");
	Grid gridCC = new Grid ("Clone");
	
	Format df  = new DecimalFormat("0.000000000000");
	Format ift = new DecimalFormat("000");
	
	
	public static void main(String[] args) {
		new Control(new CloneTest(), "Random Clones App");
	}

	public void load(Control c) {
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/test/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int) 10);
		params.add("Lattice Size", (int) 256);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Open"));
		params.add("Equil Time", 100000);
		params.add("Sim Time", 100000);
		params.add("Failure Stress (\u03C3_f)", 2.);
		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.05);
		params.add("Dissipation (\u03B1)", 0.05);
		params.add("\u03B1 width", 0.);
		params.add("Status");
		
		c.frameTogether("Stress",gridO,gridC);
		c.frameTogether("Event",gridOO,gridCC);

	}

	public void animate() {

		if(draw){
			int L = model.getL();
			
			gridO.registerData(L,L,model.getStress());
			gridC.registerData(L, L, cmodel.getStress());
			gridOO.registerData(L,L,MathUtil.bool2bin(model.getLastShower()));
			gridCC.registerData(L, L, MathUtil.bool2bin(cmodel.getLastShower()));
		}
		
		return;
	}

	public void clear() {
		
		gridO.clear();
		gridC.clear();
		gridOO.clear();
		gridCC.clear();
		return;
	}

	public void run() {
		
		// Setup model
		draw = false;
		params.set("Status", "Intializing");
		Job.animate();
		model  = new ofc2Dfast(dummyParamUtil.ofcParams(params));
		cmodel = null;
		model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".txt",params);	
		eqt   = params.iget("Equil Time");
		simt  = params.iget("Sim Time");
		ct    = 5000;
		params.set("Status", "Ready");
		Job.animate();
		
		// Equilibrate the system
		for (int jj = 0 ; jj < eqt ; jj++){
			model.evolve(jj,false);
			if(jj%500 == 0){
				params.set("Status", (jj-eqt));
				Job.animate();
			}
		}
		
		// Simulate the model without damage
		for (int jj = 0 ; jj < simt ; jj++){
			model.evolve(jj,false);
			if(jj%500 == 0){
				params.set("Status", jj);
				Job.animate();
			}
			if(jj % ct == 0 && jj > 0){
				params.set("Status", "Testing Clone");
				Job.animate();

				cmodel = CloneUtil.cloneOFC(params, model, true);
				// check 1 (make sure random numbers match up)
				double r, rc, g, gc;
				int i, ic;
				r  = model.getRand().nextDouble();
				g  = model.getRand().nextGaussian();
				i  = model.getRand().nextInt(100);
				rc = cmodel.getRand().nextDouble();
				gc = cmodel.getRand().nextGaussian();
				ic = cmodel.getRand().nextInt(100);
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Doubles");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(r));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(rc));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(r-rc));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Gaussians");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(g));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(gc));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(g-gc));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Integers");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",ift.format(i));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",ift.format(ic));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",ift.format(i-ic));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");		
				// check 2 (check stress on each site and make sure they are the same)
				double maxs = 0;
				double maxr = 0;			
				for (int kk = 0 ; kk < model.getN() ; kk++){
					if(Math.abs(model.getStress(kk) - cmodel.getStress(kk)) > maxs) maxs = model.getStress(kk) - cmodel.getStress(kk);
					if(Math.abs(model.getSr(kk) - cmodel.getSr(kk)) > maxr) maxr = model.getSr(kk) - cmodel.getSr(kk);
				}
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Max(s_i - sc_i) = ");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(maxs));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Max(sR_i - scR_i) = ");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(maxr));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","s_7 =", df.format(model.getStress(7)) ,df.format(cmodel.getStress(7)));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","sR_7 =", df.format(model.getSr(7)) ,df.format(cmodel.getSr(7)));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","GR: ");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",model.getGR(),cmodel.getGR());
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");
				// check 3 (check some model parameters)
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","N:",model.getN(),cmodel.getN());
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Sr0: ",model.getSr0(),cmodel.getSr0());
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Sf0: ",model.getSf0(),cmodel.getSf0());
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","a: ",model.getAlpha0(),cmodel.getAlpha0());
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Evolving");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");
				// evolve each model for 5000 plate updates
				draw = false;
				for (int kk = 0 ; kk < 5000 ; kk++){ // YOU HERE
					model.evolve(kk,true);
					cmodel.evolve(kk,true);
					if(kk % 500 == 0){
						params.set("Status", "Testing Clone" + ift.format(100*kk/5000));
						Job.animate();
					}
				}
				draw = false;
				// re-check systems
				// check 1 (make sure random numbers still match up)
				r  = model.getRand().nextDouble();
				g  = model.getRand().nextGaussian();
				i  = model.getRand().nextInt(100);
				rc = cmodel.getRand().nextDouble();
				gc = cmodel.getRand().nextGaussian();
				ic = cmodel.getRand().nextInt(100);
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Doubles");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(r));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(rc));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(r-rc));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Gaussians");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(g));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(gc));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(g-gc));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Integers");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",ift.format(i));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",ift.format(ic));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",ift.format(i-ic));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");		
				// check 2 (check stress on each site and make sure they are the same)
				maxs = 0;
				maxr = 0;			
				for (int kk = 0 ; kk < model.getN() ; kk++){
					if(Math.abs(model.getStress(kk) - cmodel.getStress(kk)) > maxs) maxs = model.getStress(kk) - cmodel.getStress(kk);
					if(Math.abs(model.getSr(kk) - cmodel.getSr(kk)) > maxr) maxr = model.getSr(kk) - cmodel.getSr(kk);
				}
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Max(s_i - sc_i) = ");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(maxs));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Max(sR_i - scR_i) = ");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(maxr));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","s_7 =", df.format(model.getStress(7)) ,df.format(cmodel.getStress(7)));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","sR_7 =", df.format(model.getSr(7)) ,df.format(cmodel.getSr(7)));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",df.format(model.getSr(7) - cmodel.getSr(7)));
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","GR: ");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt",model.getGR(),cmodel.getGR());
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");
				// check 3 (check some model parameters)
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","N:",model.getN(),cmodel.getN());
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Sr0: ",model.getSr0(),cmodel.getSr0());
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Sf0: ",model.getSf0(),cmodel.getSf0());
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","a: ",model.getAlpha0(),cmodel.getAlpha0());
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","*********************************************");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","Done");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","*********************************************");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");
				PrintUtil.printlnToFile("/Users/cserino/Desktop/check.txt","---------------------------------------------");
				params.set("Status", jj);
				Job.animate();
			}
		}
		params.set("Status", "Done");
		Job.signalStop();
		Job.animate();
	}
	
	public void setupIO(){
		
		ColorGradient cg = new ColorGradient();
		double sr0, dsr, sf0, dsf, clr;
		int N = model.getN();

		gridO.setAutoScale(false);
		gridC.setAutoScale(false);
		gridOO.setAutoScale(false);
		gridCC.setAutoScale(false);
		
		sf0 = params.fget("Failure Stress (\u03C3_f)");
		dsf = params.fget("\u03C3_f width");
		sr0 = params.fget("Residual Stress (\u03C3_r)");
		dsr = params.fget("\u03C3_r width");

		gridOO.setScale(0,1);
		gridCC.setScale(0,1);
		gridO.setScale(sr0-dsr,sf0+dsf);
		gridC.setScale(sr0-dsr,sf0+dsf);
		for (int jj = 0 ; jj <= N ; jj++){
			clr = sr0-dsr + jj*(sf0+dsf-sr0+dsr)/N;
			cg.getColor(clr,sr0-dsr,sf0+dsf);
			gridO.setColors(cg);
			gridC.setColors(cg);
		}
		
		return;
	}	

}
