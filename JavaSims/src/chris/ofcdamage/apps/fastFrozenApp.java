package chris.ofcdamage.apps;

import java.io.File;
import java.text.DecimalFormat;

import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.Parameters;
import chris.ofcdamage.damage2Dfast;
import chris.util.PrintUtil;
import chris.util.dummyParamUtil;

public class fastFrozenApp extends Simulation{

	private int simt, eqt;
	private DecimalFormat pfmt = new DecimalFormat("000");
	private damage2Dfast model;
	
	// temp change
	private Histogram hes = new Histogram(1.);
	

	public static void main(String[] args) {
		new Control(new fastFrozenApp(), "Damage Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
//		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int) 10);
		params.add("Lattice Size", (int) 256);
		params.add("Pass Stress to", new ChoiceValue("Live Sites","All Sites"));
//		params.add("Boundary Condtions", new ChoiceValue("Periodic","Open"));
		params.add("Equil Time", 500000);
		params.add("Sim Time", 500000);
//		params.add("Number of Lives",(int) 1);
//		params.add("NL width", (int) 0);
		params.add("Failure Stress (\u03C3_f)", 2.);
//		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.1);
		params.add("Dissipation (\u03B1)", 0.05);
//		params.add("\u03B1 width", 0.);
		params.add("Mode");
		params.add("\u03D5");
		params.add("Status");
		
		return;
	}
	
	public void run() {
		
		double phin = 1;
		double dphi = 0.1; 
		simt = params.iget("Sim Time");
		eqt  = params.iget("Equil Time");
		
		
		if(params.sget("Pass Stress to").equals("Live Sites")){
		
			while (phin > 0){

				params.set("Mode", "Freeze");
				params.set("\u03D5",1-phin);
				Parameters dparams = dummyParamUtil.ofcParams(params);
				model  = new damage2Dfast(dparams);
				if(phin == 1) model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".log",params,model);	

				// equilibrate model			
				params.set("Mode", "Equilibrating");
				params.set("\u03D5",phin);
				Job.animate();
				for (int jj = 0 ; jj < eqt ; jj++){
					model.evolveEQ(jj,false);
					if(jj%1000 == 0){
						params.set("Status", (jj-eqt));
						Job.animate();
					}
				}

				// simulate model for data
				params.set("Mode", "Simulating");
				params.set("\u03D5",phin);
				Job.animate();
				model.setBname(params.sget("Data File")+pfmt.format(100*phin));
				for (int jj = 0 ; jj < simt ; jj++){
					model.evolveEQ(jj,true);
					if(jj%1000 == 0){
						params.set("Status", (jj));
						Job.animate();
					}
				}

				if((simt-1)%damage2Dfast.dlength != 0) model.writeData(simt);
				phin -= dphi;
				phin = (double)(Math.round(100*phin))/100;
			}
			params.set("Mode", "Done");
			Job.animate();
			
		}
		else{
			while (phin > 0){

				params.set("Mode", "Freeze");
				params.set("\u03D5",1-phin);
				Parameters dparams = dummyParamUtil.ofcParams(params);
				model  = new damage2Dfast(dparams);
				if(phin == 1) model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".log",params,model);	

				// equilibrate model			
				params.set("Mode", "Equilibrating");
				params.set("\u03D5",phin);
				Job.animate();
				for (int jj = 0 ; jj < eqt ; jj++){
					model.evolve(jj,false);
					if(jj%1000 == 0){
						params.set("Status", (jj-eqt));
						Job.animate();
					}
				}

				// simulate model for data
				params.set("Mode", "Simulating");
				params.set("\u03D5",phin);
				Job.animate();
				model.setBname(params.sget("Data File")+pfmt.format(100*phin));
				for (int jj = 0 ; jj < simt ; jj++){
					//model.evolve(jj,true);	
					// temp change
					model.evolve(jj,false);	
					hes.accum(model.getGR());
					if(jj%1000 == 0){
						params.set("Status", (jj));
						Job.animate();
					}
				}
				// temp change
				//if((simt-1)%damage2Dfast.dlength != 0) model.writeData(simt);
				printHist(phin);
				phin -= dphi;
				phin = (double)(Math.round(100*phin))/100;
			}
			params.set("Mode", "Done");
			Job.animate();
		}
		
		
		
		return;
	}
	
	public void printHist(double phi){
		PrintUtil.printHistToFile(model.getOutdir()+File.separator+model.getBname()+"_"+pfmt.format(100*phi)+".txt",hes);
		hes.clear();
	}

	public void animate() {
		
		return;
	}

	public void clear() {

		return;
	}
	
	
}
