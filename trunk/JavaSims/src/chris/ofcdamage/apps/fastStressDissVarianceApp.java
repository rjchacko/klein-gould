package chris.ofcdamage.apps;

import java.io.File;
import java.text.DecimalFormat;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.ofc2Dfast;
import chris.util.PrintUtil;

public class fastStressDissVarianceApp extends Simulation{

	private int simt, eqt, GR[];
	private ofc2Dfast model;
	private double dstress[];
	private static DecimalFormat fmtEid = new DecimalFormat("00");

	public static void main(String[] args) {
		new Control(new fastStressDissVarianceApp(), "OFC Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
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
		params.add("\u03C3_r width", 0.025);
		params.add("Dissipation (\u03B1)", 0.025);
		params.add("\u03B1 width", 0.);
		params.add("Status");
		
//		params.add("Data Directory",new DirectoryValue("~"));
//		params.add("Data File", "default");
//		params.add("Random Seed", (int) 0);
//		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
//		params.add("Interaction Radius (R)", (int) 45);
//		params.add("Lattice Size", (int) 768);
//		params.add("Boundary Condtions", new ChoiceValue("Periodic","Open"));
//		params.add("Equil Time", 500000);
//		params.add("Sim Time", 500000);
//		params.add("Failure Stress (\u03C3_f)", 1.);
//		params.add("\u03C3_f width", 0.);
//		params.add("Residual Stress (\u03C3_r)", 0.25);
//		params.add("\u03C3_r width", 0.25);
//		params.add("Dissipation (\u03B1)", 0.);
//		params.add("\u03B1 width", 0.);
//		params.add("Status");
	}
	
	public void run() {
		
		// Setup model
		//double alpha[] = new double[]{0.010000000 ,0.012915497 ,0.016681005 ,0.021544347 ,0.027825594 ,0.035938137 ,0.046415888 ,0.059948425 ,0.077426368 ,0.100000000};
		double a;
		int Na = 30;
		double alpha[] = new double[Na];
		int cycle[]    = new int[Na];
		for (int jj = 0 ; jj < Na ; jj++){
			alpha[jj] = Math.pow(10,-3+jj*(Math.log10(5)+2)/(Na-1)); // Na numbers log distr on [10^-3 , 10^{log(2)-1} = 2 x 10^-1]
			cycle[jj] = jj;
		}
		
		for (int jj = 0 ; jj < Na ; jj++){
			a = alpha[jj];
			params.set("Dissipation (\u03B1)", a);
			params.set("Random Seed", params.iget("Random Seed") + 1);
			params.set("Status", "Intializing");
			Job.animate();
			model = new ofc2Dfast(params);
			if(a == alpha[0]){ 
				model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".log",params,this);	
				PrintUtil.printVectorsToFile(model.getOutdir()+File.separator+model.getBname()+"_AlphaValues.txt",cycle,alpha);
			}
			eqt     = params.iget("Equil Time");
			simt    = params.iget("Sim Time");
			dstress = new double[simt]; 
			GR      = new int[simt];
			params.set("Status", "Ready");
			Job.animate();

			// Equilibrate the system
			for (int kk = 0 ; kk < eqt ; kk++){
				model.evolve(kk,false);
				if(kk%500 == 0){
					params.set("Status", (kk-eqt));
					Job.animate();
				}
			}

			// Simulate the model without damage
			for (int kk = 0 ; kk < simt ; kk++){
				model.evolve(kk,false);
				dstress[kk] = model.getStressDissipated();
				GR[kk] = model.getGR();
				if(kk%500 == 0){
					params.set("Status", kk);
					Job.animate();
				}
			}

			PrintUtil.printVectorsToFile(model.getOutdir()+File.separator+model.getBname()+"_"+fmtEid.format(jj)+".txt", GR, dstress);
			PrintUtil.printHistToFile(model.getOutdir()+File.separator+model.getBname()+"_"+fmtEid.format(jj)+".his", model.getHfail());

			params.set("Status", "Done");
			Job.animate();
			model = null;

		}
		
		return;
	}

	public void animate() {
		
		return;
	}

	public void clear() {
		
		return;
	}
	
}
