package chris.ofcdamage.apps;

import java.io.File;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.FileValue;
import scikit.jobs.params.Parameters;
import chris.ofcdamage.ofc2Dfast;
import chris.util.PrintUtil;
import chris.util.ReadInUtil;

public class generateCatalogueApp extends Simulation{

	private ofc2Dfast model;
	private Parameters p2;
	
	public static void main(String[] args) {
		new Control(new generateCatalogueApp(), "OFC Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Parameter File",new FileValue("/Users/cserino/Documents/Catalogue/Params_run1.log"));
//		params.add("Data File");
//		params.add("Random Seed");
//		params.add("Interaction Shape");
//		params.add("Interaction Radius (R)");
//		params.add("Lattice Size");
//		params.add("Boundary Condtions");
//		params.add("Equil Time");
		params.add("Sim Time", (int)(1e6));
//		params.add("Failure Stress (\u03C3_f)");
//		params.add("\u03C3_f width");
//		params.add("Residual Stress (\u03C3_r)");
//		params.add("\u03C3_r width");
		params.add("Dissipation (\u03B1)");
//		params.add("\u03B1 width");
		params.add("Status");
	}

	public void run() {

		// read in params from file 
		params.set("Status", "Configuring");
		Job.animate();
		
		int tsim       = params.iget("Sim Time");
		ReadInUtil riu = new ReadInUtil(params.sget("Parameter File"));
		p2             = riu.getOFCparams();
		double[] av    = new double[]{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
		riu            = null;
		
		while(true){
			for(double alpha : av){
				params.set("Status", "Intializing");
				params.set("Dissipation (\u03B1)",alpha);
				Job.animate();
				model = new ofc2Dfast(p2);
				model.setClength(tsim);
				// read in stress from file
				riu   = new ReadInUtil(model.getOutdir()+File.separator+model.getBname()+"_Stress_"+(int)(100*alpha)+".txt");
				model.setStress(riu.getData(0));				
				
				// simulate <i>tsim</i> time steps
				for(int tt = 0 ; tt < tsim ; tt++){
					model.evolve(tt, false);
					if(tt%1000 == 0){
						params.set("Status",tt);
						Job.animate();
					}
				}

				// append data to data file
				model.appendCdata(model.getOutdir()+File.separator+model.getBname()+"_Catalogue_"+(int)(100*alpha)+".txt",tsim);
		
				// replace old stress file with new stress file
				PrintUtil.overwriteFile(model.getOutdir()+File.separator+model.getBname()+"_Stress_"+(int)(100*alpha)+".txt",model.getStress());
				
				// print summary of events in log file
				editLogFile(tsim, alpha);
				
				// update seed
				p2.set("Random Seed", p2.iget("Random Seed")+1);
			}
		}
	}
	
	private void editLogFile(int tsim, double alpha){

		DateFormat dateFormat = new SimpleDateFormat("dd MM yyyy HH:mm:ss");
        Date date = new Date();	
		PrintUtil.printlnToFile(model.getOutdir()+File.separator+model.getBname()+"_Catalogue.log","Added "+tsim+" events to the alpha = "+alpha+" catalogue on"+dateFormat.format(date));
		return;
	}
	
	public void animate() {return;}
	public void clear() {return;}
	
}
