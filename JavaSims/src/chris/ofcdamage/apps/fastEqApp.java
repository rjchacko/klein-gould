package chris.ofcdamage.apps;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.ofc2Dfast;

public class fastEqApp extends Simulation{

	private int simt, eqt;
	private ofc2Dfast model;
	private int eventID;
	private boolean failed[];
	private static DecimalFormat fmtEid = new DecimalFormat("00");

	public static void main(String[] args) {
		new Control(new fastEqApp(), "OFC Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int) 10);
		params.add("Lattice Size", (int) 256);
		params.add("Boundary Condtions", new ChoiceValue("Open","Periodic"));
		params.add("Equil Time", 100000);
		params.add("Sim Time", 100000);
		params.add("Failure Stress (\u03C3_f)", 2.);
		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.025);
		params.add("Dissipation (\u03B1)", 0.025);
		params.add("\u03B1 width", 0.);
		params.add("Status");
		
	}
	
	public void run() {
		
		// Setup model
		eventID = 0;
		params.set("Status", "Intializing");
		Job.animate();
		model = new ofc2Dfast(params);
		model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".log",params);	
		failed = new boolean[model.getN()];
		eqt    = params.iget("Equil Time");
		simt   = params.iget("Sim Time");
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
			model.evolve(jj,true);
			if(jj%500 == 0){
				params.set("Status", jj);
				Job.animate();
			}
		}

		if((simt-1)%ofc2Dfast.dlength != 0) model.writeData(simt);
		
//		PrintUtil.printHistToFile(model.getOutdir()+File.separator+"hfail_"+model.getBname()+".txt", model.hfail);
//		PrintUtil.printHistToFile(model.getOutdir()+File.separator+"hstress_"+model.getBname()+".txt", model.hstrs);
//		
		params.set("Status", "Done");
		Job.animate();
		
		return;
	}

	public void animate() {
		
		return;
	}

	public void clear() {
		
		return;
	}
	
	public void writeEvent(int time){
		failed = model.getLastShower();
		
		try{
			File file = new File(model.getOutdir()+File.separator+model.getBname()+"_event_"+fmtEid.format(eventID)+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println("Time = " + time + ", Events Size = " + model.getGR());
			for(int jj = 0 ; jj < model.getN(); jj++){
				if(failed[jj])
					pw.println(jj);
			}
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		return;
	}

	
}
