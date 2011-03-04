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
import chris.util.PrintUtil;

public class fastEqApp extends Simulation{

	private int simt, eqt;
	private ofc2Dfast model;
	private int eventID;
	private boolean failed[];
	private int grts[];
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
		double alpha[] = new double[]{0.010000000 ,0.012915497 ,0.016681005 ,0.021544347 ,0.027825594 ,0.035938137 ,0.046415888 ,0.059948425 ,0.077426368 ,0.100000000};
		
		for (int a = 0; a < 10 ; a++){
			
			params.set("Dissipation (\u03B1)", alpha[a]);
			params.set("Random Seed", params.iget("Random Seed") + 1);
			eventID = 0;
			params.set("Status", "Intializing");
			Job.animate();
			model = new ofc2Dfast(params);
			if(a == 0) 
				model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".log",params,this);	
			failed = new boolean[model.getN()];
			eqt    = params.iget("Equil Time");
			simt   = params.iget("Sim Time");
			grts   = new int[simt];
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
				//model.evolve(jj,true);
				model.evolve(jj,false);
				grts[jj] = model.getGR();
				if(jj%500 == 0){
					params.set("Status", jj);
					Job.animate();
				}
			}

			PrintUtil.printlnToFile(model.getOutdir()+File.separator+model.getBname()+"_"+fmtEid.format(a)+".txt", "alpha = " + fmtEid.format(alpha[a])+" .");
			PrintUtil.printVectorToFile(model.getOutdir()+File.separator+model.getBname()+"_"+fmtEid.format(a)+".txt", grts);
			//if((simt-1)%ofc2Dfast.dlength != 0) model.writeData(simt);

			//		PrintUtil.printHistToFile(model.getOutdir()+File.separator+"hfail_"+model.getBname()+".txt", model.hfail);
			//		PrintUtil.printHistToFile(model.getOutdir()+File.separator+"hstress_"+model.getBname()+".txt", model.hstrs);
			//		
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
