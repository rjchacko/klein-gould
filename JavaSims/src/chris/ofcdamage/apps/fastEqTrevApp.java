package chris.ofcdamage.apps;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import javax.imageio.ImageIO;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.ofc2Dtrevfast;
import chris.util.DirUtil;

public class fastEqTrevApp extends Simulation{

	private int simt, eqt;
	private boolean anim, rec;
	private ofc2Dtrevfast model;
	private String picdir;
	private Grid grid = new Grid ("Stress");
	private DecimalFormat fmt = new DecimalFormat("0000000");
	
	public static void main(String[] args) {
		new Control(new fastEqTrevApp(), "OFC Parameters");
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
		params.add("\u03C3_r width", 0.05);
		params.add("Dissipation (\u03B1)", 0.05);
		params.add("\u03B1 width", 0.);
		params.add("Animation", new ChoiceValue("Off","On","Record"));
		params.add("Status");
		
		c.frame(grid);
	}
	
	public void run() {
		
		// Setup model
		params.set("Status", "Intializing");
		Job.animate();
		model = new ofc2Dtrevfast(params);
		model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".txt",params);	
		eqt   = params.iget("Equil Time");
		simt  = params.iget("Sim Time");
		rec   = params.sget("Animation").equals("Record");
		if(rec) setPicDir();
	    
		params.set("Status", "Ready");
		Job.animate();
		
		// Equilibrate the system
		anim = false;
		for (int jj = 0 ; jj < eqt ; jj++){
			model.evolve(jj,false);
			if(jj%500 == 0){
				params.set("Status", (jj-eqt));
				Job.animate();
			}
		}
		
		// Simulate the model without damage
		anim = !(params.sget("Animation").equals("Off"));
		for (int jj = 0 ; jj < simt ; jj++){
			model.evolve(jj,true);
			if(jj%500 == 0){
				params.set("Status", jj);
			}
			localAnimate(jj);
		}

		if((simt-1)%ofc2Dtrevfast.dlength != 0) model.writeData(simt);
		
		params.set("Status", "Done");
		Job.animate();
		
		return;
	}

	public void localAnimate(int t){
		
		if(anim){
			// invert color map
			double[] tmp = model.getStress();
			for (int jj = 0 ; jj < model.getN() ; jj++){
				tmp[jj] = model.getSmin()+model.getSmax()-tmp[jj];
			}
			grid.registerData(model.getL(), model.getL(), tmp);			
			Job.animate();
			if(rec){
				takePicture(grid, t);
			}
		}
		else{
			Job.animate();
		}
		
	}
	
	public void animate() {
		
		return;
	}

	public void clear() {
		
		grid.clear();
		return;
	}
	
	public void setPicDir(){
		
		picdir = model.getOutdir()+"/Movie";
		//grid.setAutoScale(false);
		//grid.setScale(model.getSmin(),model.getSmax());
		DirUtil.MkDir(picdir);
	}
	
	private void takePicture(Grid grid, int t){

		String SaveAs = picdir + File.separator + grid.getTitle()+fmt.format(t)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}

		return;
	}

	
}
