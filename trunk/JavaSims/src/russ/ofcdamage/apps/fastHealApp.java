package russ.ofcdamage.apps;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import russ.ofcdamage2.damage2Dfast;
import scikit.graphics.ColorGradient;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.heal2Dfast;
import chris.util.DirUtil;

public class fastHealApp extends Simulation{

	private int simt, eqt, L;
	private double phi;
	@SuppressWarnings("unused")
	private boolean draw, record;
	private String PicDir;
	private Grid gridS, gridD;
	private ColorGradient cg = new ColorGradient();
	private DecimalFormat fmt = new DecimalFormat("0000000");
	private DecimalFormat fmtd = new DecimalFormat("0.0000");
	private heal2Dfast model;
	
	public static void main(String[] args) {
		new Control(new fastHealApp(), "Damage Parameters");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int) 10);
		params.add("Lattice Size", (int) 256);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Open"));
		params.add("Equil Time", 500000);
		params.add("Heal Time",(int) 1);
		params.add("HT Width", (int) 0);
		params.add("Failure Stress (\u03C3_f)", 2.);
		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.2);
		params.add("Dissipation (\u03B1)", 0.1);
		params.add("\u03B1 width", 0.);
		params.add("Animate", new ChoiceValue("Off","Draw","Record"));
		params.add("Status");
		params.add("Order Parameter (\u03A6)");
		
		gridS = new Grid("Stress");
		gridD = new Grid("Damage");
		c.frameTogether("Damage Model",gridS,gridD);
		
		return;
	}
	
	public void run() {
				
		// Setup model
		params.set("Status", "Intializing");
		params.set("Order Parameter (\u03A6)", "-");
		Job.animate();
		L      = params.iget("Lattice Size");
		model  = new heal2Dfast(params);
		model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".txt",params);	
		eqt    = params.iget("Equil Time");
		simt   = 0;
		phi    = 1;
		params.set("Order Parameter (\u03A6)", fmtd.format(phi));
		params.set("Status", "Ready");
		Job.animate();
		draw = false;
		
		
		// Equilibrate the system
		for (int jj = 0 ; jj < eqt ; jj++){
			model.evolve(jj,false);
			if(jj%500 == 0){
				params.set("Status", (jj-eqt));
				Job.animate();
			}
		}
		
		
		// Setup I/O
		if(params.sget("Animate").equals("Draw")){
			setupIO(false);
		}
		else if (params.sget("Animate").equals("Record")){
			setupIO(true);
		}
		
		
		// Simulate the model with healing damage
		while(phi > 0){
			phi = model.evolveH(simt,true);
			if(simt%500 == 0) {
				params.set("Status", simt);
				params.set("Order Parameter (\u03A6)", fmtd.format(phi));
				//Job.animate();
			}
//			if(record) takePicture(simt);
			simt++;
		}
		if((simt-1)%damage2Dfast.dlength != 0) model.writeData(simt);
		params.set("Status", "Done");
		Job.animate();
		
		return;
	}

	public void animate() {
		
		if(!draw) return;
		
		gridS.registerData(L, L, model.getStress());
		gridD.registerData(L, L, model.getDorA());
		takePicture(simt,model.getGR());
		return;
	}

	public void clear() {
			
		gridS.clear();
		gridD.clear();
		return;
	}

	public void setupIO(boolean rec){
		
		double sr0, dsr, sf0, dsf, clr, N;
		
		draw = true;
		N    = L*L;
		
		if(rec){
			PicDir = params.sget("Data Directory") + "/GridPics"; 
			DirUtil.MkDir(PicDir);
			record = true;
		}
		
		gridS.setAutoScale(false);
		gridD.setAutoScale(false);
		
		sf0 = params.fget("Failure Stress (\u03C3_f)");
		dsf = params.fget("\u03C3_f width");
		sr0 = params.fget("Residual Stress (\u03C3_r)");
		dsr = params.fget("\u03C3_r width");

		gridD.setScale(0,1);
		gridS.setScale(sr0-dsr,sf0+dsf);
		for (int jj = 0 ; jj <= N ; jj++){
			clr = sr0-dsr + jj*(sf0+dsf-sr0+dsr)/N;
			cg.getColor(clr,sr0-dsr,sf0+dsf);
			gridS.setColors(cg);
		}
		
		return;
	}	
	
	public void takePicture(int time){
		
		String SaveAs = PicDir + File.separator + gridS.getTitle()+fmt.format(time)+".png";
		try {
			ImageIO.write(gridS.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		SaveAs = PicDir + File.separator + gridD.getTitle()+fmt.format(time)+".png";
		try {
			ImageIO.write(gridD.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
		return;
	}
	
	public void takePicture(int time, int iter){
		
		String SaveAs = PicDir + File.separator + gridS.getTitle()+fmt.format(time)+"_"+fmt.format(iter)+".png";
		try {
			ImageIO.write(gridS.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		SaveAs = PicDir + File.separator + gridD.getTitle()+fmt.format(time)+"_"+fmt.format(iter)+".png";
		try {
			ImageIO.write(gridD.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
		return;
	}
	
}
