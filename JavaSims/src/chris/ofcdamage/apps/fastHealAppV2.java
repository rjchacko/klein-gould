package chris.ofcdamage.apps;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import scikit.graphics.ColorGradient;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.damage2Dfast;
import chris.ofcdamage.heal2DfastV2;
import chris.util.DirUtil;

public class fastHealAppV2 extends Simulation{

	private int simt, eqt, dmt, Ndead, L, N;
	private boolean draw, record;
	private String PicDir;
	private Grid gridS, gridD;
	private ColorGradient cg = new ColorGradient();
	private DecimalFormat fmt = new DecimalFormat("0000000");
	private heal2DfastV2 model;
	
	public static void main(String[] args) {
		new Control(new fastHealAppV2(), "Damage Parameters");
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
		params.add("Sim Time", 500000);
		params.add("Number of Lives",(int) 1);
		params.add("NL width", (int) 0);
		params.add("Heal Time",(int) 1);
		params.add("HT Width", (int) 0);
		params.add("Failure Stress (\u03C3_f)", 2.);
		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.05);
		params.add("Dissipation (\u03B1)", 0.05);
		params.add("\u03B1 width", 0.);
		params.add("Animate", new ChoiceValue("Draw","Off","Record"));
		params.add("Status");
		params.add("Dead Sites");
		
//		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/test"));
//		params.add("Data File", "default");
//		params.add("Random Seed", (int) 0);
//		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
//		params.add("Interaction Radius (R)", (int) 2);
//		params.add("Lattice Size", (int) 100);
//		params.add("Boundary Condtions", new ChoiceValue("Periodic","Open"));
//		params.add("Equil Time", 0);
//		params.add("Sim Time", 0);
//		params.add("Number of Lives",(int) 1);
//		params.add("NL width", (int) 0);
//		params.add("Heal Time",(int) 1);
//		params.add("HT Width", (int) 0);
//		params.add("Failure Stress (\u03C3_f)", 2.);
//		params.add("\u03C3_f width", 0.);
//		params.add("Residual Stress (\u03C3_r)", 1.);
//		params.add("\u03C3_r width", 0.2);
//		params.add("Dissipation (\u03B1)", 0.2);
//		params.add("\u03B1 width", 0.);
//		params.add("Animate", new ChoiceValue("Draw","Off","Record"));
//		params.add("Status");
//		params.add("Dead Sites");
		
		
		
		gridS = new Grid("Stress");
		gridD = new Grid("Damage");
		c.frameTogether("Damage Model",gridS,gridD);
		
		return;
	}
	
	public void run() {
				
		// Setup model
		params.set("Status", "Intializing");
		params.set("Dead Sites", "-");
		Job.animate();
		L      = params.iget("Lattice Size");
		N      = L*L;
		model  = new heal2DfastV2(params);
		model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".txt",params);	
		eqt    = params.iget("Equil Time");
		simt   = params.iget("Sim Time");
		dmt    = simt;
		Ndead  = 0;
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
		
		// Simulate the model without damage
		for (int jj = 0 ; jj < simt ; jj++){
			model.evolve(jj,true);
			if(jj%500 == 0){
				params.set("Status", jj);
				Job.animate();
				if(record) takePicture(jj);
			}
		}
		
		// Simulate the model with healing damage
		while(Ndead < N){
			Ndead = model.evolveH(dmt,true);
			if(dmt%500 == 0){
				params.set("Status", (dmt));
			}
			params.set("Dead Sites", Ndead);
			Job.animate();
			if(record) takePicture(dmt);
			dmt++;
		}

		if((dmt-1)%damage2Dfast.dlength != 0) model.writeData(dmt);
		
		params.set("Status", "Done");
		Job.animate();
		
		return;
	}

	public void animate() {
		
		if(!draw) return;
		
		gridS.registerData(L, L, model.getStress());
		gridD.registerData(L, L, model.getDorA());

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
	
}
