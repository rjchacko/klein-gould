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
import chris.ofcdamage.damage2Dfast;
import chris.util.CloneUtil;
import chris.util.DirUtil;
import chris.util.MathUtil;

public class fastSubsetDamageApp extends Simulation{
	
	private int simt, eqt, dmt, Ndead, L, N;
	private boolean draw, record, doa[], c1;
	private String PicDir;
	private Grid gridS, gridD;
	//private ColorGradient cg = new ColorGradient();
	private DecimalFormat fmt = new DecimalFormat("0000000");
	private damage2Dfast model, restore;
	
	public static void main(String[] args) {
		new Control(new fastSubsetDamageApp(), "Damage Parameters");
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
		params.add("Number of Lives",(int) 5);
		params.add("NL width", (int) 0);
		params.add("Failure Stress (\u03C3_f)", 2.);
		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.05);
		params.add("Dissipation (\u03B1)", 0.05);
		params.add("\u03B1 width", 0.);
		params.add("Animate", new ChoiceValue("Off","Draw","Record"));
		params.add("Status");
		params.add("Dead Sites");
		params.add("Mode");
		params.set("Mode","Cycle 1");


		
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
		doa    = new boolean[N];
		model  = new damage2Dfast(params);
		model.PrintParams(model.getOutdir()+File.separator+"Params_"+model.getBname()+".txt",params);	
		eqt    = params.iget("Equil Time");
		simt   = params.iget("Sim Time");
		dmt    = simt;
		Ndead  = 0;
		c1     = true;
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
		
		// clone system and save it for second run through (state 1)
		restore = CloneUtil.cloneDamage(params, model, true);
		
		// Simulate the model without damage
		for (int jj = 0 ; jj < simt ; jj++){
			model.evolve(jj,false);
			if(jj%500 == 0){
				params.set("Status", jj);
			}
			Job.animate();
			if(record) takePicture(jj);
		}

		// Simulate the model with damage
		while(Ndead < N){
			doa   = model.getDorAbool();
			Ndead = model.evolveD(dmt,false);
			if(dmt%500 == 0){
				params.set("Status", (dmt));
			}
			params.set("Dead Sites", Ndead);
			Job.animate();
			if(record) takePicture(dmt);
			dmt++;
		}
		
		// re-run simulation
		c1    = false;
		model = null; // clear from memory
		dmt   = simt;
		Ndead = 0;
		int nss = 0;
		for (int jj = 0 ; jj < N ; jj++){
			nss += MathUtil.bool2bin(doa[jj]);
		}
		restore.setNss(nss);
		params.set("Dead Sites", Ndead);
		params.set("Mode","Cycle 2");
		Job.animate();

		// Simulate the model without damage
		for (int jj = 0 ; jj < simt ; jj++){
			restore.evolve(jj,doa);
			if(jj%500 == 0){
				params.set("Status", jj);
			}
			Job.animate();
			if(record) takePicture(jj);
		}


		
		// Simulate the model with damage
		while(Ndead < N){
			Ndead = restore.evolveD(dmt,doa);
			if(dmt%500 == 0){
				params.set("Status", (dmt));
			}
			params.set("Dead Sites", Ndead);
			Job.animate();
			if(record) takePicture(dmt);
			dmt++;
		}

		if((dmt-1)%damage2Dfast.dlength != 0) restore.writeData(dmt);
		
		params.set("Status", "Done");
		Job.animate();
		
		return;
	}

	public void animate() {
		
		if(!draw) return;
		if(c1){
			gridS.registerData(L, L, model.getStress());
		}
		else{
			gridS.registerData(L, L, restore.getStress());
			}
		//gridD.registerData(L, L, model.getDorA());

		return;
	}

	public void clear() {
			
		gridS.clear();
		gridD.clear();
		return;
	}

	public void setupIO(boolean rec){
		
		double sr0, dsr, sf0, dsf;
		
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
//		for (int jj = 0 ; jj <= N ; jj++){
//			clr = sr0-dsr + jj*(sf0+dsf-sr0+dsr)/N;
//			cg.getColor(clr,sr0-dsr,sf0+dsf);
//			gridS.setColors(cg);
//		}
		
		return;
	}	
	
	public void takePicture(int time){
		
		String SaveAs = PicDir + File.separator + gridS.getTitle()+fmt.format(time)+".png";
		try {
			ImageIO.write(gridS.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
//		SaveAs = PicDir + File.separator + gridD.getTitle()+fmt.format(time)+".png";
//		try {
//			ImageIO.write(gridD.getImage(), "png", new File(SaveAs));
//		} catch (IOException e) {
//			System.err.println("Error in Writing File" + SaveAs);
//		}
		
		return;
	}
	
}
