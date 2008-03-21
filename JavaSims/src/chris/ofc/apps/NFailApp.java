package chris.ofc.apps;


import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.ofc.NfailDamage2D;
import chris.util.PrintUtil;

public class NFailApp extends Simulation {

	Grid grid1 = new Grid ("Stress Lattice");
	Grid grid2 = new Grid ("Failed Sites");

	NfailDamage2D model;
	double ScMax, rgyr;
	DecimalFormat fmt = new DecimalFormat("0000000");
	DecimalFormat fmts = new DecimalFormat("0000");
	
	ColorPalette palette1;
	ColorGradient smooth;
	
	int EEchk = 0;
	

	public static void main(String[] args) {
		new Control(new NFailApp(), "OFC Model");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/CurrentSemester/Research/Data/"));
		params.add("Random Seed",0);
		params.add("Animation", new ChoiceValue("On","Off"));
		//params.addm("Auto Scale", new ChoiceValue("Yes", "No"));
		params.add("Lattice Size",1<<9);
		params.add("Number of Lives",1);
		params.add("Life Style", new ChoiceValue("Constant","Flat","Gaussian"));
		params.add("Nlives Width",0.1);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Bordered"));
		params.add("Stress Distribution", new ChoiceValue("Flat","Hammer Blow"));
		params.add("Hammer Size",1);	
		params.add("Critical Stress (\u03C3_c)",4.0);
		params.add("\u03C3_c Noise", new ChoiceValue("Off","On"));	
		params.add("\u03C3_c width",Math.sqrt(Math.sqrt(0.4)));
		params.add("Residual Stress (\u03C3_r)",2.0);
		params.add("\u03C3_r Noise", new ChoiceValue("Off","On"));
		params.add("\u03C3_r width",Math.sqrt(Math.sqrt(2)));
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond"));
		params.add("Interaction Radius (R)",(int)(50));
		params.add("Minimum Interaction Radius (r)",0);
		params.add("Dissipation (\u03B1)",new DoubleValue(0.2,0,1));
		params.add("\u03B1 Noise", new ChoiceValue("On","Off"));
		params.add("\u03B1 Width", 0.05);
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("Number of Resets");
		params.add("Number of Showers");
			
		c.frameTogether("Stress Lattice in an OFC Model with Damage", grid1, grid2);
		
	}
	
	public void animate() {
						
		params.set("Number of Resets",model.time);
		params.set("Number of Showers",model.showernumber);
		
		if (model.showering && model.ShowGrid){
		
			int[] foo = new int[model.N];
	
				
			for (int i=0 ; i<model.N ; i++){
				smooth.getColor(model.stress[i],-2,ScMax);
				foo[i]=model.alive[i];
			}
				
			grid1.setColors(smooth);
			grid1.registerData(model.L,model.L,model.stress);
			grid2.registerData(model.L, model.L, foo);
				
			if (params.sget("Record").equals("On") && model.ShowGrid) GridPicture("both"); 
		
		}

	}

	public void clear() {

		grid1.clear();
		grid2.clear();
	}

	public void run() {
		
		model = new NfailDamage2D(params);
		
		String anmt = params.sget("Animation");
		
		if (anmt.equals("On")){
			model.ShowGrid=true;
		}
		else{
			model.ShowGrid=false;
		}
		
		PrintUtil.printlnToFile(model.outdir+File.separator+"Params.txt",params.toString());
				
		model.Initialize(params.sget("Stress Distribution"));
		
		ScMax=model.GetMax(model.Sc);		

		// Set up color scheme
		
		palette1 = new ColorPalette();
		smooth   = new ColorGradient();
		grid2.setColors(palette1);
		
		int max = model.GetMax(model.alive);

		for (int i = 0 ; i <= max ; i++){
			palette1.setColor(i,smooth.getColor(i, 0, max));
		}
		
		// Set up file
			
		PrintUtil.printlnToFile(model.outfile,"Time","N_avlnchs","N_dead","Rgyr","Omega","<FS_stress>","rho_FS");
		
		while(!(model.crack)) {
			
			model.Avalanche();

			TakeData();
			
		}
		
	}
	
	public void TakeData(){
	
		// FIX THIS EVERYWHERE OR BETTER YET
		// MOVE IT TO NFAILDAMAGE2D
		
		int[] LS = model.LiveSites(); 
		int LSlength = LS.length;
		if(LSlength>0){
			rgyr=model.radiusGyration(LS[model.rand.nextInt(LSlength)]);
		}
		else{
			rgyr=0;
		}
		
		PrintUtil.printlnToFile(model.outfile,model.time,model.Nshowers,model.NdeadS,rgyr,model.EFmetric(),model.GetAve(model.SonFS,model.SonFSindex),model.DaboutFS);
				
		return;
	}
	
	public void GridPicture(String which){
		
		
		String imgName1 = model.PicDir+File.separator+"Alive-"+fmt.format(model.time)+"-"+fmts.format(model.showernumber)+".png";
		String imgName2 = model.PicDir+File.separator+"Stress-"+fmt.format(model.time)+"-"+fmts.format(model.showernumber)+".png";
		
		if(which.equals("alive")){
			
			try {
				ImageIO.write(grid1.getImage(), "png", new File(imgName2));
			} catch (IOException e) {
				System.err.println("Error in Writing File" + imgName2);
			}
		}
		else if (which.equals("stress")){
			
			try {
				ImageIO.write(grid2.getImage(), "png", new File(imgName1));
			} catch (IOException e) {
				System.err.println("Error in Writing File" + imgName1);
			}
		}
		else if (which.equals("both")){
		
			try {
				ImageIO.write(grid2.getImage(), "png", new File(imgName1));
				ImageIO.write(grid1.getImage(), "png", new File(imgName2));
			} catch (IOException e) {
				System.err.println("Error in Writing File" + imgName1 +"and" + imgName2);
			}
		}
		else {
			System.err.println("Error. String " + which + " not found.");
		}
		
		return;
	}

}
