package chris.tests;



import java.awt.Color;

import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.ofc.NfailDamage2D;
import chris.util.LatticeNeighbors;


public class GridTests extends Simulation{

	Grid grid1 = new Grid ("Sites");

	NfailDamage2D model;

	ColorPalette palette1;

	LatticeNeighbors LN;
	
	double displaycounter;
	
	int x0,y0,i0;
	int[] foo;
	
	public static void main(String[] args) {
		new Control(new GridTests(), "OFC Model");
	}
	
	public void load(Control c) {		

		/**
		 * 	NEED THESE
		 */
		
		params.add("Input Directory",new DirectoryValue("/Users/cserino/CurrentSemester/Research/Data/"));
		params.add("Data Directory",new DirectoryValue("/Users/cserino/CurrentSemester/Research/"));
		params.add("Random Seed",0);
		params.add("Animation", new ChoiceValue("On","Off"));
		params.add("Lattice Size",500);
		params.add("Number of Lives",1);
		params.add("Life Style", new ChoiceValue("Constant","Flat","Gaussian"));
		params.add("Nlives Width",0.1);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Bordered"));	
		params.add("Critical Stress (\u03C3_c)",4.0);
		params.add("\u03C3_c Noise", new ChoiceValue("Off","On"));	
		params.add("\u03C3_c width",Math.sqrt(Math.sqrt(0.4)));
		params.add("Residual Stress (\u03C3_r)",2.0);
		params.add("\u03C3_r Noise", new ChoiceValue("Off","On"));
		params.add("\u03C3_r width",Math.sqrt(Math.sqrt(2)));
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","Ellipse"));
		params.add("Eccentricity", new DoubleValue(0,0,1));
		params.add("Interaction Radius (R)",(int)(50));
		params.add("Minimum Interaction Radius (r)",30);
		params.add("Dissipation (\u03B1)",new DoubleValue(0.2,0,1));
		params.add("\u03B1 Noise", new ChoiceValue("On","Off"));
		params.add("\u03B1 Width", 0.05);
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("Number of Resets");
		params.add("Number of Showers");
			
		/**
		 * 	TEST PARAMETERS
		 */
		
		params.add("x0", 250);
		params.add("y0", 250);
		
		c.frame(grid1);
	}
	
	public void animate() {
		
		if (params.sget("Animation").equals("On")) grid1.registerData(model.L,model.L,foo);
		params.set("Number of Resets",displaycounter);

	}

	public void clear() {
		
		grid1.clear();

	}

	public void run() {

		// Get functionality of NfailDamage2D
		
		model = new NfailDamage2D(params);
		model.Initialize("Flat");
		
		// Set up Lattice Neighbors
		
		if (params.sget("Boundary Condtions").equals("Bordered")){
			
			if(params.sget("Interaction Shape").equals("Circle")){
				LN = new LatticeNeighbors(model.L,model.L,model.rmin,model.R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Circle);
			}
			else if(params.sget("Interaction Shape").equals("Square")){
				LN = new LatticeNeighbors(model.L,model.L,model.rmin,model.R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Square);
			}
			else if(params.sget("Interaction Shape").equals("Diamond")){
				LN = new LatticeNeighbors(model.L,model.L,model.rmin,model.R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Diamond);
			}
			else if(params.sget("Interaction Shape").equals("Ellipse")){
				LN = new LatticeNeighbors(model.L,model.L,params.fget("Eccentricity"),model.rmin,model.R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Ellipse);
			}
		}
		else{
			
			if(params.sget("Interaction Shape").equals("Circle")){
				LN = new LatticeNeighbors(model.L,model.L,model.rmin,model.R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Circle);
			}
			else if(params.sget("Interaction Shape").equals("Square")){
				LN = new LatticeNeighbors(model.L,model.L,model.rmin,model.R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Square);
			}
			else if(params.sget("Interaction Shape").equals("Diamond")){
				LN = new LatticeNeighbors(model.L,model.L,model.rmin,model.R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Diamond);
			}
			else if(params.sget("Interaction Shape").equals("Ellipse")){
				LN = new LatticeNeighbors(model.L,model.L,params.fget("Eccentricity"),model.rmin,model.R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Ellipse);
			}
		}
		
		// Set up simple color scheme
		
		palette1 = new ColorPalette();
		palette1.setColor(0,Color.BLACK);
		palette1.setColor(1,Color.WHITE);
		grid1.setColors(palette1);
		
		
		// Get x0 and y0
		x0 = params.iget("x0");
		y0 = params.iget("y0");
		i0 = y0*model.L + x0;
		
		// Initialize foo
		
		foo = new int[model.N];
		
		// Set time = 0
		
		displaycounter = 0;
		
		/**
		 *	Test the ellipse method for lattice neighbors
		 */
		
//		int[] thenbs = LN.get(i0);
//		
//		for (int ii = 0 ; ii < model.N ; ii++){
//			foo[ii]=0;
//		}
//		
//		for (int ii = 0 ; ii < thenbs.length; ii++){
//			foo[thenbs[ii]]=1;
//		}
//		
//		Job.animate();
		
		/**
		 * 			E-N-D
		 * 	
		 *	Test the ellipse method for lattice neighbors
		 */
		
		////////////////////////////////////////////////////////////////////
		
		/**
		 * 	Test getStressLines(int center, int Nlines) method	
		 */
		
//		palette1.setColor(2,Color.RED);
//		grid1.setColors(palette1);
//		
//			while(true){
//				
//				displaycounter++;
//				
//				int[] thenbs = LN.get(i0);
//				
//				for (int ii = 0 ; ii < model.N ; ii++){
//					foo[ii]=0;
//				}
//				
//				for (int ii = 0 ; ii < thenbs.length; ii++){
//					foo[thenbs[ii]]=1;
//				}	
//				
//				int[] SLs = model.getStressLines(i0, 1);
//				
//				for (int ii = 0 ; ii < SLs.length ; ii++){
//					foo[SLs[ii]]=2;
//				}
//				
//				
//				Job.animate();
//			}
//		
	
	/**
	 * 			E-N-D
	 * 
	 * 	Test getStressLines(int center, int Nlines) method	
	 */

		
		////////////////////////////////////////////////////////////////////

		
	/**
	 *	Test method printArrayToFile(String fileName, int[] array, int m, int n) 
	 */
		

//	for (int ii = 0 ; ii < model.N ; ii++){
//		foo[ii]=ii+1;
//	}
//	
//	PrintUtil.printArrayToFile(model.outdir+File.separator+"ArrayTest.txt", foo, model.L, model.L);
//	
//	Job.animate();
	
	/**
	 * 			E-N-D
	 * 
	 *	Test method printArrayToFile(String fileName, int[] array, int m, int n) 
	 */
	
	
	////////////////////////////////////////////////////////////////////

	
	/**
	 *	Test method Initialize(Parameters prms)
	 */
		
	
	model.Initialize(params);
	
	displaycounter = 77;
	
	Job.animate();
		
	/**
	 * 			E-N-D
	 * 	Test method Initialize(Parameters prms)
	 */
		
		
		
		
	}	
}




