package chris.tests;



import java.awt.Color;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.old.ofc.clusters.ClustersV4;
import chris.old.ofc.old.NfailDamage2D;
import chris.old.ofc.old.NfailDamage2DV2;
import chris.util.LatticeNeighbors;


public class GridTestsV2 extends Simulation{

	Grid grid1 = new Grid ("Grid 1");
	Grid grid2 = new Grid ("Grid 2");
	Grid grid3 = new Grid ("Grid 3");
	Grid grid4 = new Grid ("Grid 4");
	Grid grid5 = new Grid ("Grid 5");
	
	NfailDamage2DV2 model;
	NfailDamage2D   OLDmodel;
	ClustersV4 dummy; 

	ColorPalette palette1;
	ColorPalette palette2;
	ColorGradient smooth;

	LatticeNeighbors LN;
	
	double displaycounter;
	
	int x0,y0,i0;
	int[] foo;
	double[] display;
	int[] parent, foo2, foo3, foo4;
	
	public static void main(String[] args) {
		new Control(new GridTestsV2(), "OFC Model");
	}
	
	public void load(Control c) {		

		/**
		 * 	NEED THESE
		 */
		
		params.add("Input Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Random Seed",0);
		params.add("Animation", new ChoiceValue("On","Off"));
		params.add("Lattice Size",11);
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
		params.add("Interaction Radius (R)",(int)(2));
		params.add("Minimum Interaction Radius (r)",0);
		params.add("Dissipation (\u03B1)",new DoubleValue(0.2,0,1));
		params.add("\u03B1 Noise", new ChoiceValue("On","Off"));
		params.add("\u03B1 Width", 0.05);
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("Number of Resets");
		params.add("Number of Showers");
			
		/**
		 * 	TEST PARAMETERS
		 */
		
		params.add("x0", 5);
		params.add("y0", 5);
		
		c.frameTogether("Test",grid1,grid2,grid3,grid4,grid5);
	}
	
	public void animate() {
		
		if (params.sget("Animation").equals("On")){
			grid5.registerData(model.L,model.L,foo);
			grid2.registerData(model.L,model.L,foo2);
//			grid3.registerData(model.L,model.L,foo3);
//			grid4.registerData(model.L,model.L,foo4);
//			grid1.registerData(model.L,model.L,display);
		}
		params.set("Number of Resets",displaycounter);

	}

	public void clear() {
		
		grid1.clear();
		grid2.clear();
		grid3.clear();
		grid4.clear();
		grid5.clear();

	}

	public void run() {

		// Get functionality of NfailDamage2D
		
		model = new NfailDamage2DV2(params);
		model.Initialize("Flat");
		
		OLDmodel = new NfailDamage2D(params);
		OLDmodel.Initialize("Flat");
		
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
		
		grid5.setColors(palette1);
		
		// Get x0 and y0
		x0 = params.iget("x0");
		y0 = params.iget("y0");
		i0 = y0*model.L + x0;
		
		// Initialize foo
		
		foo = new int[model.N];
		
		for (int jj = 0 ; jj < model.N ; jj++){
			foo[jj] = 0;
		}
		
		
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
		
	
//	model.Initialize(params);
//	
//	displaycounter = 77;
//	
//	model.CloneSim(params);
//	
//	displaycounter = 7777;
//	
//	params.set("Animation","Off");
//	
//	Job.animate();

		
	/**
	 * 			E-N-D
	 * 	Test method Initialize(Parameters prms)
	 */
		
		////////////////////////////////////////////////////////////////////

		
		/**
		 *	Test class Clusters(Parameters prms)
		 */	
//	
//		
//	    palette1 = new ColorPalette();
//	    palette2 = new ColorPalette();
//	    palette1.setColor(0,Color.WHITE);
//	    palette1.setColor(1,Color.BLACK);
//	    smooth   = new ColorGradient();
//	    grid5.setColors(palette1);
//
//	    Color[] Carray = new Color[10];
//	    
//	    Carray[0] = Color.RED;
//	    Carray[1] = Color.ORANGE;
//	    Carray[2] = Color.YELLOW;
//	    Carray[3] = Color.GREEN;
//	    Carray[4] = Color.BLUE;
//	    Carray[5] = Color.GRAY;
//	    Carray[6] = Color.PINK;
//	    Carray[7] = Color.MAGENTA;
//	    Carray[8] = Color.CYAN;
//	    Carray[9] = Color.PINK;
//	    
//	    palette2.setColor(0,Color.WHITE);
//	    for (int i = 1 ; i <= model.N ; i++){
//	    	palette2.setColor(i,Carray[i%10]);
//	    	palette2.setColor(-i,(Carray[9-i%10]));
//	    }
//	    grid1.setColors(palette2);
//	    grid2.setColors(palette2);
//	    grid3.setColors(palette2);
//	    grid4.setColors(palette2);	    
//		
//		while(true){
//			
//			
//			///////// YOU ARE TESTING ClusterV4 !!!!!!!!!!!!
//			
//			
//			dummy = new ClustersV4(model.L, "Periodic");
//			
//			foo2 = new int[model.N];
//			foo3 = new int[model.N];
//			foo4 = new int[model.N];
//			
//			for (int jj = 0 ; jj < model.N ; jj++){
//				foo[jj] = 0;
//			}
//			
//			int[] order = new int[model.N];
//			
//			for(int s = 0;s<model.N;s++) {
//				order[s] = s;
//			}
//			for(int s = 0;s<model.N-1;s++) {
//				int r = s+(int) (model.rand.nextDouble()*(model.N-s));
//				int temp = order[s];
//				order[s] = order[r];
//				order[r] = temp;
//			}
////			
////			int mnst = 25;
////			
////			order[0] = mnst+1;
////			order[1] = mnst-1;
////			order[2] = mnst+model.L;
////			order[3] = mnst-model.L;
////			order[4] = mnst;
////			order[5] = mnst+model.L+1;
////			
////			mnst=48;
////			
////			order[6] = mnst+1;
////			order[7] = mnst-1;
////			order[8] = mnst+model.L;
////			order[9] = mnst-model.L;
////			order[10] = mnst;
////			order[11] = mnst-model.L-1;
//			int count = 0;
//			display = new double[model.N];
//			for (int ii = 0 ; ii < model.N ; ii++){
//				
////				if (model.rand.nextDouble() < model.alpha){
////					dummy.addSite(order[ii]);
////					foo[order[ii]] = 1;
////				}
//				
//				
//				
//				boolean percolate = dummy.addSite(order[ii]);
//				foo[order[ii]] = 1;
//				for(int s = 0;s<model.N;s++) {
//					display[s] = dummy.getClusterSize(s);
//					foo2[s]    = dummy.getClusterNumber(s);
//					foo3[s]	   = dummy.getDist(s);
//					foo4[s]	   = dummy.getDist(s+model.N);
//				}
//				count = ii;
//				if(percolate){
//					System.out.println("System has percolated!");
//					break;
//				}
//				
//			}
//			
//			foo2[order[count]] = 4;
//			Job.animate();
//			
////			for(int s = 0;s<model.N;s++) {
////				display[s] = dummy.getClusterSize(s);
////				foo2[s]    = dummy.getClusterNumber(s);
////			}
////		    
////		    displaycounter++;
////			Job.animate();
////				
////
////			clear();	
//			    
//		}
//			  
//
//			    
//			
//	}
		
		/**
		 * 			E-N-D
		 * 	Test class Clusters(Parameters prms)
		 */
		
	
			////////////////////////////////////////////////////////////////////

			
		/**
	     *	Remind myself with R is for Square Interaction
		*/	
//		
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
		 * 	Remind myself with R is for Square Interaction
		 */
			
		
		////////////////////////////////////////////////////////////////////

		
		/**
	     *	Check Lattice Neighbors Array 
		*/	
		
		while(true){
		
			int[] temp = model.retNbs(i0);
			for (int kk = 0 ; kk < temp.length ; kk++){
				foo[temp[kk]] = 1;
			}

			foo2 = new int[model.N];

			for (int jj = 0 ; jj < model.N ; jj++){
				foo2[jj] = 0;
			}

			int[] temp2 = LN.get(i0);
			for (int kk = 0 ; kk < temp2.length ; kk++){
				foo2[temp2[kk]] = 1;
			}


			Job.animate();
		}
		
		/**
		 * 			E-N-D
	     *	Check Lattice Neighbors Array 
		*/	
		
	}
}





