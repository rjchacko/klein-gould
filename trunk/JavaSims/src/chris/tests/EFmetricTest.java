package chris.tests;


import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;

import scikit.dataset.Histogram;
import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.ofc.old.NfailDamage2D;
import chris.util.PrintUtil;

public class EFmetricTest extends Simulation {

	Grid grid1 = new Grid ("Stress Lattice");
	Grid grid2 = new Grid ("Failed Sites");
	Plot plot1 = new Plot("Histogram of Number of Showers");
	Plot plot2 = new Plot("Radius of Gyration as a Function of Time");
	NfailDamage2D model;
	Histogram histRgyr;	
	double ScMax, rgyr, Omega, OmegaOld;
	DecimalFormat fmt = new DecimalFormat("0000000");
	DecimalFormat fmts = new DecimalFormat("0000");
	
	ColorPalette palette1;
	ColorGradient smooth;
	
	int Echk = 0;
	

	public static void main(String[] args) {
		new Control(new EFmetricTest(), "OFC Model");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/CurrentSemester/Research/Test/"));
		params.add("Random Seed",0);
		params.addm("Auto Scale", new ChoiceValue("Yes", "No"));
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
		//c.frameTogether("Observables (OFC Model with Damage)",plot1,plot2);
		
	}
	
	public void animate() {
		
//		if (params.sget("Auto Scale").equals("Yes")) {
//			plot1.setAutoScale(true);
//			plot2.setAutoScale(true);
//		}
//		else{
//			plot1.setAutoScale(false);
//			plot2.setAutoScale(false);
//		}
//						
		params.set("Number of Resets",model.time);
		params.set("Number of Showers",model.showernumber);
//		
//		int[] foo = new int[model.N];
//
//			
//		for (int i=0 ; i<model.N ; i++){
//			smooth.getColor(model.stress[i],-2,ScMax);
//			foo[i]=model.alive[i];
//		}
//			
//		grid1.setColors(smooth);
//		grid1.registerData(model.L,model.L,model.stress);
//		grid2.registerData(model.L, model.L, foo);
			


	}

	public void clear() {
		plot1.clear();
		plot2.clear();
		grid1.clear();
		grid2.clear();
	}

	public void run() {
		
//		model = new NfailDamage2D(params);
//				
//		OmegaOld=0;
//		model.time=1;
//		model.showernumber=0;
//		String record;
		
		// Read in Stress Matrix etc.
		// 4 columns tab delimited

//		int ii = 0;
//		
//         double[] StressSoFar = new double[model.N];
//         
//         try{
//             File f = new File(filein); 
//             FileInputStream fis = new FileInputStream(f); 
//             BufferedInputStream bis = new BufferedInputStream(fis); 
//             DataInputStream dis = new DataInputStream(bis);  
//
//             
//             // skip first line?
//             
//             record=dis.readLine();
//             
//             
//             while ( (record=dis.readLine()) != null ) { 
//            	 
//            	 int posDelimiter = record.indexOf('\t');            	 
//            	 //System.out.println(record.substring(0, posDelimiter));
//         		 model.stress[ii]=Double.parseDouble(record.substring(0, posDelimiter));
//            	 String temp1 = record.substring(posDelimiter + 1);
//         		 posDelimiter = temp1.indexOf('\t');
//         		 //System.out.println(temp1.substring(0, posDelimiter));
//         		 model.alive[ii]=(int)(Double.parseDouble(temp1.substring(0, posDelimiter)));
//         		 String temp2 = temp1.substring(posDelimiter + 1);
//         		 posDelimiter = temp2.indexOf('\t');
//         		 //System.out.println(temp2.substring(0, posDelimiter)); 
//         		 model.alive[ii+model.N]=(int)(Double.parseDouble(temp2.substring(0, posDelimiter)));	 
//         		 String temp3 = temp2.substring(posDelimiter + 1);
//         		 //System.out.println(temp3); 
//         		 StressSoFar[ii++]=Double.parseDouble(temp3);
//            	 
//            	 //Job.animate();
//            	 
//             }
//        	 
//         }
//         catch (IOException e) { 
//             System.out.println("ERROR!" + e.getMessage()); 
//
//         }
//         
//         
//         
//  		String ffout = model.outdir + File.separator+"DoesThisWork.txt";
//          
//          for (int kk = 0; kk < model.N ; kk++){
//         	 PrintUtil.printlnToFile(ffout,model.stress[kk],(double)model.alive[kk],(double)model.alive[kk+model.N],StressSoFar[kk]);
//          }
//          
//          System.out.println("Done!");
 		
         
         // Put testing code here (set time by hand?)
         
		
		model = new NfailDamage2D(params);
		
		OmegaOld=0;
		model.time=1;
		model.showernumber=0;
		
		String record;
		int fn=0;
         
         for (int nfiles = 0 ; nfiles < 26 ; nfiles++){
        	 
        	 int ii=0;
        	 String datafile = model.outdir+File.separator+"DamageCHK"+fmts.format(fn)+".txt";
        	 double[] StressSoFar = new double[model.N];
        	 
        	 try{
		           File f = new File(datafile); 
		           FileInputStream fis = new FileInputStream(f); 
		           BufferedInputStream bis = new BufferedInputStream(fis); 
		           BufferedReader bir = new BufferedReader(new InputStreamReader(bis));		
		           
		           // skip first line
		           record=bir.readLine();
		           
		           
		           while ( (record=bir.readLine()) != null ) { 
		          	 
		          	 int posDelimiter = record.indexOf('\t');            	 
		          	 //System.out.println(record.substring(0, posDelimiter));
		       		 model.stress[ii]=Double.parseDouble(record.substring(0, posDelimiter));
		          	 String temp1 = record.substring(posDelimiter + 1);
		       		 posDelimiter = temp1.indexOf('\t');
		       		 //System.out.println(temp1.substring(0, posDelimiter));
		       		 model.alive[ii]=(int)(Double.parseDouble(temp1.substring(0, posDelimiter)));
		       		 String temp2 = temp1.substring(posDelimiter + 1);
		       		 posDelimiter = temp2.indexOf('\t');
		       		 //System.out.println(temp2.substring(0, posDelimiter)); 
		       		 model.alive[ii+model.N]=(int)(Double.parseDouble(temp2.substring(0, posDelimiter)));	 
		       		 String temp3 = temp2.substring(posDelimiter + 1);
		       		 //System.out.println(temp3); 
		       		 StressSoFar[ii++]=Double.parseDouble(temp3);
		          	 
		          	 //Job.animate();
		          	 
		           }
 
        	 }
        	 catch (IOException e) { 
        		 System.out.println("ERROR!" + e.getMessage()); 
        	 }
        	 
        	 for (int jj = 0 ; jj < model.N ; jj++){
        		 model.SsoFar[jj]=StressSoFar[jj]-model.stress[jj];
        	 }
        	 
        	 double ckOmega = model.EFmetric();
        	 
         	 PrintUtil.printlnToFile(model.outdir+File.separator+"Recalc.txt",(double) fn++, ckOmega);

         	 System.out.println("Done with file number " + fn);
         	 
         }
         
      	 System.out.println("Finished");
         
	}
}

         
         