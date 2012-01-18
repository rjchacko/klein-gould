package kang.ising;

import kang.ising.BasicStructure.BasicTools;
import kang.ising.BasicStructure.SPpercolation;
import kang.ising.BasicStructure.IsingStructure;


import java.awt.Color;
import java.text.DecimalFormat;

import chris.util.PrintUtil;
import chris.util.Random;

import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

public class spinodaldroplet extends Simulation
{
	Grid grid1=new Grid("grid1");
	Grid grid2=new Grid("grid2");
	public int L,M,R,deadsites,Dseed,Bseed,Sseed;
	public double percent,biaspercent,NJ;
	public double T,H;
	public IsingStructure IS;
	public IsingStructure Istemp;
	public BasicTools Tools;
	public SPpercolation SPP;
	
	private DecimalFormat fmt = new DecimalFormat("000");

	
	public void animate()
	{
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);      //up spin
		ising.setColor(-1, Color.WHITE);     //down spin
		ising.setColor(0, Color.RED);        //normal dilution
		ising.setColor(2, Color.BLUE);       //clusters
		ising.setColor(-2, Color.GREEN);     //
		ising.setColor(3, Color.darkGray);    // the centers of the clusters
		
		
		grid1.setColors(ising);
		grid1.registerData(IS.L1, IS.L2, IS.spin);
		grid2.setColors(ising);
		grid2.registerData(Istemp.L1, Istemp.L2, Istemp.spin);
		
		params.set("magnetization", Istemp.magnetization);
	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
	}
	
	public static void main (String[] spinodaldroplet){
		new Control(new spinodaldroplet(), "Kang Liu's percolation mapping for the spinodal nucleating droplet" );
	}
	
	public void load(Control spinodaldroplet)
	{
		//spinodaldroplet.frame (grid1);
		//spinodaldroplet.frame (grid2);
		spinodaldroplet.frameTogether("Display", grid1,grid2);

		params.add("L", 200);
		params.add("NJ",-4.0);
	    params.add("deadsites");

		params.add("percent", 0.333);
		
		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));

		//params.add("Dseed",1);    //seed for dilution configuration
		//params.add("Sseed",1);    //seed for spin flip
		
		
		
		params.addm("T", 0.402);
		params.addm("H", 0.0);
		
		
		
		params.add("MCS");
		params.add("runs");
		params.add("copies");
		params.add("magnetization");
		

	}
	
	public void run(){
		
		
		L = (int)params.fget("L");
		M = L * L;
		
		
		NJ = params.fget("NJ");

		percent=params.fget("percent");
		String dynamics= params.sget("Dynamics");
		
		
		//Dseed = (int)params.fget("Dseed");
		//Sseed = (int)params.fget("Sseed");


		
		
	    IS=new IsingStructure(L,L,0,NJ,percent,percent,"square");   
	    Istemp=new IsingStructure(L,L,0,NJ,percent,percent,"square");
	    Tools=new BasicTools();
	    T=params.fget("T");
	    
	    /*{ 
	    	IS.Dinitialization(Dseed, Dseed, 10, 10);
	    	params.set("deadsites",IS.deadsites);
	    	IS.Sinitialization(0, Sseed);
	    	Istemp=IS.clone();
	    }*/
	    
	    Job.animate();
	    
	
      
        
	    
	    
	    
	    
	    
	    

	}
	
	
}