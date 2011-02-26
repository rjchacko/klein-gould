package kang.oldising;

import kang.ising.BasicStructure.IsingBasic;
import kang.ising.BasicStructure.IsingStructure;


import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import chris.util.PrintUtil;
import chris.util.Random;


import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DoubleValue;

public class test extends Simulation{
	
	Grid grid1 = new Grid("Grid 1");
	Grid grid2 = new Grid("Grid 2");
	Grid grid3 = new Grid("Grid 3");
	Grid grid4 = new Grid("Grid 4");
	
	
	public int L1,L2,M,R,deadsites,Dseed,Bseed,Sseed;
	public double percent,biaspercent,NJ;
	public double T,H;
	public IsingStructure IStest;
	public IsingBasic IBtest;
	
	

    
	public void animate()
	{
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);
		ising.setColor(-1, Color.WHITE);
		ising.setColor(0, Color.RED);
		ising.setColor(2, Color.BLUE);
		ising.setColor(-2, Color.GREEN);
		
		//ColorGradient heatmap = new ColorGradient();
		
		grid1.setColors(ising);
		grid1.registerData(IBtest.IS.L1, IBtest.IS.L2, IBtest.IS.spin);
		grid2.setColors(ising);
		grid2.registerData(IBtest.IS.L1, IBtest.IS.L2, IBtest.IS.biaslabel);
		grid3.setColors(ising);
		grid3.registerData(IStest.L1, IStest.L2, IStest.spin);
		grid4.setColors(ising);
		grid4.registerData(IStest.L1, IStest.L2, IStest.biaslabel);

	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
		grid3.clear();
		grid4.clear();
	
	}
	
	public static void main (String[] test){
		new Control(new test(), "Kang Liu's test program" );
	}
	
	public void load(Control test){
		test.frame (grid1);
		test.frame (grid2);
		test.frame (grid3);
		test.frame (grid4);

		params.add("L1", 200);
		params.add("L2", 200);
		params.add("R", 10);
		params.add("NJ",1.0);	
		params.add("percent", new DoubleValue(0,0,1).withSlider());
		params.add("biaspercent", new DoubleValue(0,0,1).withSlider());
		params.add("deadsites");	
		params.add("Dseed",1);
		params.add("Bseed",1);
		params.add("Sseed",1);
		
		params.addm("T", new DoubleValue(9, 0, 10).withSlider());
		params.addm("H", new DoubleValue(0, 0, 2).withSlider());
		

		params.add("Totalintenergy function");
		params.add("Totalintenergy value");		
		params.add("Totalspin function");
		params.add("Totalspin value");
		params.add("Totalenergy function");
		params.add("Totalenergy value");
		
		params.add("deadsites");
		params.add("MCS");

	}
	
    public void run(){
		
		
		L1 = (int)params.fget("L1");
		L2 = (int)params.fget("L2");
		M = L1 * L2 ;
		
		R = (int)params.fget("R");
		NJ = params.fget("NJ");

		percent=params.fget("percent");
		biaspercent=params.fget("biaspercent");
		
		Dseed = (int)params.fget("Dseed");
		Bseed = (int)params.fget("Bseed");
		Sseed = (int)params.fget("Sseed");
		
	
	    IStest=new IsingStructure(L1,L2,R,NJ,percent,biaspercent);
	    IStest.Dinitialization(Dseed, Bseed, 10, 10);
	    params.set("deadsites",IStest.deadsites);
	    IStest.Sinitialization(0, Sseed);
	  
	    
	    T=params.fget("T");
	    H=params.fget("H");
	    
	    IBtest=new IsingBasic(IStest,T,H,1);
	    Job.animate();
		for(int step=0; step<1000000;step++)
		{	
		    
		    IBtest.T=params.fget("T");
		    IBtest.H=params.fget("H");
			IBtest.MCS(IBtest.IS, IBtest.flip, 1);
			params.set("Totalintenergy function", IBtest.IS.TotalIntEnergy());
		    params.set("Totalintenergy value", IBtest.IS.totalintenergy);
		    params.set("Totalspin function", IBtest.IS.TotalSpin());
		    params.set("Totalspin value", IBtest.IS.totalspin);
		    params.set("Totalenergy function", IBtest.IS.TotalEnergy(H));
		    params.set("Totalenergy value", IBtest.IS.totalenergy);
		    params.set("deadsites",IStest.deadsites);
		    params.set("MCS",step);
		    Job.animate();
		}
		
    }
	
}