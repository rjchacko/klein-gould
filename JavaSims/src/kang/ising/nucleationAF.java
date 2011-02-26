package kang.ising;

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

public class nucleationAF extends Simulation{
	
	Grid grid1 = new Grid("Grid 1");
	//Grid grid2 = new Grid("Grid 2");

	private DecimalFormat fmt = new DecimalFormat("00000");
	
	public int L1,L2,M,R,deadsites,Dseed,Bseed,Sseed;
	public double percent,biaspercent,NJ;
	public double T,H;
	public IsingStructure ISaf;
	public IsingBasic IBaf;
	
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
		grid1.registerData(IBaf.IS.L1, IBaf.IS.L2, IBaf.IS.spin);
		//grid2.setColors(ising);
		//grid2.registerData(IBaf.IS.L1, IBaf.IS.L2, IBaf.IS.biaslabel);


	}
	
	public void movie(Grid grid, int number, int copynumber)   //function to capture the grid
	{
		
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/AFnucleation/pic_"+fmt.format(copynumber)+"_"+fmt.format(number)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
	}
	
	public void clear()
	{
		grid1.clear();
		//grid2.clear();

	
	}
	
	public static void main (String[] nucleationAF){
		new Control(new nucleationAF(), "Kang Liu's AF nucleation" );
	}
	
	public void load(Control nucleationAF){
		nucleationAF.frame (grid1);
		//nucleationAF.frame (grid2);


		params.add("L1", 128);
		params.add("L2", 128);
		params.add("R", 23);
		params.add("NJ",1.0);	
		params.add("percent", new DoubleValue(0,0,1).withSlider());
		params.add("biaspercent", new DoubleValue(0,0,1).withSlider());
		params.add("deadsites");	
		//params.add("Dseed",1);
		//params.add("Bseed",1);
		//params.add("Sseed",1);
		
		params.add("Ti", 9);
		params.add("Hi", 0.95);
		params.add("Tf", 0.0965466);
		params.add("Hf", 0.95);
		
		params.addm("T", new DoubleValue(2, 0, 10).withSlider());
		params.addm("H", new DoubleValue(0, 0, 2).withSlider());
		
		params.add("MCS");
		params.add("magnetization");

	}
	
    public void run(){
		
		
		L1 = (int)params.fget("L1");
		L2 = (int)params.fget("L2");
		M = L1 * L2 ;
		
		R = (int)params.fget("R");
		NJ = params.fget("NJ");

		percent=params.fget("percent");
		biaspercent=params.fget("biaspercent");
		
		Dseed = 1;//(int)params.fget("Dseed");
		Bseed = 1;//(int)params.fget("Bseed");
		Sseed = 1;//(int)params.fget("Sseed");
		
	
	    ISaf=new IsingStructure(L1,L2,R,NJ,percent,biaspercent);
	    ISaf.Dinitialization(Dseed, Bseed, 10, 10);
	    params.set("deadsites",ISaf.deadsites);
	    ISaf.Sinitialization(0, Sseed);
	  
	    
	    T=params.fget("Ti");
	    H=params.fget("Hi");
	    
	    IBaf=new IsingBasic(ISaf,T,H,1);
	    IBaf.QuenchT=params.fget("Tf");
	    IBaf.QuenchH=params.fget("Hf");
	    
	    Job.animate();
	    
	    
	    params.set("T",IBaf.InitialT);
	    params.set("H",IBaf.InitialH);
	    for(int prestep=0; prestep<50000; prestep++)
	    {
		    IBaf.T=params.fget("T");
		    IBaf.H=params.fget("H");
			IBaf.IS.MCS(IBaf.T, IBaf.H,IBaf.flip, 1);
			Job.animate();
			params.set("MCS", 50-prestep);
			params.set("magnetization", IBaf.IS.totalspin/IBaf.IS.M);
	    }
	    
	    params.set("T",IBaf.QuenchT);
	    params.set("H",IBaf.QuenchH);
	    
	    for(int step=0; step<1000000; step++)
	    {
	    	IBaf.T=params.fget("T");
		    IBaf.H=params.fget("H");
			IBaf.IS.MCS(IBaf.T, IBaf.H,IBaf.flip, 1);
			Job.animate();
			params.set("MCS", step);
			params.set("magnetization", IBaf.IS.totalspin/IBaf.IS.M);
			if(step/10==0)
				movie(grid1,step,0);
		    
	    }
	    
	    
	    
    }

	
	
	
	
	
	
	
	
	
	
	
	
	
	
}