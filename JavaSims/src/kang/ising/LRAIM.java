package kang.ising;

import java.awt.Color;

import java.text.DecimalFormat;



import kang.util.PrintUtil;
import chris.util.Random;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;



import kang.ising.BasicStructure.IsingStructure;
import kang.ising.BasicStructure.BasicTools;


public class LRAIM extends Simulation{
	
	Grid grid1=new Grid("dilution configuration");     // the map to display the dilution configuration
	Grid grid2=new Grid("simulation");     // the map to display the simulation
	Grid grid3=new Grid("structure factor");
	
	
	public IsingStructure IS;
	public IsingStructure Istemp;

	public Random Erand;
	
	public BasicTools Tools;
	
	private DecimalFormat fmt = new DecimalFormat("000");
	private DecimalFormat bmt = new DecimalFormat("0000");
	
	public double sfactor[];
	
	//initialization parameters
	public int L,la,lb, Lp;
	public double M;
	public double NJ;
	public int Dseed, Bseed, Sseed;
	public double percent;
	public double biaspercent;
	public int deadsite;
	public String dynamics;

	//dynamic parameters
	public double T, H;
	
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
		grid1.registerData(L, L, IS.spin);
		grid2.setColors(ising);
		grid2.registerData(L, L, Istemp.spin);
		
		
		
		ColorGradient heatmap = new ColorGradient();
		grid3.setColors(heatmap);
		grid3.registerData(Lp, Lp, sfactor);
		
	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
		grid3.clear();

	}
	
	public static void main (String[] LRAIM){
		new Control(new LRAIM(), "Kang Liu's long range antiferromagnetic ising model" );
	}
	
	
	
	public void load(Control LRAIM)
	{

		LRAIM.frameTogether("Display", grid1 ,grid2, grid3);

		params.add("L", 200);
		params.add("Lp", 200);
		params.add("la",10);    // scale of the bias dilution region
		params.add("lb",10); 
		
		params.add("NJ", 4.0);
	    params.add("deadsites");

		params.add("percent", 0.111);
		params.add("biaspercent", 0.111);
		
		 		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));

	
		params.addm("T", 0.826);
		params.addm("H", 0.18);
		params.add("Emcs");    //MCS time for evolution
	
		    
		params.add("magnetization");

	}
	
	public void testrun(IsingStructure ising)
	{
		Random trand= new Random(1);
		for(int tstep=0; tstep<9999999; tstep++)
		{
			T=params.fget("T");
			H=params.fget("H");
			ising.MCS(T, H, trand, 1, dynamics);
			
			{
				ising.DilutionSF();
				for(int i=0; i<ising.M; i++)
				{
					sfactor[i]=ising.SFdilution.sFactor[i];
				}
			}
			
			Job.animate();
			params.set("Emcs", tstep);
			params.set("magnetization", ising.magnetization);
		}
	}
	
	public void run(){
		
		
		L = (int)params.fget("L");
		Lp = (int)params.fget("Lp");
		la = (int)params.fget("la");
		lb = (int)params.fget("lb");
		M = L * L;
		NJ = params.fget("NJ");

		percent=params.fget("percent");
		biaspercent=params.fget("biaspercent");
		dynamics= params.sget("Dynamics");
		
		Dseed = 1;
		Bseed = 1;
		Sseed = 1;

		
	    IS=new IsingStructure(L,L,0,NJ,percent,biaspercent,"square");   
	    Istemp=new IsingStructure(L,L,0,NJ,percent,biaspercent,"square");

	    sfactor= new double[Lp*Lp];
	   
	    
	    
	    Tools=new BasicTools();
	    T=params.fget("T");
	    H=params.fget("H");
	    
	    {//initialization
	    	
	    	IS.Dinitialization(Dseed, Bseed, la, lb);
	    	params.set("deadsites",IS.deadsites);
	    	IS.Sinitialization(1, Sseed);
	    	Istemp=IS.clone();
	    	
	    }
	    
	    Job.animate();
	    
	    Random rand=new Random(Sseed);
	    
	    testrun(Istemp);
	    
	  
        
	    

	}
	
}

