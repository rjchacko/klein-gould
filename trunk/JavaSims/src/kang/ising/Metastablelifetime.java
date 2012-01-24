package kang.ising;

import java.awt.Color;
import java.text.DecimalFormat;

import kang.ising.BasicStructure.IsingStructure;
import kang.ising.BasicStructure.BasicTools;


import chris.util.PrintUtil;
import chris.util.Random;

import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

public class Metastablelifetime extends Simulation
{
	Grid grid1=new Grid("grid1");     // the map to display the dilution configuration
	Grid grid2=new Grid("grid2");     // the map to display the simulation
	
	public IsingStructure IS;
	public IsingStructure Istemp;
	public BasicTools Tools;
	
	private DecimalFormat fmt = new DecimalFormat("000");
	private DecimalFormat qmt = new DecimalFormat("00000");
	
	
	//initialization parameters
	public int L;
	public double M;
	public double NJ;
	public int Dseed, Sseed;
	public double percent;
	public int deadsite;
	public int totalruns, totalcopies;
	public double threshold;
	public boolean lowerfield;
	
	
	//dynamic parameters
	public double T, H;
	public double MaxH, MinH,increment;
	
	//measurement quantities
	public double meantau;         //mean life time of the meta-stable states
	public double deltatau;         // stand deviation of tau
	public double mediantau;       //median life time  haven't put in the function to do average over different dilution realizations
	
	public double taudata[];
	public double tauSD[];
	public double taudatatemp;
	public double tauSDtemp;
	public double mediantautemp;
	
	
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
		
	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
	}
	
	public static void main (String[] Metastablelifetime){
		new Control(new Metastablelifetime(), "Kang Liu's nearest neighbor site-diluted ising model's metastable state life time mearsurement" );
	}
	
	public void load(Control Metastablelifetime)
	{
		//Metastablelifetime.frame (grid1);
		//Metastablelifetime.frame (grid2);
		Metastablelifetime.frameTogether("Display", grid1,grid2);

		params.add("L", 200);
		//params.add("M");
		params.add("NJ",-4.0);
	    params.add("deadsites");
	    //params.add("livingsites");
		params.add("percent", 0.333);
		 
		params.add("MaxH",1.0);
		params.add("MinH",0.01);
		params.add("increment",0.02);
		
		params.add("totalcopies",1);
		params.add("totalruns",50);
		
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
	
	public void Scanfield(IsingStructure ising, double T, double MaxH, double MinH, double increment, int steplimit, String dynamics)
	{
		String path = "/Users/liukang2002507/Desktop/simulation/Metalifetime/"+dynamics+"/data"+"<T="+fmt.format(T*1000)+", L= "+fmt.format(L)+">.txt";
		lowerfield=true;
		for(double field=MaxH; (field>MinH)&(lowerfield); field-=increment)
		{
			params.set("H", field);
			H=params.fget("H");
			Multicopies(ising, T, field, totalcopies, dynamics);
			if(meantau>steplimit)
				lowerfield=false;
			
			PrintUtil.printlnToFile(path , field , mediantau, meantau, deltatau);
			
		}
	}
	
	public void Multicopies(IsingStructure ising, double T, double H, int copies, String dynamics)
	{
		
		taudata= new double[copies];
		tauSD= new double[copies];
		for(int copy=0; copy<copies; copy++)
		{
			params.set("copies", copy+1);
			Dseed=copy+1;
			//params.set("Dseed", Dseed);   //set the random number seed for dilution
			ising.Dinitialization(Dseed, Dseed, 10, 10);
			ising.Sinitialization(1, Dseed);
			Istemp=ising.clone();
			params.set("deadsites",ising.deadsites);
			Job.animate();
			
			Lifetime(Istemp, T, H, totalruns, dynamics);
			taudata[copy]=taudatatemp;
			tauSD[copy]=tauSDtemp;
		}
		if (copies==1)
		{
			meantau=taudata[0];
			deltatau=tauSD[0];
			mediantau=mediantautemp;
			
		}
		else
		{
			meantau=Tools.Mean(taudata, copies);
			deltatau=Tools.SD(taudata, copies, meantau);
		}
		
	}
	
	
	public void Lifetime(IsingStructure ising, double T, double H, int runs, String dynamics)
	{
		
		String check = "/Users/liukang2002507/Desktop/simulation/Metalifetime/"+dynamics+"/check"+"<T="+fmt.format(T*1000)+", L= "+fmt.format(L)+">.txt";
		String pic="/Users/liukang2002507/Desktop/simulation/Metalifetime/"+dynamics+"/pic/<T="+fmt.format(T*1000)+", L= "+fmt.format(L)+">";
		
		int lifetimedata[]= new int[runs];
		
		for(int run=0; run< runs; run++)
		{
			params.set("runs", run+1);
			params.set("H", H);
			Random hrand= new Random(run+999);
			Random srand= new Random(run+1);
			double totalM=0;
			int totalnumber=0;
			int prelimit=3000;
			//if(meantau>100)
				//prelimit=(int)meantau*2;
			int step;
			
			/*for(int heatstep=0; heatstep<20; heatstep++)
			{
				ising.MCS(99, H, hrand, 1, dynamics);
				Job.animate();
				params.set("MCS", -99999);
			}*/
			for(int prestep=0; prestep<prelimit; prestep++)
			{
				ising.MCS(T, H, hrand, 1, dynamics);
				Job.animate();
				params.set("MCS", prestep-prelimit);
				if(prestep>(prelimit-20))
					{
					totalM+=ising.Magnetization();
					totalnumber++;
					}
				
			}
			double Ms=totalM/totalnumber;
			params.set("H", -H);
			for(step=0; ising.Magnetization()>(threshold*Ms); step++)
			{
				
				ising.MCS(T, -H, srand, 1, dynamics);
				Job.animate();
				params.set("MCS", step);
			}
			lifetimedata[run]=step;
			PrintUtil.printlnToFile(check , H , run, step);
			Tools.Picture(grid2,run+1,step, pic);
			
		}
		
		taudatatemp=Tools.Mean(lifetimedata, runs);
		
		tauSDtemp=Tools.SD(lifetimedata, runs, taudatatemp);
		mediantautemp=Tools.Median(lifetimedata, runs);
		
	
	}
	
	public void testrun(IsingStructure ising, String dynamics)
	{
		Random trand= new Random(1);
		for(int tstep=0; tstep<9999999; tstep++)
		{
			T=params.fget("T");
			H=params.fget("H");
			ising.MCS(T, H, trand, 1, dynamics);
			Job.animate();
			params.set("MCS", tstep);
			params.set("magnetization", ising.magnetization);
		}
	}
	
	public void run(){
		
		
		L = (int)params.fget("L");
		M = L * L;
		NJ = params.fget("NJ");

		percent=params.fget("percent");
		String dynamics= params.sget("Dynamics");
		//Dseed = (int)params.fget("Dseed");
		//Sseed = (int)params.fget("Sseed");
		MaxH = params.fget("MaxH");
		MinH = params.fget("MinH");
		increment =params.fget("increment");
		totalcopies=(int)params.fget("totalcopies");
		totalruns=(int)params.fget("totalruns");
		threshold=0.9;
		
	    IS=new IsingStructure(L,L,0,NJ,percent,percent,"square");   
	    Istemp=new IsingStructure(L,L,0,NJ,percent,percent,"square");
	    Tools=new BasicTools();
	    T=params.fget("T");
	    
	    /*{ 
	    	IS.Dinitialization(Dseed, Dseed, 10, 10);
	    	params.set("deadsites",IS.deadsites);
	    	IS.Sinitialization(1, Sseed);
	    	Istemp=IS.clone();
	    }*/
	    
	    Job.animate();
	    
	    //testrun(Istemp, dynamics);
	    //meantau=20000;
	    //Scanfield(IS, T, MaxH, MinH, increment, 20000, dynamics);
	    
	    
	    Scanfield(IS, T, 0.60, 0.30, 0.05, 200000, dynamics);
	    Scanfield(IS, T, 0.30, 0.20, 0.02, 200000, dynamics);
	    Scanfield(IS, T, 0.20, 0.10, 0.01, 200000, dynamics);
	    Scanfield(IS, T, 0.10, 0.01, 0.005, 200000, dynamics);

      
        
	    
	    
	    
	    
	    
	    

	}
}
