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
	public int L,la,lb,R, Lp;
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

		params.add("L", 256);
		params.add("Lp", 256);
		params.add("la",10);    // scale of the bias dilution region
		params.add("lb",10); 
		params.add("R", 92);
		
		params.add("NJ", 4.0);
	    params.add("deadsites");

		params.add("percent", 0.00);
		params.add("biaspercent", 0.00);
		
		 		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));

	
		params.addm("T", 1.0);
		params.addm("H", 0.2);
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
				ising.SpinSF();
				for(int i=0; i<ising.M; i++)
				{
					sfactor[i]=ising.SFdown.sFactor[i];
					
				}
				sfactor[Lp/2*Lp+Lp/2]=0;
				
			}
			
			Job.animate();
			params.set("Emcs", tstep);
			params.set("magnetization", ising.magnetization);
		}
	}
	
	public void temperatureScan(IsingStructure ising, double Tmax, double Tmin, double dT, int steplimit, int seed)
	{
		for(double t=Tmax; t>=Tmin; t-=dT)
		{
			singlerun(ising, t, H, steplimit, seed);
		}
	}
	
	public void evolution(IsingStructure ising, double T, double H,int steplimit, int seed)
	{
		String evolutionrun="evolution <L="+fmt.format(L)+", R="+fmt.format(ising.R)+", la= "+fmt.format(la)+", lb= "+fmt.format(lb)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String evolutionpath = "/Users/liukang2002507/Desktop/simulation/LRAIM/"+dynamics+"/"+evolutionrun+"[T="+fmt.format(T*100)+", H="+fmt.format(H*100)+"]"+".txt";
		
		Random rand= new Random(seed);
		for(int prestep=0; prestep<50; prestep++)
		{
			params.set("T", 9);
			params.set("H", H);
			ising.MCS(9, H, rand, 1, dynamics);
			Job.animate();
			params.set("Emcs", prestep-50);
			params.set("magnetization", ising.magnetization);
		}
		
		for(int step=0; step<steplimit; step++)
		{
			params.set("T", T);
			params.set("H", H);
			ising.MCS(T, H, rand, 1, dynamics);
			
			{
				ising.SpinSF();
				for(int i=0; i<ising.M; i++)
				{
					sfactor[i]=ising.SFdown.sFactor[i];
					
				}
				sfactor[Lp/2*Lp+Lp/2]=0;
				
			}
			
			Job.animate();
			params.set("Emcs", step);
			params.set("magnetization", ising.magnetization);
			int bestpeak1=ising.SFup.findBestSquareInt(ising.R);
			int bestpeak2=ising.SFdown.findBestSquareInt(ising.R);
			
			
			PrintUtil.printlnToFile(evolutionpath, step, bestpeak1, ising.SFup.squareSF[bestpeak1], bestpeak2, ising.SFdown.squareSF[bestpeak2]);
			
			
			
		}
		
		
		
		
	}
	
	
	public void singlerun(IsingStructure ising, double T, double H,int steplimit, int seed)
	{
		String singlerun="<L="+fmt.format(L)+", Lp="+fmt.format(Lp)+", la= "+fmt.format(la)+", lb= "+fmt.format(lb)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String singlepath = "/Users/liukang2002507/Desktop/simulation/LRAIM/"+dynamics+"/"+singlerun;
		
		Random rand= new Random(seed);
		for(int prestep=0; prestep<50; prestep++)
		{
			params.set("T", 9);
			params.set("H", H);
			ising.MCS(9, 0, rand, 1, dynamics);
			Job.animate();
			params.set("Emcs", prestep-50);
			params.set("magnetization", ising.magnetization);
		}
		
		for(int step=0; step<steplimit; step++)
		{
			params.set("T", T);
			params.set("H", H);
			ising.MCS(T, H, rand, 1, dynamics);
			
			{
				ising.SpinSF();
				for(int i=0; i<ising.M; i++)
				{
					sfactor[i]=ising.SFdown.sFactor[i];
					
				}
				sfactor[Lp/2*Lp+Lp/2]=0;
				
			}
			
			Job.animate();
			params.set("Emcs", step);
			params.set("magnetization", ising.magnetization);
			
		}
		
		outputSquareSF(ising, singlepath, T, H);
		
		
	}
	
	public void outputSquareSF(IsingStructure ising, String path, double T, double H)
	{
		String newpath= path+"[T="+fmt.format(T*100)+", H="+fmt.format(H*100)+"]"+".txt";
		double temp1=ising.SFup.squareSF[1];
		double temp2=ising.SFdown.squareSF[1];
		int peak1=1;
		int peak2=1;
		for(int r=0; r<Lp/2; r++)
		{
			PrintUtil.printlnToFile(newpath , r , ising.SFup.squareSF[r], ising.SFdown.squareSF[r]);
			if(r!=0)
				{
				if(ising.SFup.squareSF[r]>temp1) 
					{
					peak1=r;
					temp1=ising.SFup.squareSF[r];
					}
				if(ising.SFdown.squareSF[r]>temp2) 
					{
					peak2=r;
					temp2=ising.SFdown.squareSF[r];
					}
				}
					
		}
		PrintUtil.printlnToFile(newpath , "peak1 at r=  " , peak1,  ising.SFup.squareSF[peak1]);
		PrintUtil.printlnToFile(newpath , "peak2 at r=  " , peak2,  ising.SFdown.squareSF[peak2]);
		
		int bestpeak1=ising.SFup.findBestSquareInt(ising.R);
		int bestpeak2=ising.SFdown.findBestSquareInt(ising.R);
		
		PrintUtil.printlnToFile(newpath , "theorypeak1 at r=  " , bestpeak1,  ising.SFup.squareSF[bestpeak1]);
		PrintUtil.printlnToFile(newpath , "theorypeak2 at r=  " , bestpeak2,  ising.SFdown.squareSF[bestpeak2]);
		
		
		
		//now print to overall entry
		PrintUtil.printlnToFile(path+".txt" , T,  H, ising.SFup.squareSF[peak1], ising.SFdown.squareSF[peak2], peak1, peak2);
		
		
	}
	
	
	public void run(){
		
		
		L = (int)params.fget("L");
		Lp = (int)params.fget("Lp");
		la = (int)params.fget("la");
		lb = (int)params.fget("lb");
		R =(int)params.fget("R");
		M = L * L;
		NJ = params.fget("NJ");

		percent=params.fget("percent");
		biaspercent=params.fget("biaspercent");
		dynamics= params.sget("Dynamics");
		
		Dseed = 1;
		Bseed = 1;
		Sseed = 1;

		
	    IS=new IsingStructure(L,L,R,NJ,percent,biaspercent,"square");   
	    Istemp=new IsingStructure(L,L,R,NJ,percent,biaspercent,"square");

	    sfactor= new double[Lp*Lp];
	   
	    
	    
	    Tools=new BasicTools();
	    //T=params.fget("T");
	    H=params.fget("H");
	    
	    {//initialization
	    	
	    	IS.Dinitialization(Dseed, Bseed, la, lb);
	    	params.set("deadsites",IS.deadsites);
	    	IS.Sinitialization(1, Sseed);
	    	Istemp=IS.clone();
	    	
	    }
	    
	    Job.animate();
	   
	    
	    //testrun(Istemp);
	    
	    //singlerun(Istemp, 0.7, 0, 200, 1);
	  
	    //temperatureScan(Istemp, 1.20, 0.10, 0.02, 1000, 1);
        
	    evolution(Istemp, 0.85, 0, 100, 1);

	}
	
}

