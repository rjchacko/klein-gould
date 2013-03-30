package kang.ising;

import java.awt.Color;
import java.text.DecimalFormat;

import chris.util.PrintUtil;
import chris.util.Random;

import kang.ising.BasicStructure.FCIsing;
import kang.ising.BasicStructure.IsingStructure;
import kang.ising.BasicStructure.BasicTools;

import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;


public class DilutionCorrelation extends Simulation
{
	Grid Gspin=new Grid("spin");
	Grid GE=new Grid("Energy Fluctuation");
	Grid GM=new Grid("Spin Fluctuation");
	Grid Gdilution=new Grid("dilution");
	
	public int L1,L2,M,R,deadsites,Dseed,Bseed,Sseed;
	public double percent,biaspercent,NJ;
	public double T,H;
	public IsingStructure IS;
	public IsingStructure Istemp;
	public BasicTools Tools;
	
	public int progress;
    public int usefulruns=0;
    
    
    public double tempE[];
    public double tempM[];
    public double totalE[];
    public double totalM[];
    public double totalE2[];
    public double totalM2[];
    
    
    public double dilutionmap[];
    public double deltaenergy[];
    public double deltaspin[];
    
	private DecimalFormat fmt = new DecimalFormat("000");
	public double startH=0;
	
	
	
	public void animate()
	{
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);      //up spin
		ising.setColor(-1, Color.WHITE);     //down spin
		ising.setColor(0, Color.RED);        //normal dilution
		ising.setColor(2, Color.BLUE);       //clusters
		ising.setColor(-2, Color.GREEN);     //
		ising.setColor(3, Color.darkGray);    // the centers of the clusters
		
		Gspin.setColors(ising);
		Gspin.registerData(Istemp.L1, Istemp.L2, Istemp.spin);
		
		
		ColorGradient heatmap = new ColorGradient();
		
		Gdilution.setColors(heatmap);
		Gdilution.registerData(Istemp.L1, Istemp.L2, dilutionmap);
	
		GE.setColors(heatmap);
		GE.registerData(Istemp.L1, Istemp.L2, deltaenergy);
		
		GM.setColors(heatmap);
		GM.registerData(Istemp.L1, Istemp.L2, deltaspin);
		
		
		params.set("magnetization", Istemp.magnetization);
		params.set("intenergy", Istemp.totalintenergy);
	}

	public void clear()
	{
		Gspin.clear();
		Gdilution.clear();
		GE.clear();
		GM.clear();
	}
	
	public static void main (String[] DilutionCorrelation){
		new Control(new DilutionCorrelation(), "Kang Liu's dilution correlation" );
	}

	public void load(Control DilutionCorrelation){
		
		DilutionCorrelation.frameTogether ("Display", Gspin, Gdilution, GE, GM);
		

		params.add("L1", 100);
		params.add("L2", 100);
		params.add("R", 5);
		params.add("NJ",-4.0);	
		params.add("percent", 0.10);
		params.add("biaspercent", 0.10);
		params.add("deadsites");	
		params.add("Dseed",2);
		params.add("Bseed",2);
		params.add("Sseed",1);
		
		params.addm("T", 1.600);
		params.addm("H", 0.0);
		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));
		
		params.add("MCS");
		params.add("copies");
		params.add("magnetization");
		params.add("intenergy");
	}
	
	public void scanfield(IsingStructure ising, double t, double maxh, double minh, double dh, int presteplimit, int steplimit, int seed, String dynamics, Boolean print)
	{
		String SaveScan= "/Users/liukang2002507/Desktop/simulation/DilutionCorrelation/"+dynamics+"/Correlation data <q=0."+fmt.format(ising.percent*1000)+", biasq=0."+fmt.format(ising.biaspercent*1000)+", L="+fmt.format(ising.L1)+", R="+fmt.format(ising.R)+", seed="+fmt.format(seed)+">.txt";
	    for(double h=maxh; h>minh; h-=dh)
	    {
	    	double correlation[]=new double[2];
	    	params.set("H", h);
	    	correlation=spinodalrun(ising, t, h, presteplimit, steplimit, seed, dynamics, false);
	    	
	    	PrintUtil.printlnToFile(SaveScan, h, correlation[0], correlation[1]);
	    }
	
	
	}
	
	public double[] spinodalrun(IsingStructure ising, double t, double field, int presteplimit, int steplimit, int seed, String dynamics, Boolean print)
	{
		String path= "/Users/liukang2002507/Desktop/simulation/DilutionCorrelation/"+dynamics+"/SP data <q=0."+fmt.format(ising.percent*1000)+", biasq=0."+fmt.format(ising.biaspercent*1000)+", L="+fmt.format(ising.L1)+", R="+fmt.format(ising.R)+", H="+fmt.format(field*1000)+", seed="+fmt.format(seed)+">.txt";
		double[] correlation=new double[2];
		
		
		Istemp=ising.clone();
		Random cflip=new Random(seed);
	
        tempE=new double[Istemp.M];
        tempM=new double[Istemp.M];
        totalE=new double[Istemp.M];
        totalM=new double[Istemp.M];
        totalE2=new double[Istemp.M];
        totalM2=new double[Istemp.M];
        
		for(int ini=0;ini<Istemp.M; ini++)
		{
			tempE[ini]=0;
			tempM[ini]=0;
			totalE[ini]=0;
			totalE2[ini]=0;
			totalM[ini]=0;
			totalM2[ini]=0;
		}
        
        
		
		
		for(int heat=0; heat<5; heat++)
		{
			Istemp.MCS(9, field, cflip, 1, dynamics);
			Job.animate();
			params.set("MCS", -9999);

		}
	
		for(int prestep=0; prestep< presteplimit; prestep++)
		{
			Istemp.MCS(t, field, cflip, 1, dynamics);
			Job.animate();
			params.set("MCS", prestep-presteplimit);
		}
		params.set("H",-field);
		
		int count=0;
		for(int step=0; step<steplimit; step++)
		{
			
			Istemp.MCS(t, -field, cflip, 1, dynamics);
			Job.animate();
			params.set("MCS", step);
			params.set("magnetization", Istemp.Magnetization()/Istemp.M);
			
			if((step>500)&&(step<=1500))
			{
				count++;
				tempE=Istemp.LocalEnergy(-field);
				for(int jj=0; jj<Istemp.M; jj++)
				{
					
						totalE[jj]+=tempE[jj];
						totalE2[jj]+=tempE[jj]*tempE[jj];
						totalM[jj]+=Istemp.spin[jj];
						totalM2[jj]+=Istemp.spin[jj]*Istemp.spin[jj];
					
					
				}
			}
			
		}
		
		for(int jj=0; jj<Istemp.M; jj++)
		{
			if(Istemp.spin[jj]!=0)
			{
				deltaenergy[jj]=totalE2[jj]/count-(totalE[jj]/count)*(totalE[jj]/count);
				deltaspin[jj]=totalM2[jj]/count-(totalM[jj]/count)*(totalM[jj]/count);
				dilutionmap[jj]=Istemp.dilutionratioSquare(Istemp.R, jj);
			}
			
		}
		Job.animate();
		
		if(print)
			output(path,Istemp.M);
		
		correlation[0]=Tools.DCorrelation(dilutionmap, deltaspin, Istemp.M);
		correlation[1]=Tools.Correlation(dilutionmap, deltaspin, Istemp.M);

		return correlation;
			
		
	}
	
	public void output(String path, int total)
	{
		for(int i=0; i<total; i++)
		{
			PrintUtil.printlnToFile(path, i, deltaenergy[i], deltaspin[i], dilutionmap[i]);
		}
	}
	
	public void run(){
		
		
		L1 = (int)params.fget("L1");
		L2 = (int)params.fget("L2");
		M = L1 * L2;
		
		R = (int)params.fget("R");
		NJ = params.fget("NJ");

		percent=params.fget("percent");
		biaspercent=params.fget("biaspercent");
		String dynamics= params.sget("Dynamics");
		
		Dseed = (int)params.fget("Dseed");
		Bseed = (int)params.fget("Bseed");
		Sseed = (int)params.fget("Sseed");
		
	    IS=new IsingStructure(L1,L2,R,NJ,percent,biaspercent,"square");
	    Istemp=new IsingStructure(L1,L2,R,NJ,percent,biaspercent,"square");
	    Tools=new BasicTools();
	    
	    dilutionmap=new double[M];
	    deltaenergy=new double[M];
	    deltaspin=new double[M];
	    
	    IS.Dinitialization(Dseed, Bseed, 10, 10);
	    params.set("deadsites",IS.deadsites);
	    IS.Sinitialization(0, Sseed);
	    
	    Job.animate();

        
	    double[] ctemp=new double[2];
	    
	   {
	    	T=params.fget("T");
	    	
	    	//ctemp=spinodalrun(IS, T, 0.9*(1-percent), 200, 2500, 1, dynamics, true);
	    	
	    	scanfield(IS, T, 1.05*(1-percent), 0.10*(1-percent), 0.05*(1-percent), 200, 2500, 1, dynamics, false);
	    
	   }
	    
	    
	   

	    
	    

	}
	



}