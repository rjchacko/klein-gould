package kang.ising;

import java.awt.Color;
import java.text.DecimalFormat;

import kang.ising.BasicStructure.BasicTools;
import kang.ising.BasicStructure.IsingStructure;
import kang.util.PrintUtil;
import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import chris.util.Random;


public class ClusterScaling extends Simulation{
	
	Grid grid1=new Grid("dilution configuration");     // the map to display the dilution configuration
	Grid grid2=new Grid("simulation");     // the map to display the simulation
	Grid grid3=new Grid("largest cluster");
	
	
	public IsingStructure IS;
	public IsingStructure Istemp;

	public Random Erand;
	
	public BasicTools Tools;
	
	private DecimalFormat fmt = new DecimalFormat("000");
	private DecimalFormat bmt = new DecimalFormat("0000");
	
	
	
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
	
	//cluster parameters
	public int[] clustersize;
	
	
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
		
		grid3.setColors(ising);
		grid3.registerData(L, L, Istemp.largestcluster);
		
	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
		grid3.clear();

	}
	
	public static void main (String[] ClusterScaling){
		new Control(new ClusterScaling(), "Kang Liu's cluster scaling" );
	}
	
	public void printdata(String path, int[] data)
	{
		int j=0;
		while(data[j]>=0)
		{
			PrintUtil.printlnToFile(path, j+1, data[j]);
			j++;
		}
	}
	
	public void singlerunTc(IsingStructure ising, double T, int steplimit, boolean keeplargest,int seed)
	{
		String singlerun="<L="+fmt.format(L)+", R="+fmt.format(R)+", T="+fmt.format(T*100)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String singlepath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling/"+dynamics+"/"+singlerun+".txt";
		
		Random rand= new Random(seed);
		for(int prestep=0; prestep<50; prestep++)
		{
			params.set("T", 999);
			params.set("H", 0);
			ising.MCS(999, 0, rand, 1, dynamics);
			Job.animate();
			params.set("Emcs", prestep-50);
			params.set("magnetization", ising.magnetization);
		}
		
		for(int step=0; step<steplimit; step++)
		{
			params.set("T", T);
			params.set("H", H);
			ising.MCS(T, H, rand, 1, dynamics);
			
			Job.animate();
			params.set("Emcs", step);
			params.set("magnetization", ising.magnetization);
			
		}
		double pb=1-Math.exp(2*ising.J/T);
		int direction=1;
		if(ising.magnetization<0)
			direction=-1;
			
		clustersize=ising.Clustergrowth(ising.spin, direction, pb, seed, seed, keeplargest);
		printdata(singlepath, clustersize);
		
	}
	
	
	
	public void load(Control ClusterScaling)
	{

		ClusterScaling.frameTogether("Display", grid1 ,grid2, grid3);

		params.add("L", 400);
		
		params.add("la",10);    // scale of the bias dilution region
		params.add("lb",10); 
		params.add("R", 20);
		
		params.add("NJ", -4.0);
	    params.add("deadsites");

		params.add("percent", 0.00);
		params.add("biaspercent", 0.00);
		
		 		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));

	
		params.addm("T", 10.00);
		params.addm("H", 0.0);
		params.add("Emcs");    //MCS time for evolution
	
		    
		params.add("magnetization");

	}
	
	
	public void run(){
		
		
		L = (int)params.fget("L");

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

	    
	    Tools=new BasicTools();
	    T=params.fget("T");
	    H=params.fget("H");
	    
	    {//initialization
	    	
	    	IS.Dinitialization(Dseed, Bseed, la, lb);
	    	params.set("deadsites",IS.deadsites);
	    	IS.Sinitialization(1, Sseed);
	    	Istemp=IS.clone();
	    	Istemp.largestcluster=new int[IS.M];
	    	
	    }
	    
	    Job.animate();
	   
	    singlerunTc(Istemp, T, 500,true, 1);
	    Job.animate();

	}
	
}
