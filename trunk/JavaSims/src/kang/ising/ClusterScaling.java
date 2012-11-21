package kang.ising;

import java.awt.Color;
import java.text.DecimalFormat;

import kang.ising.BasicStructure.BasicTools;
import kang.ising.BasicStructure.IsingStructure;
import kang.ising.BasicStructure.IsingStructure3D;
import kang.util.PrintUtil;
import kang.ising.BasicStructure.BasicTools;

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
	
	public int adddata(String path, int[] data, int start)
	{
		int j=0;
		int end=0;
		while(data[j]>=0)
		{
			PrintUtil.printlnToFile(path, j+start, data[j]);
			j++;
		}
		
		end=j+start;
		return end;
	}
	
	public double meansize(int[] data)
	{
		double totalsize=0;
		double msize=0;
		int j=0;
		while(data[j]>=0)
		{
			totalsize+=(data[j]*data[j]);
			j++;
		}
		if(j>0)
		{
			msize=totalsize;
		}
		
		return msize;
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
			params.set("H", 0);
			ising.MCS(T, 0, rand, 1, dynamics);
			
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
	
	public void multiplerunsTc(IsingStructure3D ising, double T, int steplimit, boolean keeplargest,int seed, int copies)// average over multiple spin configurations, positve represent all up spin clusters
	{
		String multirun="multi <L="+fmt.format(L)+", R="+fmt.format(R)+", T="+fmt.format(T*100)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String positivepath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling/"+dynamics+"/"+multirun+"positive.txt";
		String negativepath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling/"+dynamics+"/"+multirun+"negative.txt";
		
		int spincopies[][]=new int [copies][ising.M];
		
		
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
			params.set("H", 0);
			ising.MCS(T, 0, rand, 1, dynamics);
			
			Job.animate();
			params.set("Emcs", step);
			params.set("magnetization", ising.magnetization);
			
		}
		
		double pb=1-Math.exp(2*ising.J/T);
		
		for(int astep=0; astep<(10*copies); astep++)
		{
			params.set("T", T);
			params.set("H", 0);
			ising.MCS(T, 0, rand, 1, dynamics);
			
			Job.animate();
			params.set("Emcs", steplimit+astep);
			params.set("magnetization", ising.magnetization);
			
			if(astep%10==0)
			{
				for(int jj=0; jj<ising.M; jj++)
				{
					spincopies[astep/10][jj]=ising.spin[jj];
				}
			}
		}
		
		int Pstart=0;
		int Nstart=0;
        for(int cc=0; cc<copies; cc++)
        {
        	
        	clustersize=ising.Clustergrowth(spincopies[cc], 1, pb, seed, seed, keeplargest);
    		int center=ising.ClusterInfo(ising.largestcluster)[1];
    		
    		Job.animate();
    		Pstart=adddata(positivepath, clustersize, Pstart);
    		
    		clustersize=ising.Clustergrowth(spincopies[cc], -1, pb, seed, seed, keeplargest);
    	    center=ising.ClusterInfo(ising.largestcluster)[1];
    		
    		Job.animate();
    		Nstart=adddata(negativepath, clustersize, Nstart);
        	
        }
	
			
		
	}
	
	
	public void singlerunHs(IsingStructure ising, double T, double H, int steplimit, boolean keeplargest,int seed)
	{
		String singlerun="Hs <L="+fmt.format(L)+", R="+fmt.format(R)+", T="+fmt.format(T*100)+", H="+bmt.format(H*1000)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String singlepath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling/"+dynamics+"/"+singlerun+".txt";
		
		Random rand= new Random(seed);
        double totalm=0;
        int mnumber=0;
		
		for(int step=0; step<steplimit; step++)
		{
			params.set("T", T);
			params.set("H", -H);
			ising.MCS(T, -H, rand, 1, dynamics);
			
			Job.animate();
			params.set("Emcs", step);
			params.set("magnetization", ising.magnetization);
			if(step>steplimit-100)
			{
				totalm+=ising.magnetization;
				mnumber++;
			}
			
			
		}
		
		double mag=totalm/mnumber;
	    double pb=1-Math.exp((1+mag/(1-ising.percent))*ising.J/T);
		int direction=-1;
			
		clustersize=ising.Clustergrowth(ising.spin, direction, pb, seed, seed, keeplargest);
		printdata(singlepath, clustersize);
		
	}
	
	public void chi(IsingStructure ising, double T, double Hmin, double Hmax, double dH, int steplimit, int seed)  // measure the mean size of clusters as a function of the field
	{
		String chirun="chi data <L="+fmt.format(L)+", R="+fmt.format(R)+", T="+fmt.format(T*100)+"p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String chipath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling/"+dynamics+"/"+chirun+".txt";
		
		double chitemp=0;
		double totalms=0;
		
		for(double h=Hmin; h<Hmax; h+=dH)
		{
		  
			totalms=meanclustersize(ising, T, h, steplimit, seed);
			chitemp=2*totalms/(1-percent-ising.magnetization);
			PrintUtil.printlnToFile(chipath, h, chitemp);
			
		}
		
	}
	
	public void approachHs(IsingStructure ising, double T, double Hmin, double Hmax, double dH, int steplimit, int seed)
	{
		String SPCrun="SPC data <L="+fmt.format(L)+", R="+fmt.format(R)+", T="+fmt.format(T*100)+"p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String SPCpath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling/"+dynamics+"/"+SPCrun+".txt";
		double pb=1-Math.exp((1+0.745/(1-ising.percent))*ising.J/T);;     //need to use the bond probability at spinodal
		
		int spantemp=0;
		for(double h=Hmin; h<Hmax; h+=dH)
		{
		  
			spantemp=findlargestcluster(ising, T, h, steplimit, pb, seed);
			PrintUtil.printlnToFile(SPCpath, h, spantemp);
			
		}
		
	}
	
	public double meanclustersize(IsingStructure ising, double T, double H, int steplimit, int seed)
	{
		
		double meansize=0;
		
		ising.Sinitialization(1, Sseed);
		Random rand= new Random(seed);
		for(int step=0; step<steplimit; step++)
		{
			params.set("T", T);
			params.set("H", -H);
			ising.MCS(T, -H, rand, 1, dynamics);
			
			Job.animate();
			params.set("Emcs", step);
			params.set("magnetization", ising.magnetization);
		}
		
		double pb=1-Math.exp((1+ising.magnetization/(1-ising.percent))*ising.J/T);
		
		clustersize=ising.Clustergrowth(ising.spin, -1, pb, seed, seed, true);
		Job.animate();
		meansize=meansize(clustersize);
		
		return meansize;
	}
	
	
	public int findlargestcluster(IsingStructure ising, double T, double H, int steplimit, double pb, int seed)
	{
		
		int largestsize=0;
		
		String FLCrun="FLC <L="+fmt.format(L)+", R="+fmt.format(R)+", T="+fmt.format(T*100)+", H="+bmt.format(H*1000)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String FLCpath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling/"+dynamics+"/"+FLCrun+".txt";
		String FLCpic = "/Users/liukang2002507/Desktop/simulation/ClusterScaling/"+dynamics+"/FLCpic/cluster";
		
		ising.Sinitialization(1, Sseed);
		Random rand= new Random(seed);
		for(int step=0; step<steplimit; step++)
		{
			params.set("T", T);
			params.set("H", -H);
			ising.MCS(T, -H, rand, 1, dynamics);
			
			Job.animate();
			params.set("Emcs", step);
			params.set("magnetization", ising.magnetization);
		}
		clustersize=ising.Clustergrowth(ising.spin, -1, pb, seed, seed, true);
		Job.animate();
		printdata(FLCpath, clustersize);
		Tools.Picture(grid3, (int)(H*1000), (int)(percent*1000), FLCpic);
		
		int[] clinfotemp=new int [2];
		clinfotemp=ising.ClusterInfo(ising.largestcluster);
		largestsize=clinfotemp[0];
		
		return largestsize;
	}
	
	public void load(Control ClusterScaling)
	{

		ClusterScaling.frameTogether("Display", grid1 ,grid2, grid3);

		params.add("L", 600);
		
		params.add("la",10);    // scale of the bias dilution region
		params.add("lb",10); 
		params.add("R", 30);
		
		params.add("NJ", -4.0);
	    params.add("deadsites");

		params.add("percent", 0.00);
		params.add("biaspercent", 0.00);
		
		 		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));

	
		params.addm("T", 1.778);
		params.addm("H", 1.26);
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
	   
	    //singlerunTc(Istemp, T, 2000, true, 1);
	    //singlerunHs(Istemp, T, H, 200, true, 1);
	    
	    //approachHs(Istemp, T, 1.26, 1.265, 0.001, 100, 1);
	    chi(Istemp, T, 1.20*(1-percent), 1.24*(1-percent), 0.002, 50, 1);
	    
	    
	    Job.animate();

	}
	
}
