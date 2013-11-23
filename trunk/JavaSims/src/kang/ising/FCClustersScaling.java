package kang.ising;

import java.text.DecimalFormat;

import chris.util.PrintUtil;
import chris.util.Random;

import kang.ising.BasicStructure.FCIsing;
import kang.ising.BasicStructure.BasicTools;
import kang.ising.BasicStructure.IsingStructure;
import kang.ising.BasicStructure.FCPercolation;

import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.Control;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

//plotting tools
import scikit.graphics.dim2.Plot;
import scikit.dataset.Accumulator;


public class FCClustersScaling extends Simulation{
	
	public int L,N,deadsites;

	public double percent,NJ;
	public double T,H;
	public FCIsing IS;
	public FCIsing Istemp;
	public BasicTools Tools;
	public FCPercolation FCP;
	public FCPercolation FCPtemp;
	public double pbN;
	
	//public int progress;
    //public int usefulruns=0;
    
    //public double varianceX;
    //public double varianceC;
    //public double Cvdata[];
    //public double Chidata[];
    
	private DecimalFormat fmt = new DecimalFormat("000");
	private DecimalFormat pmt = new DecimalFormat("0");
	//public double startH=0;
	
	
	public void animate()
	{	
		params.set("magnetization", Istemp.m);
		params.set("up", (int)Istemp.Nu);
		params.set("down", (int)Istemp.Nd);
		params.set("first cluster", (int)FCPtemp.clusters[0]);
		params.set("pbN", pbN);
			
	}
	
	public void clear()
	{
		
	}
	
	public static void main (String[] FCClustersScaling){
		new Control(new FCClustersScaling(), "Kang Liu's fully connected clusters scaling" );
	}
	
	public void load(Control FCCriticalpoint)
	{
		params.add("N");
		params.add("L", 100);
		params.add("NJ",-4.0);	
		params.add("percent", 0.40);
		params.add("livesites");	
	
		params.addm("T", 1.422);
		params.addm("H", 1.25);
		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));
		params.add("MCS");
		params.add("copies");
		
		params.add("magnetization");
		params.add("up");
		params.add("down");
		
		params.add("pbN", 0.0);
		params.add("first cluster");
	}
	
	public void FCPtestrun(FCPercolation fcp, double minpbN, double maxpbN, double dpN)
	{
		for (double pN=minpbN; pN<maxpbN; pN+=dpN)
		{
			pbN=pN;
			fcp.SetPb(pbN);
			fcp.GenerateClusters(1);
			Job.animate();
		}
	}
	
    public void run()
    {		
		L = (int)params.fget("L");
		N = L * L;
	    params.set("N", N);
		NJ = params.fget("NJ");
		percent=params.fget("percent");
		String dynamics= params.sget("Dynamics");
		
	
	    IS=new FCIsing(N);
	    Istemp=new FCIsing(N);
	    FCP=new FCPercolation(N);
	    FCPtemp=new FCPercolation(N);
	    
	    Tools=new BasicTools();
	    
	    IS.dilute(percent);
	    IS.setJ(NJ);
	    
	    IS.initializeDilution(0);
	    params.set("livesites",IS.livesites);
	    
	    
	    //Job.animate();
	    Istemp=IS.clone();
	    params.set("up", Istemp.Nu);
	    params.set("down", Istemp.Nd);
	    
	    
	    
	    FCPtestrun(FCPtemp, 0.1, 2, 0.1);
	    
	}
	
	
}