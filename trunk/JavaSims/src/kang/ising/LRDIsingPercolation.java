package kang.ising;

import java.awt.Color;
import java.text.DecimalFormat;

import chris.util.PrintUtil;
import chris.util.Random;

import kang.ising.BasicStructure.IsingStructure;

import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.params.DoubleValue;

public class LRDIsingPercolation extends Simulation
{
	Grid grid1=new Grid("grid1");
	Grid grid2=new Grid("grid2");
	
	public int Lmin,Lmax,M,R,deadsites,Dseed,Bseed,Sseed;
	public double percent,biaspercent,NJ;
	public double T,H;
	public IsingStructure IS;
	public IsingStructure Istemp;
	public double d[];  //the sequence of the scanning dilution percentage
	public double qmin,qmax;
	public int qnumber;
	public double reducedM;

	private DecimalFormat fmt = new DecimalFormat("0000");
	private DecimalFormat qmt = new DecimalFormat("000000");
	
	
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
	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
	}
	
	public static void main (String[] LRDIsingPercolation){
		new Control(new LRDIsingPercolation(), "Kang Liu's long range site-diluted ising model's percolation problem" );
	}
	
	public void load(Control LRDIsingPercolation){
		LRDIsingPercolation.frame (grid1);
		LRDIsingPercolation.frame (grid2);

		params.add("Lmin", 512);
		params.add("Lmax", 512);
		params.add("qmin",0);
		params.add("qmax",1.00);
		params.add("qnumber",50);
		
		params.add("L");
		params.add("R",2);
		params.add("NJ",-4.0);	
		params.add("percent", 0.00);
		params.add("biaspercent", 0.00);
		params.add("deadsites");	
		params.add("Dseed",1);
		params.add("Bseed",1);
		params.add("Sseed",1);
		
		params.addm("T", new DoubleValue(99, 0, 100).withSlider());
		params.addm("H", new DoubleValue(0, -2, 2).withSlider());
		
		params.add("MCS");
		params.add("copies");
		params.add("magnetization");
		params.add("reduced M");
	}
	
	public void SaturateM(IsingStructure Ising, int R, int L, int Tq, int copies, int steplimit)
	{
		int run=0;
		String path="/Users/liukang2002507/Desktop/simulation/LRDIP/R="+fmt.format(R)+"-L="+fmt.format(L)+"-run="+fmt.format(run)+".txt";
		
		for(run=1; run<=copies; run++)
		{
				IS.Dinitialization(run, run, 10, 10);
				params.set("copies", run);
	            params.set("deadsites",IS.deadsites);
	            IS.Sinitialization(0, run);
	            Istemp=IS.clone();

	            Job.animate();
	            
	            for(int prestep=0; prestep<50; prestep++)
	            {
	            	params.set("T",99);
	            	params.fget("T");
	            	params.fget("H");
	            	Random heat=new Random(run);
	            	Istemp.MCS(T,H,heat,1);
	            	params.set("MCS", prestep-50);
	            	Job.animate();
	            }
	            params.set("T",0.0001);
	            for(int step=0; step<steplimit; step++)
	            {
	            	params.fget("T");
	            	params.fget("H");
	            	Random flip=new Random(run);
	            	Istemp.MCS(T,H,flip,1);
	            	params.set("MCS", step);
	            	Job.animate();
	            }
	            
		}
		
	}
	
	public void Singlerun(IsingStructure Ising, int R, int L, double Tq, int steplimit, double q)
	{
		int run=0;
		double dilution=Ising.percent;
		String path="/Users/liukang2002507/Desktop/simulation/LRDIP/singlerun-R="+fmt.format(R)+"-L="+fmt.format(L)+"-q=0."+qmt.format(dilution*1000000)+".txt";
		String overall="/Users/liukang2002507/Desktop/simulation/LRDIP/L="+fmt.format(L)+".txt";
		
		for(run=1; run<=1; run++)
		{
				IS.Dinitialization(run, run, 10, 10);
				params.set("copies", run);
	            params.set("deadsites",IS.deadsites);
	            IS.Sinitialization(0, run);
	            Istemp=IS.clone();
	            double mag=0;
	            int positive=0;
	            int negative=0;

	            Job.animate();
	            int presteplimit=50;
	            for(int prestep=0; prestep<presteplimit; prestep++)
	            {
	            	params.set("T",99);
	            	T=params.fget("T");
	            	H=params.fget("H");
	            	Random heat=new Random(run);
	            	Istemp.MCS(T,H,heat,1);
	            	params.set("MCS", prestep-presteplimit);
	            	Job.animate();
	            	mag=Istemp.Magnetization();
	            	reducedM=mag/(1-Istemp.percent);
	            	PrintUtil.printlnToFile(path, prestep-presteplimit, mag, reducedM);
	            	params.set("magnetization",mag);
	            	params.set("reduced M", reducedM);
	            }
	            params.set("T",Tq);
	            for(int step=0; step<steplimit; step++)
	            {
	            	T=params.fget("T");
	            	H=params.fget("H");
	            	Random flip=new Random(run);
	            	Istemp.MCS(T,H,flip,1);
	            	params.set("MCS", step);
	            	mag=Istemp.Magnetization();
	            	reducedM=mag/(1-Istemp.percent);
	            	
	            	params.set("magnetization",mag);
	            	params.set("reduced M", reducedM);
	            	Job.animate();
	            	PrintUtil.printlnToFile(path, step, mag, reducedM);
	            	if(steplimit-step<=100)
	            	{
	            		if(reducedM<0)
	            			negative++;
	            		if(reducedM>0)
	            			positive++;
	            	}
	            	
	                	
	            }
	            PrintUtil.printlnToFile(overall, q, mag, Math.abs(reducedM),positive, negative);
	            
		}
		
	}
	
	
	public void run(){
		
		
		Lmin = (int)params.fget("Lmin");
		Lmax = (int)params.fget("Lmax");
		R=(int)params.fget("R");
		qmin=params.fget("qmin");
		qmax=params.fget("qmax");
		qnumber=(int)params.fget("qnumber");
		
		NJ = params.fget("NJ");
		
		Dseed = (int)params.fget("Dseed");
		Bseed = (int)params.fget("Bseed");
		Sseed = (int)params.fget("Sseed");
		
		d=new double[qnumber];
		
		for(int l=0; l<qnumber; l++)
		{
			d[l]=((qmax-qmin)/qnumber)*l+qmin;
		}
		
		for(int L=Lmin; L<=Lmax; L=L*2)
		{
		    for(int j=0; j<qnumber; j++)
		    {
		    	percent=d[j];
		    	biaspercent=percent;
		    	params.set("percent",percent);
		    	params.set("biaspercent",biaspercent);
		    	IS=new IsingStructure(L,L,R,NJ,percent,biaspercent,"diamond");
	            Istemp=new IsingStructure(L,L,R,NJ,percent,biaspercent,"diamond");

	            Singlerun(IS, R, L, 0.001*(1-percent), 20000, percent);
	            //SaturateM(IS,R,L,0.001,20,1000);
		    }
			
	    
	     


			
			
			clear();
		}

        
	   
	    

	}
	
	
}

