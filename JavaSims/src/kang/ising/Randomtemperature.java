package kang.ising;

import java.awt.Color;
import java.text.DecimalFormat;

import chris.util.PrintUtil;
import chris.util.Random;

//import kang.ising.BasicStructure.IsingStructure;

import kang.ising.BasicStructure.MeanFieldIsingStructure;

import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.params.DoubleValue;


public class Randomtemperature extends Simulation
{
	Grid grid1=new Grid("grid1");
	Grid grid2=new Grid("grid2");
	public int L1,L2,M,Sseed;
	public int Tseed;
	public double NJ;
	public double T,H;
	public MeanFieldIsingStructure MFIS;
	public MeanFieldIsingStructure Istemp;

	private DecimalFormat fmt = new DecimalFormat("000");

	public double SpecificHeat(MeanFieldIsingStructure ising, double T, double H, int presteplimit, int number, int copies)  //calculate the specific heat of a given system at T,
	{
	    double tempE[];
	    tempE= new double [number];
	    double totalC=0;
	    
		for(int c=0; c<copies;c++)
		{
			Istemp=ising.clone();
			Random cflip=new Random(c);
			params.set("copies", c);
			for(int heat=0; heat<5; heat++)
			{
				Istemp.MCS(9, H, cflip, 1);
				Job.animate();
				params.set("MCS", -9999);

			}
			for(int prestep=0; prestep< presteplimit; prestep++)
			{
				Istemp.MCS(T, H, cflip, 1);
				Job.animate();
				params.set("MCS", prestep-presteplimit);
			}
			for(int step=0; step<number; step++)
			{
				
				tempE[step]=Istemp.TotalIntEnergy();
				PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/RandomT/progress.txt", T, c, step);
				Istemp.MCS(T, H, cflip, 1);
				Job.animate();
				params.set("MCS", step);
			}
			totalC+=Istemp.Fluctuation(tempE, number);
		}
		return totalC/copies;
		
	}
	
 	public void findTc(MeanFieldIsingStructure ising,double Tmin, double Tmax, double dT)
 	{
 		double Cv=0;
 		for(double t=Tmin; t<Tmax; t+=dT)
 		{
			params.set("T", t);
 			Cv=SpecificHeat(ising,t,0,5000,2000,10);
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/RandomT/Cv.txt";
 			PrintUtil.printlnToFile(SaveAs, t, Cv);
 		}
 	}
	
 	
 	public void Heatup(MeanFieldIsingStructure ising, int steplimit, Random flip)
 	{
	    for(int heat=0; heat<steplimit; heat++)
	    {
	    	MFIS.MCS(99, 0, flip, 1);
	    	Job.animate();
			params.set("MCS", heat-steplimit);
	    	params.set("magnetization", MFIS.Magnetization());
	    }
 	}
	
	
	public void EvolutionwithNoise(MeanFieldIsingStructure ising, double T, double dt, int copies, int steplimit)
	{
		//double Mevolution[][]=new double[copies][steplimit+1];
		//double RandomT[][]=new double[copies][steplimit+1];
		
		for(int r=0; r<copies; r++)
		{
			Tseed=r+1;
			Random Trand=new Random (Tseed);
			Random Srand=new Random (Tseed);
			Istemp=ising.clone();
			params.set("copies",r+1);
			double m=Istemp.Magnetization();
			
			PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/RandomT/noise1/run#"+fmt.format(r+1)+".txt","T=    ",(double)T );
			PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/RandomT/noise1/run#"+fmt.format(r+1)+".txt","dt=    ",(double)dt );
		    PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/RandomT/noise1/run#"+fmt.format(r+1)+".txt","initialm=    ",(double)m );
			
		    for(int step=0; step<steplimit; step++)
			{
				double t=T+dt*(Trand.nextDouble()-0.5);
				Istemp.MCS(t, 0, Srand, 1);
				m=Istemp.Magnetization();
				//Mevolution[r][step]=m;
				//RandomT[r][step]=t;
				
				Job.animate();
				
				params.set("MCS", step+1);
				params.set("magnetization", m);
				params.set("T",t);
				
				PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/RandomT/noise1/run#"+fmt.format(r+1)+".txt", step ,(double)m, (double)t);
			}
			
		}
	}
	
	public void EvolutionwithNoise2(MeanFieldIsingStructure ising, double T, double dt, int copies, int steplimit)
	{
		//double Mevolution[][]=new double[copies][steplimit+1];
		//double RandomT[][]=new double[copies][steplimit+1];
		
		for(int r=0; r<copies; r++)
		{
			Tseed=r+1;
			Random Trand=new Random (Tseed);
			Random Srand=new Random (Tseed);
			Istemp=ising.clone();
			params.set("copies",r+1);
			double m=Istemp.Magnetization();
			
			PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/RandomT/noise2/run#"+fmt.format(r+1)+".txt","T=    ",(double)T );
			PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/RandomT/noise2/run#"+fmt.format(r+1)+".txt","dt=    ",(double)dt );
		    PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/RandomT/noise2/run#"+fmt.format(r+1)+".txt","initialm=    ",(double)m );
			
		    for(int step=0; step<steplimit; step++)
			{
		    	double t=0;
				if(Trand.nextDouble()>=0.5)
					t=T+dt/2;
				else
					t=T-dt/2;
				Istemp.MCS(t, 0, Srand, 1);
				m=Istemp.Magnetization();
				//Mevolution[r][step]=m;
				//RandomT[r][step]=t;
				
				Job.animate();
				
				params.set("MCS", step+1);
				params.set("magnetization", m);
				params.set("T",t);
				
				PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/RandomT/noise2/run#"+fmt.format(r+1)+".txt", step ,(double)m, (double)t);
			}
			
		}
	}
	
	
	
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
		grid1.registerData(MFIS.L1, MFIS.L2, MFIS.spin);
		grid2.setColors(ising);
		grid2.registerData(Istemp.L1, Istemp.L2, Istemp.spin);
	}

	public void clear()
	{
		grid1.clear();
		grid2.clear();
	}
	
	public static void main (String[] Randomtemperature){
		new Control(new Randomtemperature(), "Kang Liu's Random temperature Ising model" );
	}

	public void load(Control Randomtemperature){
		Randomtemperature.frame (grid1);
		Randomtemperature.frame (grid2);

		params.add("L1", 800);
		params.add("L2", 800);
		params.add("NJ",-4.0);	

		params.add("Sseed",1);
		
		params.addm("T", new DoubleValue(3, 0, 20).withSlider());
		params.addm("H", new DoubleValue(0, -2, 2).withSlider());
				
		params.add("copies");
		params.add("MCS");

		params.add("magnetization");
	}

	
	public void run(){
		
		
		L1 = (int)params.fget("L1");
		L2 = (int)params.fget("L2");
		M = L1 * L2;
		
		NJ = params.fget("NJ");

		Sseed = (int)params.fget("Sseed");
		
	    MFIS=new MeanFieldIsingStructure(L1,L2,NJ);
	    for(int ii=0; ii<MFIS.M; ii++)
	    {
	    	MFIS.spin[ii]=1;
	    }
	    Istemp=new MeanFieldIsingStructure(L1,L2,NJ);
	    Random flip= new Random(1);
	    

	    
	    MFIS.Sinitialization(0, Sseed);
	    
	    //findTc(MFIS,3.955, 3.965, 0.001);
	    Job.animate();
	    
	    Heatup(MFIS, 100, flip);
	    

	    //EvolutionwithNoise(MFIS, 3.9, 0.6, 20, 1000);
	    EvolutionwithNoise2(MFIS, 3.9, 0.6, 20, 1000);
	    
	    Job.animate();

    
	    

	}
	
	
	
	
	
	
	
}

