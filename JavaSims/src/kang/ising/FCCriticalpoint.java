package kang.ising;


import java.text.DecimalFormat;

import chris.util.PrintUtil;
import chris.util.Random;

import kang.ising.BasicStructure.FCIsing;
import kang.ising.BasicStructure.BasicTools;

import scikit.jobs.Job;
import scikit.jobs.Simulation;

import scikit.jobs.Control;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

public class FCCriticalpoint extends Simulation
{
	
	public int L,N,deadsites;

	public double percent,NJ;
	public double T,H;
	public FCIsing IS;
	public FCIsing Istemp;
	public BasicTools Tools;
	
	public int progress;
    public int usefulruns=0;
    
    public double varianceX;
    public double varianceC;
    public double Cvdata[];
    
	private DecimalFormat fmt = new DecimalFormat("000");
	public double startH=0;
	
	/*public double XforTc(FCIsing ising, double T, double H, int presteplimit, int number, int copies, String dynamics)  //calculate the specific heat of a given system at T,
	{
	    double tempM[];
	    tempM= new double [number];
	    double totalX=0;
	    
		for(int c=0; c<copies;c++)
		{
			Istemp=ising.clone();
			Random cflip=new Random(c);
			params.set("copies", c);
			for(int heat=0; heat<5; heat++)
			{
				Istemp.MCS(dynamics, cflip, cflip,9, H, 1);
				Job.animate();
				params.set("MCS", -9999);

			}
			for(int prestep=0; prestep< presteplimit; prestep++)
			{
				Istemp.MCS(dynamics, cflip, cflip, T, H,  1);
				Job.animate();
				params.set("MCS", prestep-presteplimit);
			}
			for(int step=0; step<number; step++)
			{
				
				tempM[step]=Istemp.m;
				
				Istemp.MCS(dynamics, cflip, cflip, T, H,  1);
				Job.animate();
				params.set("MCS", step);
			}
			totalX+=IS.Fluctuation(tempM, number);
		}
		return totalX/copies;
		
	}
	*/


	public double Susceptibility(FCIsing ising, double T, double H, int presteplimit, int number, int copies, String dynamics)  //calculate the susceptibility of a given system at T H,
	{
	    double tempM[];
	    double usedM[];
	    tempM= new double [number];
	    usedM= new double [number-200-200];
	    double totalX=0;
        usefulruns=0;
        double tempX[];
        tempX= new double[copies];
	    
		for(int c=0; c<copies;c++)
		{
			Istemp=ising.clone();
			Random cflip=new Random(c);
			params.set("copies", c);
			
			int endofstep=-1;
			int lifetime=2100;
			
			for(int heat=0; heat<5; heat++)
			{
				Istemp.MCS(dynamics, cflip, cflip, 99, H,  1);
				Job.animate();
				params.set("MCS", -9999);
			}
		
			for(int prestep=0; prestep< presteplimit; prestep++)
			{
				Istemp.MCS(dynamics, cflip, cflip, T, H,  1);
				
				Job.animate();
				params.set("MCS", prestep-presteplimit);
			}
			params.set("H",-H);
			double field=-H;

			for(int step=0; (endofstep<0)&(step<number); step++)
			{
				
				tempM[step]=Istemp.m;
				if(tempM[step]<0)
					{
					endofstep=1;
					lifetime=step;
					}
				Istemp.MCS(dynamics, cflip, cflip, T, field,  1);
				
				Job.animate();
				params.set("MCS", step);
				params.set("magnetization", tempM[step]);
			}
			
			if(lifetime>2000)
			{
				for(int t=0;t<number-200-200;t++)
				{
					usedM[t]=tempM[t+200];
				}
				tempX[usefulruns]=IS.Fluctuation(usedM, number-200-200);
				totalX+=tempX[usefulruns];
				PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/FCHs/usefulrunsq=0."+fmt.format(percent*1000)+".txt",H, c);
				usefulruns++;
			}
		}
		double averageX=totalX/usefulruns;
		double totalX2=0;
		for(int u=0;u<usefulruns; u++)
		{
			totalX2+=(tempX[u]-averageX)*(tempX[u]-averageX);
		}
		varianceX=Math.sqrt(totalX2/usefulruns);
		
		return averageX;
		
	}
	
 	
 	public void findTcviaX(FCIsing ising,double Tmin, double Tmax, double dT, String dynamics)
 	{
 		double Cv=0;
 		for(double t=Tmin; t<Tmax; t+=dT)
 		{
			params.set("T", t);
 			Cv=XforTc(ising,t,0,3000,1000,5, dynamics);
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/FCTc/q=0."+fmt.format(percent*1000)+".txt";
 			PrintUtil.printlnToFile(SaveAs, t, Cv);
 		}
 	}
 	
 	public void findHs(FCIsing ising,double Hmin, double Hmax, double dH, String dynamics)
 	{
 		double X=0;
 		for(double h=Hmax; h>Hmin; h-=dH)
 		{
			params.set("H", h);
 			X=Susceptibility(ising,T,h,200,2000,10, dynamics);
 			PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/FCHs/usefulrunsq=0."+fmt.format(percent*1000)+".txt","H=    ",(double)h );
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/FCHs/q=0."+fmt.format(percent*1000)+".txt";
 			PrintUtil.printlnToFile(SaveAs, h, X, varianceX, usefulruns);
 		} 
 	}
 	
 	public void scanHs(FCIsing ising, double min, double max, double dh, String dynamics)    //for the 1st search
 	{   
 		double X=0;
 		int end=-1;
 		
 		for(double h=max; end<0; h-=dh)
 		{
 			params.set("H", h);
 			X=Susceptibility(ising,T,h,100,2000,10, dynamics);
 			if(usefulruns>=5)
 				{
 				startH=h;
 				end=1;
 				}
 			if(h<min)
 				end=1;
 		}
 	}
 	
 	public void testrun(FCIsing ising, String dynamics)
 	{
 		Istemp=ising.clone();
 		Job.animate();
 		for(int t=0; t<999999999;t++)
 		{
 			H=params.fget("H");
 			T=params.fget("T");
 			
 			params.set("MCS",t+1);
 			Random rand=new Random(9999);
 			Istemp.MCS(dynamics, rand, rand, T, H, 1);
 			PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/FCHs/testrun.txt",t, T, Istemp.Nu, Istemp.Nd );
 			Job.animate();
 		}
 	}
	
	public void animate()
	{	
		params.set("magnetization", Istemp.m);
		params.set("up", (int)Istemp.Nu);
		params.set("down", (int)Istemp.Nd);
	}

	public void clear()
	{
		
	}
	
	public static void main (String[] FCCriticalpoint){
		new Control(new FCCriticalpoint(), "Kang Liu's fully connected critical points and spinodal" );
	}

	public void load(Control FCCriticalpoint){
		

		params.add("N");
		params.add("L", 200);
		params.add("NJ",-4.0);	
		params.add("percent", 0.0);
		params.add("livesites");	
	
		params.addm("T", 1.778);
		params.addm("H", 0.0);
		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));
		
		params.add("MCS");
		params.add("copies");
		
		params.add("magnetization");
		params.add("up");
		params.add("down");
	
	}

	
	public void run(){
		
		
		L = (int)params.fget("L");
		
		N = L * L;
		
		
	    params.set("N", N);
	    
		NJ = params.fget("NJ");

		percent=params.fget("percent");
		
		String dynamics= params.sget("Dynamics");
		
	
		
	    IS=new FCIsing(N);
	    Istemp=new FCIsing(N);
	    Tools=new BasicTools();
	    
	    IS.dilute(percent);
	    IS.setJ(NJ);
	    
	    IS.initializeDilution(0);
	    params.set("livesites",IS.livesites);
	    
	    
	    //Job.animate();
	    Istemp=IS.clone();
	    
	    params.set("up", Istemp.Nu);
	    params.set("down", Istemp.Nd);
	    
	    
	    
	    testrun(IS, dynamics);
	
	    //findTcviaX(IS,3.90,4.10,0.005, dynamics);
        
	    //T=params.fget("T");
	    //scanHs(IS,1,1.27,0.01, dynamics);
	    //startH=1.010;
	    //findHs(IS,startH-0.5,startH,0.005, dynamics);
      
		//CriticalpointsCv(IS, 4.00, 3.70, 0.01, 4, 2000, 2000, 5, dynamics);
	    
	    
	   

	}
	



}