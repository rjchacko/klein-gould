package kang.ising;


import java.text.DecimalFormat;

import chris.util.PrintUtil;
import chris.util.Random;

import kang.ising.BasicStructure.FCIsing;
import kang.ising.BasicStructure.BasicTools;
import kang.ising.BasicStructure.IsingStructure;

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


	public double Susceptibility(FCIsing ising, double T, double h, int presteplimit, int number, int copies, String dynamics)  //calculate the susceptibility of a given system at T H,
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
				Istemp.MCS(dynamics, cflip, cflip, 99, h,  1);
				Job.animate();
				params.set("MCS", -9999);
			}
		
			for(int prestep=0; prestep< presteplimit; prestep++)
			{
				Istemp.MCS(dynamics, cflip, cflip, T, h,  1);
				
				Job.animate();
				params.set("MCS", prestep-presteplimit);
			}
			params.set("H",-h);
			double field=-h;

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
				PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/FCHs/usefulrunsq=0."+fmt.format(percent*1000)+".txt",h, c);
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
	
 	

 	
 	public void findHs(FCIsing ising,double Hmin, double Hmax, double dH, String dynamics)
 	{
 		double X=0;
 		for(double h=Hmax; h>Hmin; h-=dH)
 		{
			params.set("H", h);
 			X=Susceptibility(ising,T,h,200,2000,20, dynamics);
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
 			X=Susceptibility(ising,T,h,100,2000,20, dynamics);
 			if(usefulruns>=10)
 				{
 				startH=h;
 				end=1;
 				}
 			if(h<min)
 				end=1;
 		}
 	}
 	
	public int count(FCIsing ising, double t, double h, int presteplimit, String dynamics)
	{
		
		int step;
		int count=0;
		for(int c=0; c<10; c++)
		{
			step=0;
			
			Istemp=ising.clone();
			Random cflip=new Random(c);
			params.set("copies", c+1);
			Istemp.initialize(0);
			
			params.set("H", h);
			params.set("T", t);
			
			for(int prestep=0; prestep< presteplimit; prestep++)
			{
				
				Istemp.MCS(dynamics, cflip, cflip, t, h, 1 );
				
				Job.animate();
				params.set("MCS", prestep-presteplimit);
			}
			
			params.set("H", -h);
			params.set("T", t);
			
			for(step=0; (Istemp.m>0)&&(step<3500); step++)
			{
				
				Istemp.MCS(dynamics, cflip, cflip, t, -h, 1 );
				Job.animate();
				params.set("MCS", step);
			}
			if(Istemp.m>0)
				count++;
		}
		
		return count;
	}
	
	public double[] Hfirstfive(FCIsing ising, double t, double startfield, double dh, int presteplimit, String dynamics)
	{
		double count=0;
		double targetfield=0;
		double output[]=new double[2];  //output[0]=targetfield, output[1]=count
		
		double h;
		for(h=startfield; count<6; h-=dh)
		{
			params.set("H", h);
			count=count(ising, t, h, presteplimit, dynamics);
		}
		targetfield=h;
		output[0]=targetfield;
		output[1]=count;
		
		
		return output;
		
	}
	
	public void ScanHsBoundary(double pmin, double pmax, double dp, double dh, String dynamics)
	{
		double temp, startfield;
		for(double p=pmax; p>pmin; p-=dp)
		{
			params.set("percent", 1-p);
			IS=new FCIsing(N);
		    IS.dilute(1-p);
		    IS.setJ(NJ);
		    
		    IS.initializeDilution(0);
		    params.set("livesites",IS.livesites);
		    
		   
		    temp=1.778*p;
		    startfield=1.273*p;
		    
		    params.set("T",temp);
		    params.set("H",startfield);
		    
		    Job.animate();
		    
		    double output[]=new double[2];
		    output=Hfirstfive(IS, temp, startfield, dh, 200, dynamics);   //ten runs to find which magnetic field provide 5 runs with more than 3500 MCS (pass zero)
		    
		    String SaveAs = "/Users/liukang2002507/Desktop/simulation/FCHs/first5.txt";
 			PrintUtil.printlnToFile(SaveAs, p, temp, output[0], output[1]);
		
		
		}
	}
	
	
	public void HSboundary(FCIsing ising, double T, double Hmin, double Hmax, double dH, String dynamics)
	{
		double lifetime[]=new double[2];
 		for(double h=Hmax; h>Hmin; h-=dH)
 		{
			params.set("H", h);
 			lifetime=lifetime(ising,T,h,200, 10, dynamics);
 			
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/FCHs/lifetime q=0."+fmt.format(percent*1000)+".txt";
 			PrintUtil.printlnToFile(SaveAs, h, lifetime[0], lifetime[1]);
 		} 
	}
	
	
	public double[] lifetime(FCIsing ising, double t, double h, int presteplimit, int copies, String dynamics)
	{
		double lifetime[]=new double[2];
		lifetime[0]=0;
		lifetime[1]=0;
		
		int step;
		
		double temp[]=new double[copies];
		
		
		for(int c=0; c<copies; c++)
		{
			step=0;
			
			Istemp=ising.clone();
			Random cflip=new Random(c);
			params.set("copies", c);
			
			Istemp.initializeDilution(0);
			for(int prestep=0; prestep< presteplimit; prestep++)
			{
				Istemp.MCS(dynamics, cflip, cflip, t, h, 1 );
				Job.animate();
				params.set("MCS", prestep-presteplimit);
			}
			for(step=0; Istemp.m>0; step++)
			{
				Istemp.MCS(dynamics, cflip, cflip, t, -h, 1 );
				Job.animate();
				params.set("MCS", step);
			}
			temp[c]=step;
		}
		
		
		lifetime[0]=Tools.Mean(temp, copies);
		lifetime[1]=Tools.SD(temp, copies, lifetime[0]);
		return lifetime;
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
		params.add("L", 400);
		params.add("NJ",-4.0);	
		params.add("percent", 0.90);
		params.add("livesites");	
	
		params.addm("T", 0.178);
		params.addm("H", 1.25);
		
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
	    
	    
	    
	    //testrun(IS, dynamics);
	
	    //findTcviaX(IS,3.90,4.10,0.005, dynamics);
        
	    T=params.fget("T");
	    scanHs(IS,1.24*(1-percent),1.27*(1-percent),0.0002, dynamics);
	    //startH=1.010;
	    findHs(IS,startH-0.1,startH,0.0002, dynamics);
      
		//CriticalpointsCv(IS, 4.00, 3.70, 0.01, 4, 2000, 2000, 5, dynamics);
	    
	    //ScanHsBoundary(0.20, 1.0, 0.05, 0.0005, dynamics);
	    
	   

	}
	



}