package kang.ising;


import java.text.DecimalFormat;

import chris.util.PrintUtil;
import chris.util.Random;


import kang.ising.BasicStructure.IsingStructure;
import kang.ising.BasicStructure.IsingStructure1D;
import kang.ising.BasicStructure.BasicTools;

import scikit.jobs.Job;
import scikit.jobs.Simulation;

import scikit.jobs.Control;
import scikit.jobs.params.ChoiceValue;


public class Criticalpoint1D extends Simulation
{
	
	public int L,R,deadsites,Dseed,Bseed,Sseed;
	public double percent,biaspercent,NJ;
	public double T,H;
	public IsingStructure1D IS;
	public IsingStructure1D Istemp;
	public BasicTools Tools;
	
	public int progress;
    public int usefulruns=0;
    
    public double varianceX;
    public double varianceC;
    public double Cvdata[];
    public double Chidata[];
    
	private DecimalFormat fmt = new DecimalFormat("000");
	public double startH=0;
	
	public double[] DataTc(IsingStructure1D ising, double T, double H, int presteplimit, int number, int copies, String dynamics)  //calculate the specific heat of a given system at T,
	{
	    
	    
	    double data[]=new double[2];   //data[0]=Cv, data[1]=Chi
		double tempM[],tempE[];
	    tempM= new double [number];
	    tempE= new double [number];
	    Chidata= new double [copies];
	    Cvdata=new double[copies];
	    double meanX=0;
	    double meanC=0;
	    
	    double totalM=0;
	    double totalE=0;
	    double averageM=0;
	    double averageE=0;
	    String check ="/Users/liukang2002507/Desktop/simulation/1Dising/Criticalpoints/"+dynamics+"/check L="+fmt.format(ising.L)+" R="+fmt.format(ising.R)+" q="+fmt.format(ising.percent*100)+".txt";
	    
	    
	   
		for(int c=0; c<copies;c++)
		{
			Istemp=ising.clone();
			Random cflip=new Random(c);
			params.set("copies", c);
			for(int heat=0; heat<5; heat++)
			{
				Istemp.MCS(9, H, cflip, 1, dynamics);
				Job.animate();
				params.set("MCS", -9999);

			}
			for(int prestep=0; prestep< presteplimit; prestep++)
			{
				Istemp.MCS(T, H, cflip, 1, dynamics);
				Job.animate();
				params.set("MCS", prestep-presteplimit);
			}
			for(int step=0; step<number; step++)
			{
				tempM[step]=Math.abs(Istemp.TotalSpin());
				tempE[step]=Istemp.TotalIntEnergy();
				
				Istemp.MCS(T, H, cflip, 1, dynamics);
				totalM+=tempM[step];
				totalE+=tempE[step];
				Job.animate();
				params.set("MCS", step);
			}
			averageM=totalM/number;
			averageE=totalE/number;
			PrintUtil.printlnToFile(check, T, c, averageM, averageE);
			Chidata[c]=IS.Fluctuation(tempM, number);
			Cvdata[c]=IS.Fluctuation(tempE, number);
		}
		meanX=Tools.Mean(Chidata, copies);
		meanC=Tools.Mean(Cvdata, copies);
		varianceX=Tools.SD(Chidata, copies, meanX);
		varianceC=Tools.SD(Cvdata, copies, meanC);
		data[0]=meanC;
		data[1]=meanX;
		
		return data;
		
	}

	
	public void Criticalpoints(IsingStructure1D ising, double Tmax, double Tmin, double increment, double targetT, int limit, int number, int copies, String dynamics)
	{
		String path="/Users/liukang2002507/Desktop/simulation/1Dising/Criticalpoints/"+dynamics+"/Tc data L="+fmt.format(ising.L)+" R="+fmt.format(ising.R)+" q="+fmt.format(ising.percent*100)+".txt";
		
		
		for(double t=Tmax; t>Tmin; t-=increment)
		{
			int prelimit=limit;
			//prelimit=(int)(Math.sqrt((Tmax-targetT)/(t-targetT))*limit);
			
			params.set("T", t);
			params.set("H", 0);
			double data[]=new double[2];
			data=DataTc(ising, t, 0, prelimit, number, copies, dynamics);
			
			PrintUtil.printlnToFile(path, t, data[0], varianceC, data[1], varianceX);
		}
	}
	
	
	public double[] DataHs(IsingStructure1D ising, double t, double h, int presteplimit, int number, int copies, String dynamics)  //calculate the specific heat of a given system at T,
	{
	    
	    
        double data[]= new double[2];      //data[0] is Cv and data[1] is chi
		
	    double tempE[];
	    double tempM[];
	    double usedE[];
	    double usedM[];
	    tempE= new double [number];
	    tempM= new double [number];
	    usedE= new double [number-200-200];
	    usedM= new double [number-200-200];
	    double totalCv=0;
	    double totalX=0;
        usefulruns=0;
        double tempCv[];
        double tempX[];
        tempCv= new double[copies];
        tempX = new double[copies];
        
        
        
        
		for(int c=0; c<copies;c++)
		{
			Istemp=ising.clone();
			Random cflip=new Random(c);
			params.set("copies", c);
	
			
			int endofstep=-1;
			int lifetime=2100;
			
			for(int heat=0; heat<5; heat++)
			{
				Istemp.MCS(9, h, cflip, 1, dynamics);
				Job.animate();
				params.set("MCS", -9999);

			}
		
			for(int prestep=0; prestep< presteplimit; prestep++)
			{
				Istemp.MCS(t, h, cflip, 1, dynamics);
				Job.animate();
				params.set("MCS", prestep-presteplimit);
			}
			params.set("H",-h);
			double field=-h;

			
			for(int step=0; (endofstep<0)&(step<number); step++)
			{
				tempE[step]=Istemp.TotalEnergy(field);
		        tempM[step]=Istemp.TotalSpin();
		        
				if(tempM[step]<0)
					{
					endofstep=1;
					lifetime=step;
					}
				Istemp.MCS(t, field, cflip, 1, dynamics);
				Job.animate();
				params.set("MCS", step);
				params.set("magnetization", tempM[step]/L);
			}
			
			if(lifetime>2000)
			{
				for(int tt=0;tt<number-200-200;tt++)
				{
					
					usedE[tt]=tempE[tt+200];
					usedM[tt]=tempM[tt+200];
				}
				
				tempCv[usefulruns]=IS.Fluctuation(usedE, number-200-200);
				tempX[usefulruns]=IS.Fluctuation(usedM, number-200-200);
				totalCv+=tempCv[usefulruns];
				totalX+=tempX[usefulruns];
				
				
				PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/1Dising/Spinodals/"+dynamics+"/check L="+fmt.format(ising.L)+" R="+fmt.format(ising.R)+" q="+fmt.format(ising.percent*100)+".txt",h, c);
				usefulruns++;
				
			}
		}
		double averageCv=totalCv/usefulruns;
		double averageX=totalX/usefulruns;
		double totalCv2=0;
		double totalX2=0;
		for(int u=0;u<usefulruns; u++)
		{
			totalCv2+=(tempCv[u]-averageCv)*(tempCv[u]-averageCv);
			totalX2+=(tempX[u]-averageX)*(tempX[u]-averageX);
		}
		varianceC=Math.sqrt(totalCv2/usefulruns);
		varianceX=Math.sqrt(totalX2/usefulruns);
		
		
		data[0]=averageCv;
		data[1]=averageX;
		
		return data;
	}
	
 	public void Spinodals(IsingStructure1D ising, double Hmin, double Hmax, double dH, String dynamics)
 	{
 		double[] data=new double[2];
 		for(double h=Hmax; h>Hmin; h-=dH)
 		{
			params.set("H", h);
 			data=DataHs(ising,T,h,200,2000,20, dynamics);
 			PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/1Dising/Spinodals/"+dynamics+"/usefulruns q=0."+fmt.format(percent*1000)+" L="+fmt.format(ising.L)+" R="+fmt.format(ising.R)+".txt","H=    ",(double)h );
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/1Dising/Spinodals/"+dynamics+"/Hs data L="+fmt.format(ising.L)+" R="+fmt.format(ising.R)+" q="+fmt.format(ising.percent*100)+".txt";
 			PrintUtil.printlnToFile(SaveAs, h, data[0], varianceC, data[1], varianceX, usefulruns);
 		} 
 	}
	
 	public void scanHs(IsingStructure1D ising, double min, double max, double dh, String dynamics)    //for the 1st search
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
 	
	public double Susceptibility(IsingStructure1D ising, double T, double H, int presteplimit, int number, int copies, String dynamics)  //calculate the specific heat of a given system at T,
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
				Istemp.MCS(9, H, cflip, 1, dynamics);
				Job.animate();
				params.set("MCS", -9999);

			}
		
			for(int prestep=0; prestep< presteplimit; prestep++)
			{
				Istemp.MCS(T, H, cflip, 1, dynamics);
				Job.animate();
				params.set("MCS", prestep-presteplimit);
			}
			params.set("H",-H);
			double field=-H;

			
			for(int step=0; (endofstep<0)&(step<number); step++)
			{
				
				tempM[step]=Istemp.TotalSpin();
				if(tempM[step]<0)
					{
					endofstep=1;
					lifetime=step;
					}
				Istemp.MCS(T, field, cflip, 1, dynamics);
				Job.animate();
				params.set("MCS", step);
				params.set("magnetization", tempM[step]/L);
			}
			
			if(lifetime>2000)
			{
				for(int t=0;t<number-200-200;t++)
				{
					usedM[t]=tempM[t+200];
				}
				tempX[usefulruns]=IS.Fluctuation(usedM, number-200-200);
				totalX+=tempX[usefulruns];
				//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/Hs/usefulrunsq=0."+fmt.format(percent*1000)+".txt",H, c);
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
 	
 	
	public void animate()
	{
		
		params.set("magnetization", Istemp.magnetization);
		params.set("intenergy", Istemp.totalintenergy);
		
		
	}
	
	
	public void clear()
	{
		
	}
	
	public static void main (String[] Criticalpoint1D){
		new Control(new Criticalpoint1D(), "Kang Liu's critical points and spinodal in 1D" );
	}
	
	public void load(Control Criticalpoint1D){
		

		params.add("L", 4000);

		params.add("R", 200);
		params.add("NJ",-4.0);	
		params.add("percent", 0.00);
		params.add("biaspercent", 0.00);
		params.add("deadsites");	
		params.add("Dseed",2);
		params.add("Bseed",2);
		params.add("Sseed",1);
		
		params.addm("T", 1.778);
		params.addm("H", 0.0);
		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));
		
		params.add("MCS");
		params.add("copies");
		params.add("magnetization");
		params.add("intenergy");
	}
	
	public void run(){
		
		
		L = (int)params.fget("L");
	
		
		
		R = (int)params.fget("R");
		NJ = params.fget("NJ");

		percent=params.fget("percent");
		biaspercent=params.fget("biaspercent");
		String dynamics= params.sget("Dynamics");
		
		Dseed = (int)params.fget("Dseed");
		Bseed = (int)params.fget("Bseed");
		Sseed = (int)params.fget("Sseed");
		
	    IS=new IsingStructure1D(L,R,NJ,percent,biaspercent);
	    Istemp=new IsingStructure1D(L,R,NJ,percent,biaspercent);
	    Tools=new BasicTools();
	    
	    IS.Dinitialization(Dseed, Bseed, 10);
	    params.set("deadsites",IS.deadsites);
	    IS.Sinitialization(0, Sseed);
	    
	    Job.animate();
	    //findTc(IS,3.52,3.59,0.001, dynamics);
	    //findTcviaX(IS,3.90,4.10,0.005, dynamics);
        
	   {
	    	T=params.fget("T");
	    	
	    	scanHs(IS,0,1.270*(1-IS.percent),0.002, dynamics);    //startH gives 5 more runs in copies runs
	 	    
	 	    //startH=1.065;    
	 	    
	 	    //findHs(IS,startH-0.1,startH,0.002, dynamics);
	       
	 	    //findHsCv(IS,startH-0.1,startH,0.01, dynamics);
	 	    
	 	    Spinodals(IS,startH-0.2,startH,0.002, dynamics);
	    }
	    
	    
	   
	    //Criticalpoints(IS, 4.20*(1-IS.percent), 3.60*(1-IS.percent), 0.01, 4, 2000, 2000, 10, dynamics);
	    
	    
		//CriticalpointsCv(IS, 3.96, 3.905, 0.005, 4, 2000, 2000, 10, dynamics);
	    
	    //HSboundary(IS, T, 0.900, 1.090, 0.001, dynamics);
	   
	    //ScanHsBoundary(0.60, 1.00, 0.05, 0.005, dynamics);
	    
	    

	}
	
}
	