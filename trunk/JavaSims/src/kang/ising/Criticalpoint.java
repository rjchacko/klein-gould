package kang.ising;

import java.awt.Color;
import java.text.DecimalFormat;

import chris.util.PrintUtil;
import chris.util.Random;

import kang.ising.BasicStructure.IsingStructure;
import kang.ising.BasicStructure.BasicTools;

import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

public class Criticalpoint extends Simulation
{
	Grid grid1=new Grid("grid1");
	Grid grid2=new Grid("grid2");
	public int L1,L2,M,R,deadsites,Dseed,Bseed,Sseed;
	public double percent,biaspercent,NJ;
	public double T,H;
	public IsingStructure IS;
	public IsingStructure Istemp;
	public BasicTools Tools;
	
	public int progress;
    public int usefulruns=0;
    
    public double varianceX;
    public double varianceC;
    public double Cvdata[];
    
	private DecimalFormat fmt = new DecimalFormat("000");
	public double startH=0;
	
	public double XforTc(IsingStructure ising, double T, double H, int presteplimit, int number, int copies, String dynamics)  //calculate the specific heat of a given system at T,
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
				
				tempM[step]=Istemp.TotalSpin();
				//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/Tc/progress.txt", T, c, step);
				Istemp.MCS(T, H, cflip, 1, dynamics);
				Job.animate();
				params.set("MCS", step);
			}
			totalX+=IS.Fluctuation(tempM, number);
		}
		return totalX/copies;
		
	}
	
	public double SpecificHeat(IsingStructure ising, double T, double H, int presteplimit, int number, int copies, String dynamics)  //calculate the specific heat of a given system at T,
	{
	    double tempE[];
	    tempE= new double [number];
	    Cvdata= new double [copies];
	    double meanC=0;
	    
	    double totalM=0;
	    double averageM=0;
	    String check ="/Users/liukang2002507/Desktop/simulation/CriticalpointsCv/"+dynamics+"/check L="+fmt.format(ising.L1)+" q="+fmt.format(ising.percent*100)+".txt";
	    
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
				
				tempE[step]=Istemp.TotalIntEnergy();
				//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/Tc/progress.txt", T, c, step);
				Istemp.MCS(T, H, cflip, 1, dynamics);
				totalM+=Istemp.magnetization;
				Job.animate();
				params.set("MCS", step);
			}
			averageM=totalM/number;
			PrintUtil.printlnToFile(check, T, c, averageM);
			Cvdata[c]=IS.Fluctuation(tempE, number);
		}
		meanC=Tools.Mean(Cvdata, copies);
		varianceC=Tools.SD(Cvdata, copies, meanC);
		return meanC;
		
	}
	
	public void CriticalpointsCv(IsingStructure ising, double Tmax, double Tmin, double increment, double targetT, int limit, int number, int copies, String dynamics)
	{
		String path="/Users/liukang2002507/Desktop/simulation/CriticalpointsCv/"+dynamics+"/L="+fmt.format(ising.L1)+" q="+fmt.format(ising.percent*100)+".txt";
		for(double t=Tmax; t>Tmin; t-=increment)
		{
			int prelimit=limit;
			//prelimit=(int)(Math.sqrt((Tmax-targetT)/(t-targetT))*limit);
			
			params.set("T", t);
			params.set("H", 0);
			double c=SpecificHeat(ising, t, 0, prelimit, number, copies, dynamics);
			
			PrintUtil.printlnToFile(path, t, c, varianceC);
		}
	}
	
	public double Susceptibility(IsingStructure ising, double T, double H, int presteplimit, int number, int copies, String dynamics)  //calculate the specific heat of a given system at T,
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
				params.set("magnetization", tempM[step]/M);
			}
			
			if(lifetime>2000)
			{
				for(int t=0;t<number-200-200;t++)
				{
					usedM[t]=tempM[t+200];
				}
				tempX[usefulruns]=IS.Fluctuation(usedM, number-200-200);
				totalX+=tempX[usefulruns];
				PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/Hs/usefulrunsq=0."+fmt.format(percent*1000)+".txt",H, c);
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
	

 	public void findTc(IsingStructure ising,double Tmin, double Tmax, double dT, String dynamics)
 	{
 		double Cv=0;
 		for(double t=Tmin; t<Tmax; t+=dT)
 		{
			params.set("T", t);
 			Cv=SpecificHeat(ising,t,0,3000,1000,5, dynamics);
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/Tc/q=0."+fmt.format(percent*1000)+".txt";
 			PrintUtil.printlnToFile(SaveAs, t, Cv);
 		}
 	}
 	
 	public void findTcviaX(IsingStructure ising,double Tmin, double Tmax, double dT, String dynamics)
 	{
 		double Cv=0;
 		for(double t=Tmin; t<Tmax; t+=dT)
 		{
			params.set("T", t);
 			Cv=XforTc(ising,t,0,3000,1000,5, dynamics);
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/Tc/q=0."+fmt.format(percent*1000)+".txt";
 			PrintUtil.printlnToFile(SaveAs, t, Cv);
 		}
 	}
 	
	
 	public void findHs(IsingStructure ising,double Hmin, double Hmax, double dH, String dynamics)
 	{
 		double X=0;
 		for(double h=Hmax; h>Hmin; h-=dH)
 		{
			params.set("H", h);
 			X=Susceptibility(ising,T,h,200,2000,10, dynamics);
 			PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/Hs/usefulrunsq=0."+fmt.format(percent*1000)+".txt","H=    ",(double)h );
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/Hs/q=0."+fmt.format(percent*1000)+".txt";
 			PrintUtil.printlnToFile(SaveAs, h, X, varianceX, usefulruns);
 		} 
 	}
 	
 	public void scanHs(IsingStructure ising, double min, double max, double dh, String dynamics)    //for the 1st search
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
		
		params.set("magnetization", Istemp.magnetization);
	}

	public void clear()
	{
		grid1.clear();
		grid2.clear();
	}
	
	public static void main (String[] Criticalpoint){
		new Control(new Criticalpoint(), "Kang Liu's critical points and spinodal" );
	}

	public void load(Control Criticalpoint){
		Criticalpoint.frame (grid1);
		Criticalpoint.frame (grid2);

		params.add("L1", 200);
		params.add("L2", 200);
		params.add("R", 10);
		params.add("NJ",-4.0);	
		params.add("percent", new DoubleValue(0.0,0,1).withSlider());
		params.add("biaspercent", new DoubleValue(0.0,0,1).withSlider());
		params.add("deadsites");	
		params.add("Dseed",1);
		params.add("Bseed",1);
		params.add("Sseed",1);
		
		params.addm("T", new DoubleValue(1.778, 0, 10).withSlider());
		params.addm("H", new DoubleValue(0, -2, 2).withSlider());
		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));
		
		params.add("MCS");
		params.add("copies");
		params.add("magnetization");
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
	    
	    IS.Dinitialization(Dseed, Bseed, 10, 10);
	    params.set("deadsites",IS.deadsites);
	    IS.Sinitialization(0, Sseed);
	    
	    Job.animate();
	    //findTc(IS,3.52,3.59,0.001, dynamics);
	    findTcviaX(IS,3.90,4.10,0.005, dynamics);
        
	    //T=params.fget("T");
	    //scanHs(IS,0,1.27,0.01, dynamics);
	    //startH=1.010;
	    //findHs(IS,startH-0.5,startH,0.005, dynamics);
      
		//CriticalpointsCv(IS, 4.00, 3.70, 0.01, 4, 2000, 2000, 5, dynamics);
	    
	    
	    
	    
	    
	    

	}
	



}