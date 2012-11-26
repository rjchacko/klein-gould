package kang.ising;

import java.awt.Color;
import java.text.DecimalFormat;

import chris.util.PrintUtil;
import chris.util.Random;



import kang.ising.BasicStructure.IsingStructure3D;
import kang.ising.BasicStructure.BasicTools;

import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

public class Criticalpoint3D extends Simulation
{
	Grid gridx=new Grid("gridx");
	Grid gridy=new Grid("gridy");
	Grid gridz=new Grid("gridz");
	
	
	public int L,M,R,deadsites,Dseed,Bseed,Sseed;
	public double percent,biaspercent,NJ;
	public double T,H;
	public IsingStructure3D IS;
	public IsingStructure3D Istemp;
	public BasicTools Tools;
	
	public int isingx[];
	public int isingy[];
	public int isingz[];
	
	
	public int progress;
    public int usefulruns=0;
    
    public double varianceX;
    public double varianceC;
    public double Cvdata[];
    
	private DecimalFormat fmt = new DecimalFormat("000");
	public double startH=0;
	
	
	public double SpecificHeat(IsingStructure3D ising, double T, double H, int presteplimit, int number, int copies, String dynamics)  //calculate the specific heat of a given system at T,
	{
	    double tempE[];
	    tempE= new double [number];
	    Cvdata= new double [copies];
	    double meanC=0;
	    
	    double totalM=0;
	    double averageM=0;
	    String check ="/Users/liukang2002507/Desktop/simulation/CriticalpointsCv3D/"+dynamics+"/check L="+fmt.format(ising.L1)+" R="+fmt.format(ising.R)+" q="+fmt.format(ising.percent*100)+".txt";
	    
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
	
	
	public void CriticalpointsCv(IsingStructure3D ising, double Tmax, double Tmin, double increment, double targetT, int limit, int number, int copies, String dynamics)
	{
		String path="/Users/liukang2002507/Desktop/simulation/CriticalpointsCv3D/"+dynamics+"/L="+fmt.format(ising.L1)+" R="+fmt.format(ising.R)+" q="+fmt.format(ising.percent*100)+".txt";
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
	
	
	public double SpecificHeatHs(IsingStructure3D ising, double T, double H, int presteplimit, int number, int copies, String dynamics)  //calculate the specific heat of a given system at T,
	{
	    double tempE[];
	    double usedE[];
	    tempE= new double [number];
	    usedE= new double [number-200-200];
	    double totalCv=0;
        usefulruns=0;
        double tempCv[];
        tempCv= new double[copies];
	    
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
				
				tempE[step]=Istemp.TotalEnergy(field);
				if(Istemp.TotalSpin()<0)
					{
					endofstep=1;
					lifetime=step;
					}
				Istemp.MCS(T, field, cflip, 1, dynamics);
				Job.animate();
				params.set("MCS", step);
				params.set("magnetization", Istemp.TotalSpin()/M);
			}
			
			if(lifetime>2000)
			{
				for(int t=0;t<number-200-200;t++)
				{
					usedE[t]=tempE[t+200];
				}
				tempCv[usefulruns]=IS.Fluctuation(usedE, number-200-200);
				totalCv+=tempCv[usefulruns];
				PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/Hs3D/Cv3D usefulrunsq=0."+fmt.format(percent*1000)+".txt",H, c);
				usefulruns++;
			}
		}
		double averageCv=totalCv/usefulruns;
		double totalCv2=0;
		for(int u=0;u<usefulruns; u++)
		{
			totalCv2+=(tempCv[u]-averageCv)*(tempCv[u]-averageCv);
		}
		varianceC=Math.sqrt(totalCv2/usefulruns);
		
		return averageCv;
		
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
		
		for(int i=0;i<L; i++)
			for(int j=0; j<L; j++)
			{
				isingx[i*L+j]=Istemp.spin[L/2*L*L+i*L+j];
				isingy[i*L+j]=Istemp.spin[i*L*L+L/2*L+j];
				isingz[i*L+j]=Istemp.spin[i*L*L+j*L+L/2];
			}
		
		gridx.setColors(ising);
		gridx.registerData(L, L, isingx);
		gridy.setColors(ising);
		gridy.registerData(L, L, isingy);
		gridz.setColors(ising);
		gridz.registerData(L, L, isingz);
		
		params.set("magnetization", Istemp.magnetization);
		params.set("intenergy", Istemp.totalintenergy);
	}

	public void clear()
	{
		gridx.clear();
		gridy.clear();
		gridz.clear();
	}
	
	public static void main (String[] Criticalpoint3D){
		new Control(new Criticalpoint3D(), "Kang Liu's critical points and spinodal in 3D" );
	}

	public void load(Control Criticalpoint3D){


		Criticalpoint3D.frameTogether("Display", gridx, gridy, gridz);
		
		
		params.add("L", 30);

		params.add("R", 3);
		params.add("NJ",-6.0);	
		params.add("percent", 0.00);
		params.add("biaspercent", 0.00);
		params.add("deadsites");	
		params.add("Dseed",2);
		params.add("Bseed",2);
		params.add("Sseed",1);
		
		params.addm("T", 6.0);
		params.addm("H", 0.0);
		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));
		
		params.add("MCS");
		params.add("copies");
		params.add("magnetization");
		params.add("intenergy");
	}
	
	
	public void run(){
		
		
		L = (int)params.fget("L");
		
		M = L * L * L;
		
		R = (int)params.fget("R");
		NJ = params.fget("NJ");
		isingx= new int[L*L];
		isingy= new int[L*L];
		isingz= new int[L*L];

		percent=params.fget("percent");
		biaspercent=params.fget("biaspercent");
		String dynamics= params.sget("Dynamics");
		
		Dseed = (int)params.fget("Dseed");
		Bseed = (int)params.fget("Bseed");
		Sseed = (int)params.fget("Sseed");
		
	    IS=new IsingStructure3D(L,L,L,R,NJ,percent,biaspercent,"square");
	    Istemp=new IsingStructure3D(L,L,L,R,NJ,percent,biaspercent,"square");
	    Tools=new BasicTools();
	    
	    IS.Dinitialization(Dseed, Bseed, 10, 10, 10);
	    params.set("deadsites",IS.deadsites);
	    IS.Sinitialization(0, Sseed);
	    
	    Job.animate();

	    
	   
		CriticalpointsCv(IS, 6.30, 6.00, 0.02, 6, 1000, 1000, 10, dynamics);
	    
	
	    
	    

	}
	
	
	
}
	
	
