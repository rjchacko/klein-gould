package kang.ising.MutiplicativeNoise; // noise-3 is the multiplicative noise with the Gaussian distribution, this is the right noise to map ising model onto black-scholes equation

import java.text.DecimalFormat;

import kang.ising.Randomtemperature;
import kang.ising.BasicStructure.FCIsing;
import kang.ising.BasicStructure.MeanFieldIsingStructure;

import java.awt.Color;
import java.text.DecimalFormat;

import chris.util.PrintUtil;
import chris.util.Random;

import scikit.graphics.ColorPalette;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.Control;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;


public class FCIsingRandomT extends Simulation
{
	public int N,L;
	public double NJ;
	public double T,H;
	public double threshold;
	public int copies;
	public int runs;
	public double m, E;
	
	public int duration[];
	public double MeanDuration;
	public double SDDuration;
	
	
	public FCIsing FCising;
	public FCIsing Istemp;
	
	
	private DecimalFormat fmt = new DecimalFormat("000");
	private DecimalFormat smt = new DecimalFormat("000000");
	
	public double[] Gaussian(Random rand) // the function using Box-Muller transformation to generate a gaussian distributed random number
	{
		double x1, x2, w;
		double y[]=new double[2];
		
		do{
			x1=2.0*rand.nextDouble()-1.0;
			x2=2.0*rand.nextDouble()-1.0;
			w=x1*x1+x2*x2;
		}while(w>=1.0);
		
		w=Math.sqrt((-2.0)*Math.log(w)/w);
		y[0]=x1*w;
		y[1]=x2*w;
		return y;
		
	}
	
 	public void Heatup(String dynamics, FCIsing ising, int hseed, int steplimit)
 	{
	    Random heatrand=new Random(hseed);
 	    for(int heat=0; heat<steplimit; heat++)
	    {
	    	ising.MCS(dynamics, heatrand,heatrand,999,0,1);
	    	Job.animate();
			params.set("MCS", heat-steplimit);
	    }
 	}
 	
 	public void SingleRun(String dynamics, FCIsing ising, int Tseed, int Fseed, double T, double dT, int stepstart, int steplimit, int noisetype)
 	{
 		params.set("#run", 0);
 		Random Trand=new Random(Tseed);
 		Random Frand=new Random(Fseed);
 		String SaveAs = "/Users/liukang2002507/Desktop/simulation/FCRandomT/"+dynamics+"/noise"+fmt.format(noisetype)+"/L="+fmt.format(L)+" run#"+fmt.format(Tseed+1)+"<"+fmt.format(T*100)+","+fmt.format(dT*100)+">"+".txt";
 		
 		
    	PrintUtil.printlnToFile(SaveAs ,"InitialM=    ",(double)ising.m);
		
 		double t=T;
 		for(int step=0; step<steplimit; step++)
 		{
 			if(noisetype==1)
 			{
 			    t=T-dT/2+dT/2*(Trand.nextDouble());
 	 			ising.MCS(dynamics,Frand, Frand, t, 0, 1);
 			}
 			if(noisetype==2)
 			{
 				if(Trand.nextDouble()<0.5)
 					t=T-dT/2;
 				else
 					t=T+dT/2;
 	 			ising.MCS(dynamics,Frand, Frand, t, 0, 1);
 			}
	 		if(noisetype==3)
 	 		{
 	 			Istemp.MCSnoise(dynamics, Frand, Frand, Trand, T, dT, 0, 1);
 	 		}
	 		if(noisetype==0)
 	 		{
 	 			t=T;
 	 	 		Istemp.MCS(dynamics, Frand, Frand, t, 0, 1);
 	 		}
	 		
 			params.set("MCS", stepstart+step);
 			params.set("T", t);
 			m=ising.M/N;
 			Job.animate();
 			PrintUtil.printlnToFile(SaveAs, step+stepstart ,(double)m, (double)t);
 			
 		}
 	}
 	
 	public int Duration(double mag[], double threshold, int size, int avglength)  // given an array of magnetization, first find out the saturation magnetzation, and determine the time during which the system is stay near zero
 	{
 		int D=9999;
 		double totalM=0;
 		for(int i=0; i< avglength; i++)
 		{
 			totalM+=mag[size-avglength+i];
 		}
 		
 		double Ms= totalM/avglength;
 		double criticalM= Ms*threshold;
 		boolean keepon=true;
 		for(int j=0; (j< size)&(keepon); j++)
 		{
 			if(mag[j]*mag[j]>criticalM*criticalM)
 				{
 				keepon=false;
 				D=j;
 				}
 		}
 		return D;
 	}
 	
 	public double Mean(int data[], int size)
 	{
 		double total=0;
 		double mean=0;
 		for(int q=0; q<size; q++)
 		{
 			total+=data[q];
 		}
 		mean=total/size;
 		return mean;
 	}
 	
 	public double SD(int data[], int size, double mean)
 	{
 		double totalSD=0;
 		double SD=0;
 		for (int p=0; p<size; p++)
 		{
 			totalSD+=((data[p]-mean)*(data[p]-mean));
 		}
 		SD=Math.sqrt(totalSD/size);
 		
 		return SD;
 	}
 	
 	public double Mean(double data[], int size)
 	{
 		double total=0;
 		double mean=0;
 		for(int q=0; q<size; q++)
 		{
 			total+=data[q];
 		}
 		mean=total/size;
 		return mean;
 	}
 	
 	public double SD(double data[], int size, double mean)
 	{
 		double totalSD=0;
 		double SD=0;
 		for (int p=0; p<size; p++)
 		{
 			totalSD+=((data[p]-mean)*(data[p]-mean));
 		}
 		SD=Math.sqrt(totalSD/size);
 		
 		return SD;
 	}
 	
 	public double SpecificHeat(FCIsing ising, double T, int presteplimit,int steplimit, String dynamics)
 	{
 		double energy[]= new double[steplimit];
 		double Cv=0;
 		double MeanE=0;
 		//double mag[]= new double[steplimit];
 		Heatup(dynamics, ising, 1, 10);
 		Random Cvrand= new Random(999);
 		
 		
 	    for(int prestep=0; prestep<presteplimit; prestep++)
 	    {
 	    	ising.MCS(dynamics, Cvrand, Cvrand, T, 0, 1);
 	    	Job.animate();
 	    	params.set("MCS", prestep-presteplimit);
 	    }
 	    
 	    for(int step=0; step<steplimit; step++)
 	    {
 	    	ising.MCS(dynamics, Cvrand, Cvrand, T, 0, 1);
 	    	Job.animate();
 	    	params.set("MCS", step);
 	    	energy[step]=ising.Energy(0);
 	    }
 		MeanE=Mean(energy,steplimit);
 		Cv=SD(energy,steplimit,MeanE);
 		
 		return Cv;
 	}
 	
 	public void SearchforTc(int L, double MinT, double MaxT, double dT, int start, int limit, String dynamics)
 	{
    	FCising=new FCIsing(L*L);
        FCising.setJ(-4.0);
        FCising.initialize(0.0);	    
        Istemp=new FCIsing(L*L);
        Istemp=FCising.clone();
        Job.animate();
        String path="/Users/liukang2002507/Desktop/simulation/FCRandomT/"+dynamics+"/Tc/L="+smt.format(L)+".txt";
        double specificheat=0;
        params.set("L", L);
        
        for(double t=MinT; t<MaxT; t+=dT)
        {
        	params.set("T", t);
        	specificheat=SpecificHeat(Istemp,t, start, limit, dynamics);
            PrintUtil.printlnToFile(path , t , specificheat);
        }
        
        
 	}
 	
 	public void SizeEffectOnDuration(int copies, int Min, int Max, double T, double dT, int steplimit, int noisetype, String dynamics)
 	{
 		String path = "/Users/liukang2002507/Desktop/simulation/FCRandomT/"+dynamics+"/noise"+fmt.format(noisetype)+"/SizeEffect"+"<"+fmt.format(T*100)+","+fmt.format(dT*100)+">"+".txt";
 		for(int i=Min; i<=Max; i=i*2)
 		{
 			params.set("L",i);
 			MultiRuns(i, copies, T, dT, steplimit, noisetype, dynamics);
 			PrintUtil.printlnToFile(path , i , MeanDuration, SDDuration);
 		}
 	}
 	
 	public void SEwithoutANoise(int copies, int runs, int Min, int Max, double T, double dT, int steplimit, int noisetype, String dynamics)//size effect on duration without the additive noise
 	{
 		String path = "/Users/liukang2002507/Desktop/simulation/FCRandomT/NoAdditiveNoise/"+dynamics+"/noise"+fmt.format(noisetype)+"/SizeEffect"+"<"+fmt.format(T*100)+","+fmt.format(dT*100)+">"+".txt";
 		for(int i=Min; i<=Max; i=i*2)
 		{
 			params.set("L",i);
 			MRNoANoise(i, copies, runs, T, dT, steplimit, noisetype, dynamics);
 			PrintUtil.printlnToFile(path , i , MeanDuration, SDDuration);
 		}
 	}
 	
 	
 	public void AvgMagnetization(FCIsing ising, double dataM[], double noiseT[], double T, double dT, int L, int noiseseed, int runs, int steplimit, int noisetype, String dynamics)
 	{
 		double[] tempm=new double[steplimit];
 		double[] totalm=new double[steplimit];
 		for(int ii=0; ii<steplimit; ii++)
 		{
 			tempm[ii]=0;
 			totalm[ii]=0;
 		}
 		
 		for(int r=0; r<runs; r++)
 		{
 			Random spinfliprand= new Random(r+1);
 			params.set("#run", r+1);
 			
 			String SaveAs = "/Users/liukang2002507/Desktop/simulation/FCRandomT/NoAdditiveNoise/"+dynamics+"/noise"+fmt.format(noisetype)+"/MultiRuns/L="+fmt.format(L)+"<"+fmt.format(T*100)+","+fmt.format(dT*100)+">"+"noise#"+fmt.format(noiseseed)+" run#"+fmt.format(r+1)+".txt";
 			
 			Istemp=ising.clone();
 			PrintUtil.printlnToFile(SaveAs ,"InitialM=    ",(double)Istemp.m);
 			
 			for(int step=0; step<steplimit; step++)
 			{
 				Istemp.MCS(dynamics, spinfliprand, spinfliprand, noiseT[step], 0, 1);
 				m=Istemp.M/N;
 				tempm[step]=m;
 				Job.animate();
 				params.set("MCS", step);
 	 			params.set("T", noiseT[step]);
 				PrintUtil.printlnToFile(SaveAs, step,(double)m, (double)noiseT[step]);
 			}
 			
 			for(int ff=0; ff<steplimit; ff++)
 			{
 				if(tempm[steplimit-1]<0)
 					totalm[ff]-=tempm[ff];
 				else
 					totalm[ff]+=tempm[ff];
 			}
 			
 		}
 		for(int f=0; f<steplimit; f++)
 		{
 			dataM[f]=totalm[f]/runs;
 		}
 		
 		
 	}
 	
    public void MRNoANoise(int L, int copies, int runs, double T, double dT, int steplimit, int noisetype, String dynamics)
 	{
 		double[]tempM= new double[steplimit];
 		double[]noiseT= new double[steplimit];
 		
	    FCising=new FCIsing(L*L);
	    FCising.setJ(NJ);
	    Istemp= new FCIsing(L*L);
	    Istemp.initialize(0.0);
	    
	    
	    for(int n=0; n<copies; n++)
	    {
	    	FCising.initialize(0.0);
	    	Heatup(dynamics,FCising, n+1, 100);      // heat up the initial configuration
	    	params.set("#noise", n+1);
	    	Random Trand=new Random(n+1);
	    	
	    	String Savepath = "/Users/liukang2002507/Desktop/simulation/FCRandomT/NoAdditiveNoise/"+dynamics+"/noise"+fmt.format(noisetype)+"/MultiRuns/AverageM/L="+fmt.format(L)+"<"+fmt.format(T*100)+","+fmt.format(dT*100)+">"+"noise#"+fmt.format(n+1)+".txt";
	    	
	    	for(int s=0; s<steplimit;s++)   // this cycle is used to generate the multipicative noise in temperature
	    	{
 				if(noisetype==1)
 	 			{
 					noiseT[s]=T-dT/2+dT/2*(Trand.nextDouble());
 	 			}
 	 			if(noisetype==2)
 	 			{
 	 				if(Trand.nextDouble()<0.5)
 	 					noiseT[s]=T-dT/2;
 	 				else
 	 					noiseT[s]=T+dT/2;
 	 			}
 	 			if(noisetype==0)
 	 			{
 	 				noiseT[s]=T;
 	 			}
 	 			
 	 			if(noisetype==3)
 	 			{
 	 				noiseT[s]=T+dT*Gaussian(Trand)[0];
 	 				s++;
 	 				noiseT[s]=T+dT*Gaussian(Trand)[1];
 	 			}
	    	}
	    	
	    	AvgMagnetization(FCising, tempM, noiseT, T, dT, L ,n+1, runs, steplimit, noisetype, dynamics);
	    	for(int sss=0; sss<steplimit; sss++)
	    	{
	    		PrintUtil.printlnToFile(Savepath, sss, (double)tempM[sss], (double)noiseT[sss]);
	    	}
	    	duration[n]=Duration(tempM, threshold, steplimit, 100);
	 
	    }
 		MeanDuration=Mean(duration,copies);
 		SDDuration=SD(duration,copies,MeanDuration);    
	    
	    
 	}
 	
 	public void MultiRuns(int L, int copies, double T, double dT, int steplimit, int noisetype, String dynamics)
 	{
 		double[]tempM= new double[steplimit];
 
	    FCising=new FCIsing(L*L);
	    FCising.setJ(NJ);
	    FCising.initialize(0.0);
	    Istemp=new FCIsing(L*L);
	    Heatup(dynamics,FCising, 1, 100);
 		
 		for(int r=0; r< copies; r++)
 		{
 	 		params.set("#run", r+1);
 	 		Random Trand=new Random(r+1);
 	 		Random Frand=new Random(r+1);
 			Istemp=FCising.clone();
 			
 			
 	 		String SaveAs = "/Users/liukang2002507/Desktop/simulation/FCRandomT/"+dynamics+"/noise"+fmt.format(noisetype)+"/MultiRuns/L="+fmt.format(L)+" run#"+fmt.format(r+1)+"<"+fmt.format(T*100)+","+fmt.format(dT*100)+">"+".txt";
 
 	    	PrintUtil.printlnToFile(SaveAs ,"InitialM=    ",(double)Istemp.m);
 			
 			double t=T;
 			for(int step=0;step<steplimit; step++)
 			{
 				if(noisetype==1)
 	 			{
 					t=T-dT/2+dT/2*(Trand.nextDouble());
 	 	 			Istemp.MCS(dynamics, Frand, Frand, t, 0, 1);
 	 			}
 	 			if(noisetype==2)
 	 			{
 	 				if(Trand.nextDouble()<0.5)
 	 					t=T-dT/2;
 	 				else
 	 					t=T+dT/2;
 	 	 			Istemp.MCS(dynamics,Frand, Frand, t, 0, 1);
 	 			}
 	 			if(noisetype==0)
 	 			{
 	 				t=T;
 	 	 			Istemp.MCS(dynamics, Frand, Frand, t, 0, 1);
 	 			}

 	 			if(noisetype==3)
 	 			{
 	 				Istemp.MCSnoise(dynamics, Frand, Frand, Trand, T, dT, 0, 1);
 	 			}
 	 			
 	 			params.set("MCS", step);
 	 			params.set("T", t);
 	 			m=Istemp.M/N;
 	 			PrintUtil.printlnToFile(SaveAs, step,(double)m, (double)t);
 	 			tempM[step]=m;
 	 			Job.animate();	
 			}
 			duration[r]=Duration(tempM, threshold, steplimit, 100);
 		}
 		MeanDuration=Mean(duration,copies);
 		SDDuration=SD(duration,copies,MeanDuration);
 		
 	}
 	
	public void animate()
	{
		params.set("magnetization", Istemp.m);	
	}
	
	public void clear()
	{
		
	}
 	
	public static void main (String[] FCIsingRandomT){
		new Control(new FCIsingRandomT(), "Kang Liu's fully connected Ising model with multiplicative noises" );
	}

	public void load(Control FCIsingRandomT){


		params.add("L", 800);
		params.add("NJ",-4.0);
		params.add("copies" , 20);
		params.add("runs",20);
		
		params.addm("T", new DoubleValue(3, 0, 5).withSlider());
		params.addm("H", new DoubleValue(0, -2, 2).withSlider());
		params.addm("Dynamics", new ChoiceValue("Metropolis", "Glauber"));
		
		params.add("#run");
		params.add("#noise");
		params.add("MCS");

		params.add("magnetization");
		//params.add("energy");
	}
	
	public void run(){
		
		
		L = (int)params.fget("L");
		N=L*L;
		NJ = params.fget("NJ");
		copies= (int)params.fget("copies");
		runs=(int)params.fget("runs");
		threshold=0.5;
		duration= new int[copies];
		String dynamics= params.sget("Dynamics");
		
		//SizeEffectOnDuration(copies, 100, 1600, 1.78, 0, 1000, 0, dynamics);
		//SizeEffectOnDuration(copies, 100, 1600, 1.78, 0.2, 1000, 2, dynamics);
		//SizeEffectOnDuration(copies, 100, 1600, 1.78, 0.4, 1000, 2, dynamics);
		//SizeEffectOnDuration(copies, 100, 1600, 1.78, 0.6, 1000, 2, dynamics);
		
		SEwithoutANoise(copies, runs, 100, 400, 1.78, 0.2, 1000, 2, dynamics);
		SEwithoutANoise(copies, runs, 100, 400, 1.78, 0.4, 1000, 0, dynamics);
		SEwithoutANoise(copies, runs, 100, 400, 1.78, 0.0, 1000, 0, dynamics);
		
		//SearchforTc(100,3.88,3.96,0.002, 1000, 5000);  //[3.88,3.96]
		//SearchforTc(200,3.92,3.98,0.002, 1000, 5000);  //[3.92,3.98]
		//SearchforTc(400,3.94,4.00,0.002, 1000, 3000);  //[3.94,4.00]
		//SearchforTc(800,3.945,3.97,0.002, 1000, 2000);  //[3.94,4.00]-[3.945,3.97]
		//SearchforTc(1600,3.955,3.965,0.002, 1000, 1000); //[3.94,4.00]-[3.955,3.965]
		//SearchforTc(3200,3.96,3.98,0.005, 1000, 1000);  //[3.96,4.00]-[3.96,3.98]
		
	    {
	    	FCising=new FCIsing(N);
	        FCising.setJ(NJ);
	        FCising.initialize(0.0);	    
	        Istemp=new FCIsing(N); 
	        Job.animate();
	        Heatup(dynamics,FCising, 1, 100);
	        Job.animate();
	        Istemp=FCising.clone();
	    }
	    
	    
	    
	    //SingleRun(Istemp, 0, 0, 3.9, 0.2, 0, 1000, 2);
	    //SingleRun(Istemp, 0, 0, 3.9, 0.4, 500, 2000, 1);
	    
	    

    
	    

	}
 	

	
}