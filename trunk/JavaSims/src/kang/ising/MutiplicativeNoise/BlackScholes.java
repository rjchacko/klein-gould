package kang.ising.MutiplicativeNoise;

import java.text.DecimalFormat;

import kang.ising.BasicStructure.FCIsing;

import chris.util.PrintUtil;
import chris.util.Random;

import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.Control;
import scikit.jobs.params.ChoiceValue;

public class BlackScholes extends Simulation
{
	public int N,L;
	public double NJ;
	public double T,H;
	public double threshold;
	public int copies;
	public int runs;
	public double m, E;
	public double m0;    
	
	public int duration[];
	public double MeanDuration;
	public double SDDuration;
	
	
	public FCIsing FCising;
	public FCIsing Istemp;
	
	
	private DecimalFormat fmt = new DecimalFormat("000");
	//private DecimalFormat smt = new DecimalFormat("000000");
	
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

 	public void SEwithoutANoise(int copies, int runs, int Min, int Max, double T, double dT, int steplimit, int noisetype, String dynamics)//size effect on duration without the additive noise
 	{
 		String path = "/Users/liukang2002507/Desktop/simulation/BlackScholes/"+dynamics+"/noise"+fmt.format(noisetype)+"/SizeEffect"+"<"+fmt.format(T*100)+","+fmt.format(dT*100)+">"+".txt";
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
 			
 			String SaveAs = "/Users/liukang2002507/Desktop/simulation/BlackScholes/"+dynamics+"/noise"+fmt.format(noisetype)+"/noiseruns/copyruns/L="+fmt.format(L)+"<"+fmt.format(T*100)+","+fmt.format(dT*100)+">"+"noise#"+fmt.format(noiseseed)+" run#"+fmt.format(r+1)+".txt";
 			
 			Istemp=ising.clone();
 			PrintUtil.printlnToFile(SaveAs ,"InitialM=    ",(double)Istemp.m);
 			
 			for(int step=0; step<steplimit; step++)
 			{
 				Istemp.MCS(dynamics, spinfliprand, spinfliprand, noiseT[step], 0, 1);
 				m=Istemp.m;
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
	    	FCising.initialize(m0);
	    	params.set("#noise", n+1);
	    	Random Trand=new Random(n+1);
	    	
	    	String Savepath = "/Users/liukang2002507/Desktop/simulation/BlackScholes/"+dynamics+"/noise"+fmt.format(noisetype)+"/noiseruns/L="+fmt.format(L)+"<"+fmt.format(T*100)+","+fmt.format(dT*100)+">"+"noise#"+fmt.format(n+1)+".txt";
	    	
	    	for(int s=0; s<steplimit;s++)   // this cycle is used to generate the multipicative noise in temperature
	    	{
 	 			if(noisetype==0)
 	 			{
 	 				noiseT[s]=T;
 	 			}
 	 			
 	 			if(noisetype==1)
 	 			{
 	 				noiseT[s]=T+dT*(Gaussian(Trand)[0]);
 	 				s++;
 	 				noiseT[s]=T+dT*(Gaussian(Trand)[1]);
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
    
	public void animate()
	{
		params.set("magnetization", Istemp.m);	
	}

	public void clear()
	{
		
	}
 	
	public static void main (String[] BlackScholes){
		new Control(new BlackScholes(), "Kang Liu's Ising model for B-S equation" );
	}

	public void load(Control BlackScholes){


		params.add("L", 100);
		params.add("NJ",-4.0);
		params.add("copies" , 20);  //number of noises
		params.add("runs",20);   //number of lattice copies for ensemble average to get rid of the additive noise
		
		params.addm("T",3.0);
		params.addm("H",0.0);
		params.addm("m0", 0.001);
		params.addm("Dynamics", new ChoiceValue("Glauber","Metropolis"));
		
		params.add("#run");
		params.add("#noise");
		params.add("MCS");

		params.add("magnetization");
		//params.add("energy");
	}
	
	public void run(){
		
		
		NJ = params.fget("NJ");
		copies= (int)params.fget("copies");
		runs=(int)params.fget("runs");
		m0=params.fget("m0");
		
		threshold=0.2;
		duration= new int[copies];
		String dynamics= params.sget("Dynamics");
		
		SEwithoutANoise(copies, runs, 400, 400, 3.60, 0.40, 500, 1, dynamics);
		SEwithoutANoise(copies, runs, 400, 400, 3.60, 0.40, 500, 0, dynamics);
	}
	
	
}



