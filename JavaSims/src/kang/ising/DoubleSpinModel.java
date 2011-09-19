package kang.ising;

import java.text.DecimalFormat;


import kang.ising.BasicStructure.FCIsing;
import kang.ising.BasicStructure.FCIsingWithDamage;
import kang.ising.MutiplicativeNoise.FCIsingRandomT;


import java.text.DecimalFormat;

import chris.util.PrintUtil;
import chris.util.Random;


import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.Control;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

public class DoubleSpinModel extends Simulation
{
	public int N,L;
	public double NJ;
	public double q, a;    //q is the percentage of the 2nd kind of spins, a is the value of these spins
	public double T,H;
	public double Tc;

	public double m, E;
	

	
	public FCIsingWithDamage FCisingWD;
	public FCIsingWithDamage IstempWD;
	
	
	private DecimalFormat fmt = new DecimalFormat("000");
	private DecimalFormat smt = new DecimalFormat("000000");
	
 	public void Heatup(String dynamics, FCIsingWithDamage ising, int hseed, int steplimit)
 	{
	    Random heatrand=new Random(hseed);
 	    for(int heat=0; heat<steplimit; heat++)
	    {
	    	ising.MCS(dynamics, heatrand,heatrand,999,0,1);
	    	Job.animate();
			params.set("MCS", heat-steplimit);
	    }
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
 	
	public double SpecificHeat(FCIsingWithDamage ising, double T, double H, int presteplimit,int steplimit, String dynamics, int seed)
 	{
 		double energy[]= new double[steplimit];
 		double NCv=0;
 		double Cv=0;
 		double MeanE=0;
 		//double mag[]= new double[steplimit];
 		Heatup(dynamics, ising, seed, 100);
 		Random Cvrand= new Random(seed);
 		
 		
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
 	    	energy[step]=ising.Energy2(0);
 	    }
 		MeanE=Mean(energy,steplimit);
 		NCv=SD(energy,steplimit,MeanE);
 		Cv=NCv/ising.N;
 		
 		return Cv;
 	}
	
	public double Susceptibility(FCIsingWithDamage ising, double T, double H, int presteplimit,int steplimit, String dynamics, int seed)
 	{
 		double magnetization[]= new double[steplimit];
 		double NChi=0;
 		double Chi=0;
 		double MeanM=0;
 		//double mag[]= new double[steplimit];
 		Heatup(dynamics, ising, seed, 100);
 		Random Chirand= new Random(seed);
 		
 		
 	    for(int prestep=0; prestep<presteplimit; prestep++)
 	    {
 	    	ising.MCS(dynamics, Chirand, Chirand, T, 0, 1);
 	    	Job.animate();
 	    	params.set("MCS", prestep-presteplimit);
 	    }
 	    
 	    for(int step=0; step<steplimit; step++)
 	    {
 	    	ising.MCS(dynamics, Chirand, Chirand, T, 0, 1);
 	    	Job.animate();
 	    	params.set("MCS", step);
 	    	magnetization[step]=ising.m;
 	    }
 		MeanM=Mean(magnetization,steplimit);
 		NChi=SD(magnetization,steplimit,MeanM);
 		Chi=NChi/ising.N;
 		
 		return Chi;
 	}
	
	
 	public void SearchforTc(int L, double q, double a, double MinT, double MaxT, double dT, int start, int limit, int copies,String dynamics)
 	{
    	FCisingWD=new FCIsingWithDamage(L*L);
        FCisingWD.setJ(-4.0);
        FCisingWD.initialize(q,a);	    
        IstempWD=new FCIsingWithDamage(L*L);
        IstempWD=FCisingWD.clone();
        Job.animate();
        String path="/Users/liukang2002507/Desktop/simulation/DoubleSpinModel/Tc/"+dynamics+"/L="+smt.format(L)+" q="+fmt.format(q*100)+" a="+fmt.format(a*100)+".txt";
        String check="/Users/liukang2002507/Desktop/simulation/DoubleSpinModel/Tc/"+dynamics+"/check L="+smt.format(L)+" q="+fmt.format(q*100)+" a="+fmt.format(a*100)+".txt";
        double specificheat=0;
        double deltaCv=0;
        params.set("L", L);
        PrintUtil.printlnToFile(path , "Tc(q,a)=", 4*(1-q+q*a*a));
        
        for(double t=MinT; t<MaxT; t+=dT)
        {
        	params.set("T", t);
        	double[] Cvtemp= new double [copies];
        	for(int c=0; c<copies; c++)
        	{
        		Cvtemp[c]=SpecificHeat(IstempWD,t, 0, start, limit, dynamics,c+1);
        		PrintUtil.printlnToFile(check , t , c, Cvtemp[c]);
        		params.set("#run",c+1);
        	}
        	specificheat=Mean(Cvtemp,copies);
        	deltaCv=SD(Cvtemp,copies,specificheat);
            PrintUtil.printlnToFile(path , t , specificheat,deltaCv);
        }
        
        
 	}
 	
 	public void SearchforHs(int L, double q, double a, double MinH, double MaxH, double dH, int start, int limit, int copies,String dynamics)
 	{
    	FCisingWD=new FCIsingWithDamage(L*L);
        FCisingWD.setJ(-4.0);
        FCisingWD.initialize(q,a);	    
        IstempWD=new FCIsingWithDamage(L*L);
        IstempWD=FCisingWD.clone();
        Job.animate();
        String path="/Users/liukang2002507/Desktop/simulation/DoubleSpinModel/Hs/"+dynamics+"/L="+smt.format(L)+" q="+fmt.format(q*100)+" a="+fmt.format(a*100)+".txt";
        String check="/Users/liukang2002507/Desktop/simulation/DoubleSpinModel/Hs/"+dynamics+"/check L="+smt.format(L)+" q="+fmt.format(q*100)+" a="+fmt.format(a*100)+".txt";
        double susceptibility=0;
        double deltaChi=0;
        params.set("L", L);
        Tc=4*(1-q+q*a*a);
        PrintUtil.printlnToFile(path , "Tc(q,a)=", 4*(1-q+q*a*a));
        params.set("T", 4*Tc/9);
        for(double h=MinH; h<MaxH; h+=dH)
        {
        	params.set("H", h);
        	double[] Chitemp= new double [copies];
        	for(int c=0; c<copies; c++)
        	{
        		Chitemp[c]=Susceptibility(IstempWD , 4/9*Tc, h, start, limit, dynamics,c+1);
        		PrintUtil.printlnToFile(check , h , c, Chitemp[c]);
        		params.set("#run",c+1);
        	}
        	susceptibility=Mean(Chitemp,copies);
        	deltaChi=SD(Chitemp,copies,susceptibility);
            PrintUtil.printlnToFile(path , h , susceptibility,deltaChi);
        }
        
        
 	}
 	
 	
	public void animate()
	{
		params.set("magnetization", IstempWD.Magnetization());
		params.set("energy1", IstempWD.Energy1(H));
		params.set("energy2", IstempWD.Energy2(H));
		
		
		/*params.set("Ms", IstempWD.Ms);
		params.set("Mx", IstempWD.Mx);
		params.set("N", IstempWD.N);
		
		
		params.set("Nsu", IstempWD.Nsu);
		params.set("Nsd", IstempWD.Nsd);
		params.set("Nxu", IstempWD.Nxu);
		params.set("Nxd", IstempWD.Nxd);*/
		
	}
 	
	public void clear()
	{
		
	}
	
	public static void main (String[] DoubleSpinModel){
		new Control(new DoubleSpinModel(), "Kang Liu's fully connected Ising model with 2 kinds of spins" );
	}

	public void load(Control DoubleSpinModel){


		params.add("L", 800);
		params.add("NJ",-4.0);
	    
		params.add("a", new DoubleValue(0.5, 0, 1).withSlider());
		params.add("q", new DoubleValue(0.2,0,1).withSlider());
        params.addm("Dynamics", new ChoiceValue("Metropolis", "Glauber"));
		
		params.addm("T", new DoubleValue(3, 0, 50).withSlider());
		params.addm("H", new DoubleValue(0, -2, 2).withSlider());
		params.add("Tc");
		
		params.add("#run");
		params.add("MCS");

		params.add("magnetization");
		params.add("energy1");
		params.add("energy2");
		
		/*params.add("Ms");
		params.add("Mx");
		params.add("N/");
		
		params.add("Nsu");
		params.add("Nsd");
		params.add("Nxu");
		params.add("Nxd");*/
		
		
	}
 	
	public void run(){
		
		
		L = (int)params.fget("L");
		N=L*L;
		NJ = params.fget("NJ");
        String dynamics= params.sget("Dynamics");
        q=params.fget("q");
        a=params.fget("a");
        
        T=params.fget("T");
        H=params.fget("H");
        Tc=4*(1-q+q*a*a);
        params.set("Tc", Tc);
        Random spinrand=new Random(1);
        Random fliprand=new Random(2);
        
        
    	/*FCisingWD=new FCIsingWithDamage(N);
        FCisingWD.setJ(NJ);
        FCisingWD.initialize(q,a);	    
        IstempWD=new FCIsingWithDamage(N); 
        Job.animate();
        Heatup(dynamics,FCisingWD, 1, 100);
        Job.animate();
        IstempWD=FCisingWD.clone();
        
        for(int time=0; time<100000; time++)
        {
        	T=params.fget("T");
            H=params.fget("H");
        	IstempWD.MCS(dynamics, spinrand, fliprand, T, H, 1);
        	Job.animate();
        	params.set("MCS", time);
        }*/
        
        double MinT= Tc-0.2;
        double MaxT= Tc+0.2;
        double dT=0.02;
        
        
        SearchforTc(L,q,a, MinT, MaxT, dT, 1000, 2000, 10,dynamics);


		
		
	
	}
}