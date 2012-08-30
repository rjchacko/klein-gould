package kang.ising.BasicStructure;


import scikit.jobs.Job;
import chris.util.Random;

public class FCIsing {

	public int N;
	public double Nu, Nd;
	public double M;
	public double m;
	public double J;
	public double percent;   //the percent of the dilution
	public double livesites;  //the livesites=N*(1-percent)
	
	public FCIsing(int N)
	{
		this.N=N;
		this.livesites=N;
	}
	
	
	public FCIsing clone()
	{
		FCIsing copy=new FCIsing(N);
		copy.Nu=Nu;
		copy.Nd=Nd;
		copy.M=M;
		copy.m=m;
		copy.J=J;
		copy.percent=percent;
		copy.livesites=livesites;
		
		return copy;
	}
	
	public void dilute(double percent)
	{
		this.percent=percent;
		this.livesites=(int)(this.N*(1-percent));
	}
	
	public void setJ(double NJ)
	{

		this.J=NJ/(this.N-1);                     //normalize the interaction constant
	}
	
	
	
	public void initializeDilution(double m)    //initialization of the dilute system
	{
		this.Nu=Math.round((1-percent+m)*N/2);
		this.Nd=this.livesites-this.Nu;
		this.M=this.Nu-this.Nd;
	    this.m=M/N;
	    assert (this.m==m);
	    
	   
	}
	
	
	public void initialize(double m)
	{
		Nu=Math.round((1+m)*N/2);
		Nd=N-Nu;
		M=Nu-Nd;
	    this.m=M/N;
	    assert (this.m==m);
	}
	
	public void initialize(int delta)
	{
		Nu=N/2+delta/2;
		Nd=N/2-delta/2;
		M=delta;
		this.m=M/N;
	}
	
	public void spinflip(Random spinrand, Random fliprand, double T, double H)   //spinrand for choosing the up or down spin, fliprand for determining the spinflip
	{
		double Echange=0;
		if(spinrand.nextDouble()<(Nu/livesites))
		{
			Echange=-2*J*(Nu-Nd)+2*J+2*H;
			if(Echange<=0)
			{
				Nu--;
				Nd++;
				M=Nu-Nd;
			}
			else
			{
				if(fliprand.nextDouble()<=Math.exp(-Echange/T))
				{
					Nu--;
					Nd++;
					M=Nu-Nd;
				}
			}
		}
		else
		{
			Echange=2*J*(Nu-Nd)+2*J-2*H;
			if(Echange<=0)
			{
				Nu++;
				Nd--;
				M=Nu-Nd;
			}
			else
			{
				if(fliprand.nextDouble()<=Math.exp(-Echange/T))
				{
					Nu++;
					Nd--;
					M=Nu-Nd;
				}
			}
			
		}
		m=M/N;
	}
	
	
	public void spinflipGlauber(Random spinrand, Random fliprand, double T, double H)   //Glauber dynamics
	{
		double Echange=0;
		double Pglauber=0;
		if(spinrand.nextDouble()<(Nu/livesites))
		{
			Echange=-2*J*(Nu-Nd)+2*J+2*H;
			Pglauber=1/(1+Math.exp(Echange/T));

				if(fliprand.nextDouble()<Pglauber)
				{
					Nu--;
					Nd++;
					M=Nu-Nd;
				}

		}
		else
		{
			Echange=2*J*(Nu-Nd)+2*J-2*H;
			Pglauber=1/(1+Math.exp(Echange/T));
			
				if(fliprand.nextDouble()<Pglauber)
				{
					Nu++;
					Nd--;
					M=Nu-Nd;
				}
			
			
		}
		m=M/N;
	}
	
	
	public void MCS(String dynamics, Random spinrand, Random fliprand, double T, double H, double ratio)
	{
		int flipnumber=(int)(ratio*livesites);
		for(int j=0; j<flipnumber; j++)
		{
			if(dynamics=="Metropolis")
				spinflip(spinrand, fliprand, T, H);
			if(dynamics=="Glauber")
				spinflipGlauber(spinrand, fliprand, T, H);
			//Job.animate();
		}
		//Job.animate();
	}
	
	
	public void MCSnoise(String dynamics, Random spinrand, Random fliprand, Random Trand, double T, double dT, double H, double ratio)
	{
		int flipnumber=(int)(ratio*N);
		double t=T;
		for(int j=0; j<flipnumber; j++)
		{
	
 	 	      if(Trand.nextDouble()<0.5)
 	 	     		t=T-dT/2;
 	    	  else
 	 			    t=T+dT/2;

			 if(dynamics=="Metropolis")
				  spinflip(spinrand, fliprand, t, H);
			 if(dynamics=="Glauber")
				  spinflipGlauber(spinrand, fliprand, t, H);
		}
	}
	
	public double Energy(double H)
	{
		double E=0;
		E=J/2*(M*M)-J/2*livesites-H*M;
		return E;
	}
	
    public double Fluctuation(double data[], int size)
    {
    	double sum=0;
    	for(int j=0; j<size; j++)
    	{
    		sum+=data[j];
    	}
    	double avg= sum/size;
    	double sumD=0;
    	for(int k=0;k<size;k++)
    	{
    		sumD+=(data[k]-avg)*(data[k]-avg);
    	}
    	return sumD/size;
    }

	
	
}
