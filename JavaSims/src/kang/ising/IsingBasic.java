package kang.ising;


import kang.ising.IsingStructure;
import chris.util.Random;
import chris.util.PrintUtil;


public class IsingBasic{
	

	
	public IsingStructure ising;

	public double T;     //temperature
	public double QuenchT;  //the temperature after the quench
	public double QuenchH;   //the temperature after the quench
	public double H;     //field
	public double field; //field for display purpose
	public double temperature; //temperature for display purpose

	

 	

 	

 	


 	
 	
	public void MCS(int range, double constJ, int spin[], Random flip, double ratio, double temp, double field)
	{
	    int j=0;
	    int rationalMCS= (int) (ratio*M);
	    for (int f=0; f< rationalMCS; f++)
	    {
		   j=(int) (flip.nextDouble()*M);
		   spinflip(range, constJ, spin, j, flip, temp, field);
		 
	    }
	    magnetization=magnetization();
	    totalenergy=totalintenergy+field*totalspin;
	}
	
	public void temperatureset(double finaltemp)
	{
		temperature=finaltemp;
	}
	
	public void fieldset(double finalfield)
	{
		field=finalfield;
	}
	
	public void fieldflip()
	{
	    fieldset(-field);
	}
	
	public void Dinitialization(int spin[], double p, double biasp, int a, int b, int Dseed, int Bseed)//dilution initialization
	{
		Random Drand= new Random(Dseed);
		Random Brand= new Random(Bseed);
		deadsites=0;
		int biasrangelabel[]=new int[M];
		
		for (int t=0;t<M;t++)
			biasrangelabel[t]=0;
			
		
		int cx, cy; //x,y index for the center
		cx=(int)L1/2;
		cy=(int)L2/2;
		if(biasp!=p)
		    rectangle(biasrangelabel, a,b, cx,cy);
		
		for(int j=0; j<M; j++)
		{
			if (biasrangelabel[j]==1)
				if(Brand.nextDouble()<biasp)
					{
					spin[j]=0;
					deadsites++;
					}
			if (biasrangelabel[j]==0)
				if(Drand.nextDouble()<p)
					{
					spin[j]=0;
					deadsites++;
					}
		}
	
	}
	
	public void Sinitialization(int spin[], int type, int Sseed)//spin config initialization (type: 0-random 1-spin up   2-spin down)
	{
		Random spinrand= new Random(Sseed);
		totalintenergy=0;
		totalspin=0;
		
		
		if(type==0)
			for (i=0; i<M; i++)
		    {
			   if(spin[i]!=0)
			   {
				   spin[i]=-1;
			   if (spinrand.nextDouble()> 0.5)
				   spin[i]=1;
			   }
				
		    }
		if(type==1)
			for (i=0; i<M; i++)
		    {
			   if(spin[i]!=0)
				   spin[i]=1;
		    }
		
		if(type==2)
			for (i=0; i<M; i++)
		    {
			   if(spin[i]!=0)
				   spin[i]=-1;
		    }
		
		for (i=0; i<M; i++)
		{

			initialcopy[i]=spin[i];  // here, make the copy of the system
			totalspin+=spin[i];      // calculate the total spin at the beginning
		}
		
		totalintenergy= totalIntEnergy(R, NJ, spin);
		
	}
	

	

	
	public void DoubleEnergyMetric(int spin[], int Dseed)  //calculate the energy metric for a given dilution realization by running the MCS on two systems 
	{
		
	}
	
	public void DoubleSpinMetric(int spin[], int Dseed)  //calculate the double spin metric
	{
		
	}
	
	public void SpecificHeat(int range, double constJ, int spin[], double temp, double field, int presteplimit, int totalstep)  //calculate the specific heat of a given system at T,
	{
		int s[]= new int [M];
		for(int j=0; j<M; j++)
		{
			s[j]=spin[j];
		}   // in order to avoid the change to the original array
		
		Random SHrand= new Random (10);
		
		for(int prestep=0; prestep< presteplimit; prestep++)
		{
			MCS(range, constJ, s, SHrand, 1, temp, field);
		}
		
		
	}
	
	public void findTcfromCv()
	{
		
	}
	
	
	
	
}