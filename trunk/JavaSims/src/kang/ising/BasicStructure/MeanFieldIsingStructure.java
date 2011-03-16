package kang.ising.BasicStructure;

import chris.util.Random;

public class MeanFieldIsingStructure{// the class of fully connected ising model
	
	public int spin[];
	public int initialcopy[];
	
	public int L1,L2,M;
	
	public double J;
	public double NJ;  // the interaction constant before renormalization
	
	public int totalspin;
	public double totalintenergy;
	
	public MeanFieldIsingStructure(int L1, int L2, double NJ)
	{
		this.L1=L1;
		this.L2=L2;
		this.M=L1*L2;
		this.NJ=NJ;
		this.J=NJ/(M-1);     //rescale the interaction constant
		
		this.spin= new int[M];
		this.initialcopy= new int[M];
		
	}
	
	public MeanFieldIsingStructure clone()
	{
		MeanFieldIsingStructure copy= new MeanFieldIsingStructure(L1, L2, NJ);
		for(int j=0; j<M; j++)
		{
			copy.spin[j]=spin[j];
			copy.initialcopy[j]=initialcopy[j];
		}
		copy.totalspin=totalspin;
		copy.totalintenergy=totalintenergy;
		
		
		
		return copy;
	}
	
	public void Sinitialization(int type, int Sseed)//spin config initialization (type: 0-random 1-spin up   2-spin down)
	{
		Random spinrand= new Random(Sseed);
		totalintenergy=0;
		totalspin=0;
		
		
		if(type==0)
			for (int i=0; i<M; i++)
		    {
			   if(spin[i]!=0)
			   {
				   spin[i]=-1;
			   if (spinrand.nextDouble()> 0.5)
				   {
				   spin[i]=1;

				   }
			   
			   }
				
		    }
		if(type==1)
			for (int i=0; i<M; i++)
		    {
			   if(spin[i]!=0)
				   {spin[i]=1;

				   }
		    }
		
		if(type==2)
			for (int i=0; i<M; i++)
		    {
			   if(spin[i]!=0)
				   {
				   spin[i]=-1;

				   }
		    }
		
		for (int i=0; i<M; i++)
		{

			initialcopy[i]=spin[i];  // here, make the copy of the system
			totalspin+=spin[i];      // calculate the total spin at the beginning
		}
		
		
		
	}
	
	
	
	public double interactionEchange (int j)//function for interaction energy
 	{ 
		double Energy=0;
		double Energychange=0;

		Energy=spin[j]*J*(totalspin-spin[j]);
		Energychange=-2*Energy;
		return Energychange;	
    }
	
	public void MetropolisSpinflip(int j, Random flip, double temperature, double field)
	{
		
 		if(spin[j]!=0)
		{
	 		double ZeemanE=2*field*spin[j]; //the change in Zeeman's energy if we flip the spin[j]
			double IntEchange=interactionEchange(j);
	 		double Echange=ZeemanE+IntEchange;
			int tempspin= spin[j];
	 		
			if(Echange<0)
			{
				spin[j]=-spin[j];
				totalspin-=tempspin*2;
				
			}
			
			else
			{
				if(temperature>10)
				{
					spin[j]=-spin[j];
					totalspin-=tempspin*2;
				}
				else if(flip.nextDouble()<=Math.exp(-Echange/temperature))
						{
					            spin[j]=-spin[j];
								totalspin-=tempspin*2;

						}
			}
		}

		
	}
	
	public void MCS(double T, double H, Random flip, double ratio)
	{
	    int j=0;
	    
	    int rationalMCS= (int) (ratio*M);
	    for (int f=0; f< rationalMCS; f++)
	    {
		   j=(int) (flip.nextDouble()*M);
		   MetropolisSpinflip(j, flip, T, H);
	    }

	}
	
	public double Magnetization()
	{
  		double m;
  		m=((double)totalspin/M);
  		return m;
	}
	
 	public double TotalIntEnergy()
 	{
 		double TotalE=0;
 	    TotalE=(totalspin*totalspin-M)*J;
 		return TotalE;
 	}
 	
  	public double TotalEnergy(double field)
 	{
 	    return TotalIntEnergy()+field*totalspin;	
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