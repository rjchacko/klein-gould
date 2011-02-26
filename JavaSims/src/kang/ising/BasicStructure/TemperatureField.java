package kang.ising.BasicStructure;

import kang.ising.BasicStructure.IsingStructure;
import chris.util.Random;


public class TemperatureField{
	
	public IsingStructure IS;
	public IsingStructure PIS;
	public double Tfield[];
	
	public TemperatureField(IsingStructure Ising)
	{
		this.IS=Ising.clone();
		this.PIS=Ising.clone();
		this.Tfield=new double[Ising.M];
		
		for(int j=0;j<PIS.M; j++)
		{
			if(PIS.spin[j]==0)
				{
				PIS.spin[j]=1;	
				}
		}// get rid of the dilution
		//PIS.totalintenergy=PIS.TotalIntEnergy();
		PIS.totalspin=PIS.TotalSpin();
		PIS.magnetization=PIS.Magnetization();
		
	    //IS.dilutionmap(IS.R);
		
	}
	
	public void PMCS(double T, double H, Random flip, double ratio)
	{
	    int j=0;
	    
	    int rationalMCS= (int) (ratio*PIS.M);
	    for (int f=0; f< rationalMCS; f++)
	    {
		   j=(int) (flip.nextDouble()*PIS.M);
		   if(IS.spin[j]!=0)
		   {
			   PIS.MetropolisSpinflip(j, flip, T, H);
		   }
		   else  //the sites corresponding to the diluted sites always flip
		   {
			    int tempspin=PIS.spin[j];
				PIS.spin[j]=-PIS.spin[j];
				PIS.totalspin-=tempspin*2;
				PIS.magnetization=PIS.Magnetization();
		   }
	    }

	}
	
	public void GeneratingT(double T, double dilutionmap[])
	{
		for(int i=0;i<IS.M; i++)
		{
			Tfield[i]=T/(1-dilutionmap[i]);
		}
	}
	
	public void TMCS(double T, double H, Random flip, double ratio)
	{
		int j=0;
		
	    int rationalMCS= (int) (ratio*PIS.M);
	    for (int f=0; f< rationalMCS; f++)
	    {
		   j=(int) (flip.nextDouble()*PIS.M);
		   PIS.MetropolisSpinflip(j, flip, Tfield[f], H);
	    }
	}

	
	
	
	
	
	
}