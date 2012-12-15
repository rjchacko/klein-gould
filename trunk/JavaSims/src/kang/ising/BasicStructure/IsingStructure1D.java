package kang.ising.BasicStructure;


import kang.ising.BasicStructure.StructureFactor;
import chris.util.Random;
//import chris.util.PrintUtil;

public class IsingStructure1D{
	
	//structures and parameters
	public int spin[];
	public int initialcopy[];   //the array of the initial copy of the system
	public int biaslabel[];
	public double dilutionmap[];
	public int largestcluster[];
	
	
	public StructureFactor SFup;
	public StructureFactor SFdown;
	public StructureFactor SFdilution;
	
	public double upspin[];
	public double downspin[];
	public double damage[];
	
	
	
	public int L; //parameters for the lattice                                                                                                                                                                                                                                                                                                                                                                                                        
	public double J;     //interaction constant after normalization
	public double NJ;    //interaction constant before normalization
	public double percent;  //dilution percent
	public int deadsites;  //the number of diluted sites in the system
	
	
	public double totalintenergy;
	public int totalspin;
	public double magnetization;
	
	public int R;   //interaction range R=0 is NN interaction
	
	public int biasA;
	
	public double biaspercent;
	
	//the function for this class IsingStructure
	
	public IsingStructure1D(int L, int R, double NJ, double percent, double biaspercent)     //generating function
	{
		this.L=L;
		
		this.R=R;
		this.percent=percent;
		this.biaspercent=biaspercent;
	
		this.deadsites=0;

		this.NJ=NJ;
		
		
			
		this.J=NJ/(2*R);
		
			
				
		this.spin=new int [L];
		this.initialcopy=new int [L];
		this.biaslabel=new int [L];
		this.dilutionmap=new double[L];
		this.damage=new double[L];
		
	}
	
	public IsingStructure1D clone()
	{
		IsingStructure1D copy= new IsingStructure1D(L, R, NJ, percent, biaspercent);
		for(int t=0;t<L; t++)
		{
			copy.spin[t]=spin[t];
			copy.initialcopy[t]=initialcopy[t];
			copy.biaslabel[t]=biaslabel[t];
			copy.dilutionmap[t]=dilutionmap[t];
			copy.damage[t]=damage[t];

		}
		copy.biasA=biasA;
	
		copy.totalspin=totalspin;
		copy.totalintenergy=totalintenergy;
		copy.deadsites=deadsites;
	
		return copy;
	}
	
	public IsingStructure1D Dperturbation(int seed)
	{
		IsingStructure1D perturb=new IsingStructure1D(L, R, NJ, 0, 0);
		perturb.Dinitialization(1, 1, 0);
		perturb.percent=percent;
		perturb.biaspercent=biaspercent;
		
		Random DPrand=new Random(seed);
		int j;
		int i=deadsites;
	

		while(i>0)
		{
			j=(int)(DPrand.nextDouble()*L);
			if(perturb.spin[j]!=0)
			{
				perturb.spin[j]=0;
				perturb.damage[j]=1;
				perturb.deadsites++;
				i--;
			}
		}
		perturb.Sinitialization(1, 1);
		
		return perturb;
	}
	
	public void Dinitialization(int Dseed, int Bseed, int A )//dilution initialization
	{
		Random Drand= new Random(Dseed);
		Random Brand= new Random(Bseed);
		biasA=A;
		
		deadsites=0;

		
		for (int t=0;t<L;t++)
			{
			spin[t]=1;
			biaslabel[t]=0;
			dilutionmap[t]=0;
			damage[t]=0;
			}
			
		
		int cx; //x,y index for the center
		cx=(int)L/2;
		

		
		for(int j=0; j<L; j++)
		{
			if (biaslabel[j]==1)
				if(Brand.nextDouble()<biaspercent)
					{
					spin[j]=0;
					damage[j]=1;
					deadsites++;
					}
			if (biaslabel[j]==0)
				if(Drand.nextDouble()<percent)
					{
					spin[j]=0;
					damage[j]=1;
					deadsites++;
					}
		}
		

	}
	

	
	public void Sinitialization(int type, int Sseed)//spin config initialization (type: 0-random 1-spin up   2-spin down)
	{
		
		
		Random spinrand= new Random(Sseed);
		totalintenergy=0;
		totalspin=0;
		
		
		if(type==0)
			for (int i=0; i<L; i++)
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
			for (int i=0; i<L; i++)
		    {
			   if(spin[i]!=0)
				   {
				   spin[i]=1;
				   
				   }
		    }
		
		if(type==2)
			for (int i=0; i<L; i++)
		    {
			   if(spin[i]!=0)
				   {
				   spin[i]=-1;
				   
				   }
		    }
		
		for (int i=0; i<L; i++)
		{

			initialcopy[i]=spin[i];  // here, make the copy of the system
			totalspin+=spin[i];      // calculate the total spin at the beginning
		}
		
		totalintenergy= TotalIntEnergy();
		magnetization=Magnetization();
		
	}
	
	public void SpinSetup()
	{
		upspin=new double[L];
		downspin=new double[L];
		for(int j=0; j<L;j++)
		{
			upspin[j]=0;
			downspin[j]=0;
			if(spin[j]==1)
				upspin[j]=1;
			else if(spin[j]==-1)
				downspin[j]=1;	
		}
	}
	

	//basic functions
	

	


	public int X(int bx)
	{
		int realx=bx;
		if (bx>=L)
			realx=bx-L;
		if (bx<0)
			realx=bx+L;
		return realx;
	}
	

	



	
	public int SumInRange(int spin[], int j, int R)  //sum in the range of 2R+1
	{
		int S=0;
		int nx=j;
		
		int kx;

		{
			for (int m=-R; m<=R; m++)
			    {
				   kx=X(nx+m);
				  
				   S+=spin[kx];	
			    }
		}
		
		return S;
	}
		
	public void line(int label[], int a, int cx)  //draw a line of 2a at cx
	{
		int bx;
		int x;
		
		for(bx=cx-a; bx<cx+a; bx++)
			{
				x=X(bx);
				label[x]=1;
			}
	}
	

 	public double interactionEchange (int j)//function for interaction energy
 	{ 
		double Energy=0;
		double Energychange=0;

	
		{
			int S=SumInRange(spin, j, R);
			
			Energy=J*spin[j]*S-J;
			Energychange=-2*Energy;
			
		}
		
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
				magnetization=Magnetization();
				totalintenergy+=IntEchange;
				
			}
			
			else
			{
				if(flip.nextDouble()<=Math.exp(-Echange/temperature))
						{
					            spin[j]=-spin[j];
								totalspin-=tempspin*2;
								magnetization=Magnetization();
								totalintenergy+=IntEchange;	
						}
			}
		}

	}
 	
 	public void GlauberSpinflip(int j, Random flip, double temperature, double field)
	{
		
 		if(spin[j]!=0)
		{
	 		double ZeemanE=2*field*spin[j]; //the change in Zeeman's energy if we flip the spin[j]
			double IntEchange=interactionEchange(j);
	 		double Echange=ZeemanE+IntEchange;
	 		double Pglauber=1/(Math.exp(Echange/temperature));
			int tempspin= spin[j];
	 		
			
			if(flip.nextDouble()<=Pglauber)
						{
					            spin[j]=-spin[j];
								totalspin-=tempspin*2;
								magnetization=Magnetization();
								//totalintenergy+=IntEchange;	
						}
			
		}

	}
 	
	public void MCS(double T, double H, Random flip, double ratio, String dynamics)
	{
	    int j=0;
	    
	    int rationalMCS= (int) (ratio*L);
	    for (int f=0; f< rationalMCS; f++)
	    {
		   j=(int) (flip.nextDouble()*L);
		   if(dynamics=="Metropolis")
			   MetropolisSpinflip(j, flip, T, H);
		   if(dynamics=="Glauber")
			   GlauberSpinflip(j, flip, T, H);
			   
	    }

	}
	
  	public int TotalSpin()
 	{
 		int total=0;
 		for(int k=0;k<L;k++)
 		{
 			total+=spin[k];
 		}
 		return total;
 	}
 	
 	public double SpinIntEnergy(int j)
 	{
 		double SpinE=0;
 		if(spin[j]!=0)
 		{
 			
 				int S=SumInRange(spin, j, R);
 				SpinE=J*spin[j]*S-J;
 		}
 		
 		return SpinE;
 	}
 	
 	public double TotalIntEnergy()
 	{
 		double TotalE=0;
 		for(int k=0; k<L;k++)
 			TotalE+=SpinIntEnergy(k);
 		totalintenergy=TotalE;
 		return TotalE;
 	}
 	
  	public double TotalEnergy(double field)
 	{
 	    return totalintenergy+field*totalspin;	
 	}
  	
  	public double Magnetization()
  	{
  		double m;
  		m=((double)totalspin/L);
  		return m;
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