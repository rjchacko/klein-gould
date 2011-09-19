package kang.ising.BasicStructure;

import chris.util.Random;

public class FCIsingWithDamage {

	public int N;
	public int Ns;
	public int Nx;
	public double Nsu, Nsd;   //number of the up and down normal spins  
	public double Ms;     //Ms=Nsu-Nsd

	public double m;
	public double Nxu, Nxd;   //number of the up and down damaged spins
	public double Mx;     //Mx=Nxu-Nxd

	public double q;      //percentage of the damaged spins
	public double p;      //p=1-q, percentage of the normal spins
	public double a;      //the moment of the damaged spins(-a, +a)
	
	
	public double J;
	
	public FCIsingWithDamage(int N)
	{
		this.N=N;
	}
	
	public FCIsingWithDamage clone()
	{
		FCIsingWithDamage copy=new FCIsingWithDamage(N);
		copy.Ns=Ns;
		copy.Nx=Nx;
		copy.Nsu=Nsu;
		copy.Nsd=Nsd;
		copy.Ms=Ms;
		copy.Nxu=Nxu;
		copy.Nxd=Nxd;
		copy.Mx=Mx;
		copy.m=m;
		copy.J=J;
		copy.p=p;
		copy.q=q;
		copy.a=a;
		
		return copy;
	}
	
	public void setJ(double NJ)
	{

		this.J=NJ/(N-1);
	
	}
	
	public void initialize(double q, double a)
	{
	    Nx=(int)(N*q);
	    Ns=N-Nx;
	    this.a=a;
	    this.q=q;
	    
		
		
		Nsu=Math.round(Ns/2);
		Nsd=Ns-Nsu;
		Ms=Nsu-Nsd;
		
		Nxu=Math.round(Nx/2);
		Nxd=Nx-Nxu;
		Mx=Nxu-Nxd;
		
	    this.m=(Ms+Mx*a)/N;
	    
	}
	
	public double Magnetization()
	{
		double M=0;
		M=(Ms+Mx*a)/N;
		return M;
	}

	public double Energy1(double H)
	{
		double E=0;
		E=J*((Nsu*(Nsu-1)+Nsd*(Nsd-1))/2-Nsu*Nsd+a*a*(Nxu*(Nxu-1)+Nxd*(Nxd-1))/2-a*a*Nxu*Nxd+a*Nsu*Nxu+a*Nsd*Nxd-a*Nsu*Nxd-a*Nxu*Nsd)-H*(Nsu-Nsd+a*Nxu-a*Nxd);
		return E;
	}
	
	public double Energy2(double H)
	{
		double E=0;
		E=J*(Ms*Ms/2-Ns/2+a*a*Mx*Mx/2-a*a*Nx/2+a*Ms*Mx)-H*(Ms+a*Mx);
		return E;
	}
	
	public void spinflip(Random spinrand, Random fliprand, double T, double H)   //spinrand for choosing the up or down spin, fliprand for determining the spinflip
	{
		double Echange=0;
		double randomnumber=spinrand.nextDouble();
		
		if(randomnumber*N<Nsu)                  //zone 1 choose a regular up spin and flip to down
		{
			Echange=-2*J*(Nsu-Nsd)+2*J-2*J*a*(Nxu-Nxd)+2*H;
			if(Echange<0)
			{
				Nsu--;
				Nsd++;
				Ms=Nsu-Nsd;
			}
			else
			{
				if(fliprand.nextDouble()<Math.exp(-Echange/T))
				{
					Nsu--;
					Nsd++;
					Ms=Nsu-Nsd;
				}
			}
		}
		
		else if(randomnumber*N<Ns)               //zone 2 choose a regular down spin and flip to up
		{
			Echange=2*J*(Nsu-Nsd)+2*J+2*J*a*(Nxu-Nxd)-2*H;
			if(Echange<0)
			{
				Nsd--;
				Nsu++;
				Ms=Nsu-Nsd;
			}
			else
			{
				if(fliprand.nextDouble()<Math.exp(-Echange/T))
				{
					Nsd--;
					Nsu++;
					Ms=Nsu-Nsd;
				}
			}
		}
		
		else if(randomnumber*N<Ns+Nxu)               //zone 3 choose a damaged up spin and flip to down
		{
			Echange=-2*J*(Nxu-Nxd)*a*a+2*J*a*a-2*J*a*(Nsu-Nsd)+2*H*a;
			if(Echange<0)
			{
				Nxu--;
				Nxd++;
				Mx=Nxu-Nxd;
			}
			else
			{
				if(fliprand.nextDouble()<Math.exp(-Echange/T))
				{
					Nxu--;
					Nxd++;
					Mx=Nxu-Nxd;
				}
			}
		}
		
		else               //zone 4 choose a damaged down spin and flip to up
		{
			Echange=2*J*(Nxu-Nxd)*a*a+2*J*a*a+2*J*a*(Nsu-Nsd)-2*H*a;
			if(Echange<0)
			{
				Nxd--;
				Nxu++;
				Mx=Nxu-Nxd;
			}
			else
			{
				if(fliprand.nextDouble()<Math.exp(-Echange/T))
				{
					Nxd--;
					Nxu++;
					Mx=Nxu-Nxd;
				}
			}
		}
		
		m=(Ms+a*Mx)/N;
	}
	
	
	public void spinflipGlauber(Random spinrand, Random fliprand, double T, double H)   //Glauber dynamics
	{
		double Echange=0;
		double Pglauber=0;
		double randomnumber=spinrand.nextDouble();
		
		if(randomnumber*N<Nsu)                  //zone 1 choose a regular up spin and flip to down
		{
			Echange=-2*J*(Nsu-Nsd)+2*J-2*J*a*(Nxu-Nxd)+2*H;
			Pglauber=1/(1+Math.exp(Echange/T));

				if(fliprand.nextDouble()<Pglauber)
				{
					Nsu--;
					Nsd++;
					Ms=Nsu-Nsd;
				}

		}
		
		else if(randomnumber*N<Ns)               //zone 2 choose a regular down spin and flip to up
		{
			Echange=2*J*(Nsu-Nsd)+2*J+2*J*a*(Nxu-Nxd)-2*H;
			Pglauber=1/(1+Math.exp(Echange/T));
			
				if(fliprand.nextDouble()<Pglauber)
				{
					Nsu++;
					Nsd--;
					Ms=Nsu-Nsd;
				}
			
		}
		
		else if(randomnumber*N<Ns+Nxu)               //zone 3 choose a damaged up spin and flip to down
		{
			Echange=-2*J*(Nxu-Nxd)*a*a+2*J*a*a-2*J*a*(Nsu-Nsd)+2*H*a;
			Pglauber=1/(1+Math.exp(Echange/T));
			
				if(fliprand.nextDouble()<Pglauber)
				{
					Nxu--;
					Nxd++;
					Mx=Nxu-Nxd;
				}
			
		}
		
		else               //zone 4 choose a damaged down spin and flip to up
		{
			Echange=2*J*(Nxu-Nxd)*a*a+2*J*a*a+2*J*a*(Nsu-Nsd)-2*H*a;
			Pglauber=1/(1+Math.exp(Echange/T));
			
				if(fliprand.nextDouble()<Pglauber)
				{
					Nxu++;
					Nxd--;
					Mx=Nxu-Nxd;
				}
			
		}
		
		m=(Ms+a*Mx)/N;
	}
	
	
	
	
	public void MCS(String dynamics, Random spinrand, Random fliprand, double T, double H, double ratio)
	{
		int flipnumber=(int)(ratio*N);
		for(int j=0; j<flipnumber; j++)
		{
			if(dynamics=="Metropolis")
				spinflip(spinrand, fliprand, T, H);
			if(dynamics=="Glauber")
				spinflipGlauber(spinrand, fliprand, T, H);
		}
	}
	
	
	
	
}