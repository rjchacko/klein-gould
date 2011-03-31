package kang.ising.BasicStructure;

import chris.util.Random;

public class FCIsing {

	public int N;
	public double Nu, Nd;
	public double M;
	public double m;
	public double J;
	
	public FCIsing(int N)
	{
		this.N=N;
	}
	
	public FCIsing clone()
	{
		FCIsing copy=new FCIsing(N);
		copy.Nu=Nu;
		copy.Nd=Nd;
		copy.M=M;
		copy.m=m;
		copy.J=J;
		return copy;
	}
	
	public void initialize(double m)
	{
		Nu=Math.round((1+m)*N/2);
		Nd=N-Nu;
		M=Nu-Nd;
	    this.m=M/N;
	    assert (this.m==m);
	}
	
	public void spinflip(Random spinrand, Random fliprand, double T, double H)   //spinrand for choosing the up or down spin, fliprand for determining the spinflip
	{
		double Echange=0;
		if(spinrand.nextDouble()<(Nu/N))
		{
			Echange=-2*J*(Nu-Nd)+2*J+2*H;
			if(Echange<0)
			{
				Nu--;
				Nd++;
				M=Nu-Nd;
			}
			else
			{
				if(fliprand.nextDouble()<Math.exp(-Echange/T))
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
			if(Echange<0)
			{
				Nu++;
				Nd--;
				M=Nu-Nd;
			}
			else
			{
				if(fliprand.nextDouble()<Math.exp(-Echange/T))
				{
					Nu++;
					Nd--;
					M=Nu-Nd;
				}
			}
			
		}
	}
	
	public void MCS(Random spinrand, Random fliprand, double T, double H, double ratio)
	{
		int flipnumber=(int)(ratio*N);
		for(int j=0; j<flipnumber; j++)
		{
			spinflip(spinrand, fliprand, T, H);
		}
	}
	
	
}
