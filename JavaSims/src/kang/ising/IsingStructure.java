package kang.ising;

import chris.util.Random;



public class IsingStructure{
	
	//structures and parameters
	public int spin[];
	public int initialcopy[];   //the array of the initial copy of the system
	public int biaslabel[];
	
	
	public int L1, L2, M; //parameters for the lattice                                                                                                                                                                                                                                                                                                                                                                                                        
	public double J;     //interaction constant after normalization
	public double NJ;    //interaction constant before normalization
	public double percent;  //dilution percent
	public int deadsites;  //the number of diluted sites in the system
	
	//private double magnetization;
	//private double totalenergy;
	private double totalintenergy;
	private int totalspin;
	
	public int R;   //interaction range R=0 is NN interaction
	
	public int biasA;
	public int biasB;
	public double biaspercent;
	
	//the function for this class IsingStructure
	
	public IsingStructure(int L1, int L2, int R, double NJ, double percent, double biaspercent)     //generating function
	{
		this.L1=L1;
		this.L2=L2;
		this.M=L1*L2;
		this.R=R;
		this.percent=percent;
		this.biaspercent=biaspercent;

		this.NJ=NJ;
		if(R==0)
			this.J=NJ/4;
		if(R>0)
			this.J=NJ/((2*R+1)*(2*R+1)-1);
		this.spin=new int [M];
		this.initialcopy=new int [M];
		this.biaslabel=new int [M];
		
	}
	
	public IsingStructure clone(IsingStructure IS)
	{
		IsingStructure copy= new IsingStructure(IS.L1, IS.L2, IS.R, IS.NJ, IS.percent, IS.biaspercent);
		for(int t=0;t<IS.M; t++)
		{
			copy.spin[t]=IS.spin[t];
			copy.initialcopy[t]=IS.initialcopy[t];
			copy.biaslabel[t]=IS.biaslabel[t];
		}
		copy.biasA=IS.biasA;
		copy.biasB=IS.biasB;
		copy.totalspin=IS.totalspin;
		copy.totalintenergy=IS.totalintenergy;
		
		
		return copy;
	}
	
	//basic functions

	public int X(int bx)
	{
		int realx=bx;
		if (bx>L1)
			realx=bx-L1;
		if (bx<0)
			realx=bx+L1;
		return realx;
	}
	
	public int Y(int by)
	{
		int realy=by;
		if (by>L2)
			realy=by-L2;
		if (by<0)
			realy=by+L2;
		return realy;
	}
	
	public int distance (int a, int b)     // the code to calculate the distance between two points on the lattice
	{
		int dis=0;
		int ax, ay, bx, by;
		int dx2, dy2;
		ax= a/L2;
		ay= a%L2;
		bx= b/L2;
		by= b%L2;
		
		dx2=(ax-bx)*(ax-bx);
		dy2=(ay-by)*(ay-by);
		if((ax-bx+L1)*(ax-bx+L1)<(ax-bx)*(ax-bx))
			dx2=(ax-bx+L1)*(ax-bx+L1);
		if((ax-bx-L1)*(ax-bx-L1)<(ax-bx)*(ax-bx))
			dx2=(ax-bx-L1)*(ax-bx-L1);
		if((ay-by+L2)*(ay-by+L2)<(ay-by)*(ay-by))
			dy2=(ay-by+L2)*(ay-by+L2);
		if((ay-by-L2)*(ay-by-L2)<(ay-by)*(ay-by))
			dy2=(ay-by-L2)*(ay-by-L2);

		dis=dx2+dy2;
		return dis;
	}
	
	public void rectangle(int label[], int a, int b, int cx, int cy)  //draw a rectangle of 2a*2b at (cx,cy)
	{
		int bx, by;
		int x,y;
		
		for(bx=cx-a; bx<cx+a; bx++)
			for(by=cy-b; by<cy+b; by++)
			{
				x=X(bx);
				y=Y(by);
				label[x*L2+y]=1;
			}
	}
	
 	public int Nneighber(int a,int i )// function for the index of nearest neighbor
 	{
		int nx,ny; //index for neighbor
		int ni=0;
		nx=(int)i/L2;
		ny=(int)i%L2;
		
		if (a==0) {
			ni=nx*L2+ny-1;
			if  (ny==0) {
				ni=nx*L2+ny+L2-1;
		}
			
		}//(x,y-1) up
		
     	if (a==1){
			ni=(nx+1)*L2+ny;
			if  (nx==L1-1) {
				ni=(nx-L1+1)*L2+ny;
			}
			
		}//(x+1,y) right
		
		if (a==2){
			ni=nx*L2+ny+1;
			if  (ny==L2-1) {
				ni=nx*L2+ny-L2+1;
			}
			
		}//(x,y+1) down
		
		if (a==3){
			ni=(nx-1)*L2+ny;
			if  (nx==0) {
				ni=(nx+L1-1)*L2+ny;
			}
		}//(x-1,y) left
		
		return ni;
		
	}
	
 	public double interactionEchange (int j)//function for interaction energy
 	{ 
		double Energy=0;
		double Energychange=0;
		if (R==0)
		{
			int b,k;
		    for(b=0; b<4;b++)
		    {
			k=Nneighber(b,j);
			Energy=Energy+J*spin[j]*spin[k];
		    }
		    Energychange=-2*Energy;
		    
		}
		if (R!=0)
		{
			int S=0;
			int nx=j/L2;
			int ny=j%L2;
			int kx, ky;

			for (int m=-R; m<=R; m++)
				for (int n=-R; n<=R; n++)
				{
					kx=X(nx+m);
					ky=Y(ny+n);
					S+=spin[kx*L2+ky];	
				}
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
				totalintenergy+=IntEchange;
				
			}
			
			else
			{
				if(flip.nextDouble()<=Math.exp(-Echange/temperature))
						{
					            spin[j]=-spin[j];
								totalspin-=tempspin*2;
								totalintenergy+=IntEchange;	
						}
			}
		}

		
	}
 	
 	public int totalspin()
 	{
 		int total=0;
 		for(int k=0;k<M;k++)
 		{
 			total+=spin[k];
 		}
 		totalspin=total;
 		return total;
 	}
 	
 	public double magnetization()
 	{
 		double m=0;
 		m=totalspin/M;
 		return m;
 	}
 	
 	public double SpinIntEnergy(int j)
 	{
 		double SpinE=0;
 		if(spin[j]!=0)
 		{
 			if (R==0)
 			{
 
 				int b,k;
 			    for(b=0; b<4;b++)
 			    {
 				k=Nneighber(b,j);
 				SpinE+=J*spin[j]*spin[k];
 			    } 
 			}
 			
 			if (R!=0)
 			{
 				int S=0;
 				int nx=j/L2;
 				int ny=j%L2;
 				int kx, ky;
 				
 				for (int m=-R; m<=R; m++)
 					for (int n=-R; n<=R; n++)
 					{
 						kx=X(nx+m);
 						ky=Y(ny+n);
 						S+=spin[kx*L2+ky];	
 					}
 				SpinE=J*spin[j]*S-J;
 				
 			}
 			
 			
 		}
 		
 		return SpinE;
 	}
 	
 	public double totalIntEnergy()
 	{
 		double TotalE=0;
 		for(int k=0; k<M;k++)
 			TotalE+=SpinIntEnergy(k);
 		return TotalE;
 	}
 	
  	public double totalEnergy(double field)
 	{
 	    return totalintenergy+field*totalspin;	
 	}
 	
 	
 	
}