package kang.ising.BasicStructure;



import chris.util.Random;
//import chris.util.PrintUtil;

public class J1J2Structure{
	
	//structures and parameters
	public int spin[];
	public int initialcopy[];   //the array of the initial copy of the system
	public int display[];  //the array to translate the system's spin configuration to color
	
	public int biaslabel[];
	//public double dilutionmap[];
	
	//J1-ferromagnetic (J1<0)
	public double J1;     //nearest neighbor interaction constant after normalization
	public double NJ1;    //nearest neighbor interaction constant before normalization
	//J2-antiferromagnetic (J2>0)
	public double J2;     //next nearest neighbor interaction constant after normalization
	public double NJ2;    //next nearest neighbor interaction constant before normalization
	
	
	public int L1, L2, M; //parameters for the lattice  
	public double percent;  //dilution percent
	public int deadsites;  //the number of diluted sites in the system
	
	
	public double totalintenergy;
	public int totalspin;
	public int totalspinX;
	public int totalspinY;
	
	public double magnetization;
	public double mx,my;
	public double mm2;
	

	public int biasA;
	public int biasB;
	public double biaspercent;
	
	//the function for this class IsingStructure
	
	public J1J2Structure(int L1, int L2, double NJ1, double NJ2, double percent, double biaspercent)     //generating function
	{
		this.L1=L1;
		this.L2=L2;
		this.M=L1*L2;
		
		this.percent=percent;
		this.biaspercent=biaspercent;
		
		this.deadsites=0;

		this.NJ1=NJ1;
		this.NJ2=NJ2;
		
		this.J1=NJ1/4;
		this.J2=NJ2/4;
				
		this.spin=new int [M];
		this.initialcopy=new int [M];
		
		this.biaslabel=new int [M];
		this.display=new int [M];
		//this.dilutionmap=new double[M];
		
	}
	
	public J1J2Structure clone()
	{
		J1J2Structure copy= new J1J2Structure(L1, L2, NJ1, NJ2, percent, biaspercent);
		for(int t=0;t<M; t++)
		{
			copy.spin[t]=spin[t];
			copy.initialcopy[t]=initialcopy[t];
			copy.biaslabel[t]=biaslabel[t];
			//copy.dilutionmap[t]=dilutionmap[t];
		}
		copy.biasA=biasA;
		copy.biasB=biasB;
		copy.totalspin=totalspin;
		copy.totalspinX=totalspinX;
		copy.totalspinY=totalspinY;
		copy.totalintenergy=totalintenergy;
		copy.deadsites=deadsites;
		copy.mx=mx;
		copy.my=my;
		copy.mm2=mm2;
		
		
		return copy;
	}
	
	/*public IsingStructure Dperturbation(int seed)
	{
		IsingStructure perturb=new IsingStructure(L1, L2, R, NJ, 0, 0, shape);
		perturb.Dinitialization(1, 1, 0, 0);
		perturb.percent=percent;
		perturb.biaspercent=biaspercent;
		
		Random DPrand=new Random(seed);
		int j;
		int i=deadsites;
	

		while(i>0)
		{
			j=(int)(DPrand.nextDouble()*L1*L2);
			if(perturb.spin[j]!=0)
			{
				perturb.spin[j]=0;
				perturb.deadsites++;
				i--;
			}
		}
		perturb.Sinitialization(1, 1);
		
		return perturb;
	}*/
	
	
	public void Dinitialization(int Dseed, int Bseed, int A, int B )//dilution initialization
	{
		Random Drand= new Random(Dseed);
		Random Brand= new Random(Bseed);
		biasA=A;
		biasB=B;
		deadsites=0;

		
		for (int t=0;t<M;t++)
			{
			spin[t]=1;
			biaslabel[t]=0;
			//dilutionmap[t]=0;
			}
			
		
		int cx, cy; //x,y index for the center
		cx=(int)L1/2;
		cy=(int)L2/2;
		if(biaspercent!=percent)
		    rectangle(biaslabel, A,B, cx,cy);
		
		for(int j=0; j<M; j++)
		{
			if (biaslabel[j]==1)
				if(Brand.nextDouble()<biaspercent)
					{
					spin[j]=0;
					deadsites++;
					}
			if (biaslabel[j]==0)
				if(Drand.nextDouble()<percent)
					{
					spin[j]=0;
					deadsites++;
					}
		}
		
		/*or(int k=0; k<M; k++)
		{
			dilutionmap[k]=dilutionratio(R,k);
		} // draw the dilution map
	*/
	}
	
	public void Sinitialization(int type, int Sseed)//spin config initialization (type: 0-random 1-spin up   2-spin down)
	{
		Random spinrand= new Random(Sseed);
		totalintenergy=0;
		totalspin=0;
		totalspinX=0;
		totalspinY=0;
		
		
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
			totalspinX+=Xsign(i)*spin[i];
			totalspinY+=Ysign(i)*spin[i];
		}
		
		totalintenergy= TotalIntEnergy();
		magnetization=Magnetization();
		mx=Mx();
		my=My();
		mm2=MM2();
		
		
	}
	
	//basic functions
	
	
	/*public double dilutionratioCircle(int r, int i)    //calculate the dilution ratio within the circle with radius r
	{
		double ratio;
		double dilutedsite;
		dilutedsite=0;
		double totalinrange;
		totalinrange=0;
		int j;
		int disij;

		for(j=0; j<M;j++)
		{
            disij= distance(i,j);

			if(disij<=r*r)
			{
				totalinrange++;
				if(spin[j]==0)
					dilutedsite++;
			}
		}
	
		ratio=dilutedsite/totalinrange;
		return ratio;
	}*/

	public int X(int bx)
	{
		int realx=bx;
		if (bx>=L1)
			realx=bx-L1;
		if (bx<0)
			realx=bx+L1;
		return realx;
	}
	
	public int Y(int by)
	{
		int realy=by;
		if (by>=L2)
			realy=by-L2;
		if (by<0)
			realy=by+L2;
		return realy;
	}
	
	public int distance (int a, int b)     // the code to calculate the square of the distance between two points on the lattice
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

	public int CenterOfMass(String CGshape)              // return the index in the array corresponding to the center of mass (with specific coarse graining shape)
	{
		int center=0;
		
		int totalM=TotalSpin();	
		int minorityspin=0;
		// now determine which direction is the minority direction
		if(totalM>0)
			minorityspin=-1;
		if(totalM<0)
			minorityspin=1;
		// now record all the sites with minority direction
		int map[]= new int[M];
		int CGmap[]= new int[M];    // the coarse grained map[]
		int totalMinorityspins=0;
		
		for(int mi=0; mi<M; mi++)
		{
			if(spin[mi]==minorityspin)
				{
				     map[mi]=1;
				     totalMinorityspins++;     //find out how many minority spins are there in the whole lattice and then determine the range of CG
				}
			else
				map[mi]=0;
		}
		
		int CGR= (int)(Math.sqrt(totalMinorityspins)/2);   //estimate the best range for CG
		int MAX=0;
		int maxi=0;
		for(int ci=0; ci<M; ci++)
		{
			CGmap[ci]=SumInRange(map,ci,CGR, CGshape);
			if(CGmap[ci]>MAX)
				{
				     maxi=ci;
				     MAX=CGmap[ci];
				}
		}
		center=maxi;
		
		return center;
		
	}
	
	public int SumInRange(int spin[], int j, int R, String shape)  //sum in the range of (2R+1)*(2R+1) square or (R+1)^2-R^2 diamond
	{
		int S=0;
		int nx=j/L2;
		int ny=j%L2;
		int kx, ky;

		if(shape=="square")
		{
			for (int m=-R; m<=R; m++)
				for (int n=-R; n<=R; n++)
			    {
				   kx=X(nx+m);
				   ky=Y(ny+n);
				   S+=spin[kx*L2+ky];	
			    }
		}
		if(shape=="diamond")
		{
			for (int m=-R; m<=R; m++)
				for (int n=-(R-Math.abs(m)); n<=(R-Math.abs(m)); n++)
				{
					kx=X(nx+m);
					ky=Y(ny+n);
					S+=spin[kx*L2+ky];	
				}
		}
		
		return S;
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
	
 	public int NextNneighber(int a,int i )// function for the index of nearest neighbor
 	{
		int nx,ny; //index for neighbor
		int ni=0;
		nx=(int)i/L2;
		ny=(int)i%L2;
		
		if (a==0) {
			ni=X(nx+1)*L2+Y(ny-1);
		
		}
			
		//(x+1,y-1) up-right
		
     	if (a==1){
			ni=X(nx+1)*L2+Y(ny+1);
		
		}//(x+1,y+1) down-right
		
		if (a==2){
			ni=X(nx-1)*L2+Y(ny+1);
			
		}//(x-1,y+1) down-left
		
		if (a==3){
			ni=X(nx-1)*L2+Y(ny-1);
		
		}//(x-1,y-1) up-left
		
		return ni;
		
	}
 	
	public int Nneighber(int a,int i )// function for the index of next nearest neighbor
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
 	
	
	
	//the display functions to convert spin configuration into color maps
 	public int BlockSum(int j)
 	{
 		int jx=j/L2;
 		int jy=j%L2;
 		int sum=0;
 		int a=spin[jx*L2+jy];
 		int b=spin[X(jx+1)*L2+jy];
 		int c=spin[jx*L2+Y(jy+1)];
 		int d=spin[X(jx+1)*L2+Y(jy+1)];
 		sum=a+b+c+d;
 		return sum;
 	}
 	
 	public int XBlockSumX(int j)
 	{
 		int jx=j/L2;
 		int jy=j%L2;
 		int sum=0;
 		int a=spin[jx*L2+jy];
 		int b=spin[X(jx+1)*L2+jy];
 		int c=spin[jx*L2+Y(jy+1)];
 		int d=spin[X(jx+1)*L2+Y(jy+1)];
 		sum=Xsign(j)*(a-b+c-d);
 		return sum;
 	}
 	
 	public int YBlockSumY(int j)
 	{
 		int jx=j/L2;
 		int jy=j%L2;
 		int sum=0;
 		int a=spin[jx*L2+jy];
 		int b=spin[X(jx+1)*L2+jy];
 		int c=spin[jx*L2+Y(jy+1)];
 		int d=spin[X(jx+1)*L2+Y(jy+1)];
 		sum=Ysign(j)*(a+b-c-d);
 		return sum;
 	}
	
 	public void Display()    //the function to generate the color map
 	{
 		int valueX, valueY;
 		for(int j=0; j<M; j++)
 		{
 			valueX=XBlockSumX(j);
 			valueY=YBlockSumY(j);
 			if(BlockSum(j)==0)
 			{
 				if((valueX==0)&(valueY==-4))
 					display[j]=1;    //1 is red
 				else if((valueX==4)&(valueY==0))
 					display[j]=2;    //2 is blue
 				else if((valueX==-4)&(valueY==0))
 					display[j]=3;    //3 is black
 				else if((valueX==0)&(valueY==4))
 					display[j]=4;    //4 is green
 				else
 					display[j]=-2;   //-2 is the cross configuration(gray), which is not a ground state either
 			}
 			else
 				display[j]=-1;    //-1 is white
 		}
 			
 	}
 	
 	/////
  	public double interactionEchange (int j)//function for interaction energy
 	{ 
		double Energy=0;
		double Energychange=0;
		
		{
			int b,kn,nn;
		    for(b=0; b<4;b++)
		    {
			kn=Nneighber(b,j);
			nn=NextNneighber(b,j);
			
			Energy=Energy+J1*spin[j]*spin[kn]+J2*spin[j]*spin[nn];
		    }
		    Energychange=-2*Energy;
		    
		}

		
		return Energychange;	
    }

 	public void MetropolisSpinflip(int j, Random flip, double temperature, double[] field)
	{
		
 		if(spin[j]!=0)
		{
	 		double ZeemanE=2*field[j]*spin[j]; //the change in Zeeman's energy if we flip the spin[j]
			double IntEchange=interactionEchange(j);
	 		double Echange=ZeemanE+IntEchange;
			int tempspin= spin[j];
	 		
			if(Echange<0)
			{
				spin[j]=-spin[j];
				totalspin-=tempspin*2;
				totalspinX-=tempspin*2*Xsign(j);
				totalspinY-=tempspin*2*Ysign(j);
				magnetization=Magnetization();
				mx=Mx();
				my=My();
				mm2=MM2();
				totalintenergy+=IntEchange;
				
			}
			
			else
			{
				if(flip.nextDouble()<=Math.exp(-Echange/temperature))
						{
					            spin[j]=-spin[j];
								totalspin-=tempspin*2;
								totalspinX-=tempspin*2*Xsign(j);
								totalspinY-=tempspin*2*Ysign(j);
								magnetization=Magnetization();
								mx=Mx();
								my=My();
								mm2=MM2();
								totalintenergy+=IntEchange;	
								
						}
			}
		}

	}
 	
 	public void GlauberSpinflip(int j, Random flip, double temperature, double[] field)
	{
		
 		if(spin[j]!=0)
		{
	 		double ZeemanE=2*field[j]*spin[j]; //the change in Zeeman's energy if we flip the spin[j]
			double IntEchange=interactionEchange(j);
	 		double Echange=ZeemanE+IntEchange;
	 		double Pglauber=1/(Math.exp(Echange/temperature));
			int tempspin= spin[j];
	 		
			
			if(flip.nextDouble()<=Pglauber)
						{
					            spin[j]=-spin[j];
								totalspin-=tempspin*2;
								totalspinX-=tempspin*2*Xsign(j);
								totalspinY-=tempspin*2*Ysign(j);
								magnetization=Magnetization();
								mx=Mx();
								my=My();
								mm2=MM2();
								totalintenergy+=IntEchange;	
						}
			
		}

	}
 	
	public void MCS(double T, double[] H, Random flip, double ratio, String dynamics)
	{
	    int j=0;
	    
	    int rationalMCS= (int) (ratio*M);
	    for (int f=0; f< rationalMCS; f++)
	    {
		   j=(int) (flip.nextDouble()*M);
		   if(dynamics=="Metropolis")
			   MetropolisSpinflip(j, flip, T, H);
		   if(dynamics=="Glauber")
			   GlauberSpinflip(j, flip, T, H);
			   
	    }
	    Display();

	}
	
 	
  	public int TotalSpin()
 	{
 		int total=0;
 		for(int k=0;k<M;k++)
 		{
 			total+=spin[k];
 		}
 		return total;
 	}
  	
  	public int TotalSpinX()
 	{
 		int total=0;
 		for(int k=0;k<M;k++)
 		{
 			total+=(spin[k]*Xsign(k));
 		}
 		return total;
 	}
  	
  	public int TotalSpinY()
 	{
 		int total=0;
 		for(int k=0;k<M;k++)
 		{
 			total+=(spin[k]*Ysign(k));
 		}
 		return total;
 	}
  	
  	public int Xsign(int j)
  	{
  		int jx=j/L2;
  		int signX=1-2*(jx%2);
  		return signX;
  	}
  	
  	public int Ysign(int j)
  	{
  		int jy=j%L2;
  		int signY=1-2*(jy%2);
  		return signY;
  	}
  	
 	
 	public double SpinIntEnergy(int j)
 	{
 		double SpinE=0;
 		if(spin[j]!=0)
 		{
 				int b,kn,nn;
 			    for(b=0; b<4;b++)
 			    {
 				kn=Nneighber(b,j);
 				nn=NextNneighber(b,j);
 				SpinE+=(J1*spin[j]*spin[kn]+J2*spin[j]*spin[nn]);
 			    } 
 		
 		}
 		
 		return SpinE;
 	}
 	
 	public double TotalIntEnergy()
 	{
 		double TotalE=0;
 		for(int k=0; k<M;k++)
 			TotalE+=SpinIntEnergy(k);
 		totalintenergy=TotalE;
 		return TotalE;
 	}
 	
  	public double TotalEnergy(double field)
 	{
 	    return TotalIntEnergy()+field*totalspin;	
 	}
  	
  	public double Magnetization()
  	{
  		double m;
  		m=((double)totalspin/M);
  		return m;
  	}
  	
  	public double Mx()
  	{
  		double m;
  		m=((double)totalspinX/M);
        return m;
  	}
  	
  	public double My()
  	{
  		double m;
  		m=((double)totalspinY/M);
        return m;
  	}
  	
  	public double MM2()
  	{
  		double m;
  		m=mx*mx+my*my;
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