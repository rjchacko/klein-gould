package kang.ising.BasicStructure;



import chris.util.Random;
//import chris.util.PrintUtil;

public class IsingStructure3D{ 
	
	//structures and parameters
	public int spin[];
	public int initialcopy[];   //the array of the initial copy of the system
	public int biaslabel[];
	public double dilutionmap[];
	public int largestcluster[];
	
	
	public int L1, L2, L3, M; //parameters for the lattice                                                                                                                                                                                                                                                                                                                                                                                                        
	public double J;     //interaction constant after normalization
	public double NJ;    //interaction constant before normalization
	public double percent;  //dilution percent
	public int deadsites;  //the number of diluted sites in the system
	public String shape;
	
	public double totalintenergy;
	public int totalspin;
	public double magnetization;
	
	public int R;   //interaction range R=0 is NN interaction
	
	public int biasA;
	public int biasB;
	public int biasC;
	public double biaspercent;
	
	//the function for this class IsingStructure3D
	
	public IsingStructure3D(int L1, int L2, int L3, int R, double NJ, double percent, double biaspercent, String shape)     //generating function
	{
		this.L1=L1;
		this.L2=L2;
		this.L3=L3;
		this.M=L1*L2*L3;
		this.R=R;
		this.percent=percent;
		this.biaspercent=biaspercent;
		this.shape=shape;

		this.NJ=NJ;
		if(R==0)
			this.J=NJ/6;
		if(R>0)
			{
			if(shape=="square")
				this.J=NJ/((2*R+1)*(2*R+1)*(2*R+1)-1);
			}
				
		this.spin=new int [M];
		this.initialcopy=new int [M];
		this.biaslabel=new int [M];
		this.dilutionmap=new double[M];
		
	}
	
	public IsingStructure3D clone()
	{
		IsingStructure3D copy= new IsingStructure3D(L1, L2, L3, R, NJ, percent, biaspercent, shape);
		for(int t=0;t<M; t++)
		{
			copy.spin[t]=spin[t];
			copy.initialcopy[t]=initialcopy[t];
			copy.biaslabel[t]=biaslabel[t];
			copy.dilutionmap[t]=dilutionmap[t];
		}
		copy.biasA=biasA;
		copy.biasB=biasB;
		copy.biasC=biasC;
		copy.totalspin=totalspin;
		copy.totalintenergy=totalintenergy;
		
		
		return copy;
	}
	
	public void Dinitialization(int Dseed, int Bseed, int A, int B, int C )//dilution initialization
	{
		Random Drand= new Random(Dseed);
		Random Brand= new Random(Bseed);
		biasA=A;
		biasB=B;
		biasC=C;
		deadsites=0;

		
		for (int t=0;t<M;t++)
			{
			spin[t]=1;
			biaslabel[t]=0;
			dilutionmap[t]=0;
			}
			
		
		int cx, cy, cz; //x,y,z index for the center
		
		cx=(int)L1/2;
		cy=(int)L2/2;
		cz=(int)L3/2;
		if(biaspercent!=percent)
		    cube(biaslabel, A,B, C, cx,cy, cz);      // need to wrtie a function named cube which draws a cube centered at cx,cy,cz                      
		
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
	
	public void dilutionmap(int range)
	{
		for(int k=0; k<M; k++)
		{
			dilutionmap[k]=dilutionratio(R,k);
		}
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
		
		totalintenergy= TotalIntEnergy();
		magnetization=Magnetization();
		
	}
	
//basic functions
	
	
	
	public double dilutionratio(int r, int i)    //calculate the dilution ratio within the circle with radius r
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
	}
	
	public int getX(int i)
	{
		int X=0;
		X=(int)(i/(L2*L3));
		return X;

	}
	
	public int getY(int i)
	{
		int Y=0;
		Y=(int)((int)i%(L2*L3))/L3;
		return Y;

	}
	
	public int getZ(int i)
	{
		int Z=0;
		Z=(int)((int)i%(L2*L3))%L3;
		return Z;

	}
	
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
	
	public int Z(int bz)
	{
		int realz=bz;
		if (bz>=L3)
			realz=bz-L3;
		if (bz<0)
			realz=bz+L3;
		return realz;
	}
	
	public int distance (int a, int b)     // the code to calculate the square of the distance between two points on the lattice
	{
		int dis=0;
		int ax, ay, az, bx, by, bz;
		int dx2, dy2, dz2;
		ax= a/L2;
		ay= a%L2;
		ax=(int)(a/(L2*L3));
		ay=(int)((int)a%(L2*L3))/L3;
		az=(int)((int)a%(L2*L3))%L3;
		
		bx=(int)(b/(L2*L3));
		by=(int)((int)b%(L2*L3))/L3;
		bz=(int)((int)b%(L2*L3))%L3;
		
		dx2=(ax-bx)*(ax-bx);
		dy2=(ay-by)*(ay-by);
		dz2=(az-bz)*(az-bz);
		if((ax-bx+L1)*(ax-bx+L1)<(ax-bx)*(ax-bx))
			dx2=(ax-bx+L1)*(ax-bx+L1);
		if((ax-bx-L1)*(ax-bx-L1)<(ax-bx)*(ax-bx))
			dx2=(ax-bx-L1)*(ax-bx-L1);
		if((ay-by+L2)*(ay-by+L2)<(ay-by)*(ay-by))
			dy2=(ay-by+L2)*(ay-by+L2);
		if((ay-by-L2)*(ay-by-L2)<(ay-by)*(ay-by))
			dy2=(ay-by-L2)*(ay-by-L2);
		if((az-bz+L3)*(az-bz+L3)<(az-bz)*(az-bz))
			dz2=(az-bz+L3)*(az-bz+L3);
		if((az-bz-L3)*(az-bz-L3)<(az-bz)*(az-bz))
			dz2=(az-bz-L3)*(az-bz-L3);

		dis=dx2+dy2+dz2;
		return dis;
	}
	
	public int SumInRange(int spin[], int j, int R, String shape)  //sum in the range of (2R+1)*(2R+1) square or (R+1)^2-R^2 diamond
	{
		int S=0;
		int nx=(int)(j/(L2*L3));
		int ny=(int)((int)j%(L2*L3))/L3;
		int nz=(int)((int)j%(L2*L3))%L3;
		int kx, ky, kz;

		if(shape=="square")
		{
			for (int m=-R; m<=R; m++)
				for (int n=-R; n<=R; n++)
					for(int l=-R; l<=R; l++)
			    {
				   kx=X(nx+m);
				   ky=Y(ny+n);
				   kz=Z(nz+l);
				   S+=spin[kx*L2*L3+ky*L3+kz];	
			    }
		}

		
		return S;
	}
	
	public void cube(int label[], int a, int b,int c, int cx, int cy, int cz)  //draw a cube of 2a*2b*2c at (cx,cy,cz)
	{
		int bx, by, bz;
		int x,y,z;
		
		for(bx=cx-a; bx<cx+a; bx++)
			for(by=cy-b; by<cy+b; by++)
				for(bz=cz-c; bz<cz+c; bz++)
			{
				x=X(bx);
				y=Y(by);
				z=Z(bz);
				label[x*L2*L3+y*L3+z]=1;
			}
	}
	
	public int Nneighber(int a,int i ){// function for the index of nearest neighbor
		int nx,ny,nz; //index for neighbor
		int ni=0;
		nx=(int)(i/(L2*L3));
		ny=(int)((int)i%(L2*L3))/L3;
		nz=(int)((int)i%(L2*L3))%L3;
		
		if (a==0) {
			ni=nx*L2*L3+(ny-1)*L3+nz;
			if  (ny==0) {
				ni=nx*L2*L3+(ny+L2-1)*L3+nz;
		}
			
		}//(x,y-1,z) up
		
     	if (a==1){
			ni=(nx+1)*L2*L3+ny*L3+nz;
			if  (nx==L1-1) {
				ni=(nx-L1+1)*L2*L3+ny*L3+nz;
			}
			
		}//(x+1,y,z) right
		
		if (a==2){
			ni=nx*L2*L3+(ny+1)*L3+nz;
			if  (ny==L2-1) {
				ni=nx*L2*L3+(ny-L2+1)*L3+nz;
			}
			
		}//(x,y+1,z) down
		
		if (a==3){
			ni=(nx-1)*L2*L3+ny*L3+nz;
			if  (nx==0) {
				ni=(nx+L1-1)*L2*L3+ny*L3+nz;
			}
		}//(x-1,y,z) left
		
		if (a==4){
			ni=nx*L2*L3+ny*L3+nz+1;
			if  (nz==L3-1) {
				ni=nx*L2*L3+ny*L3+nz-L3+1;
			}
		}//(x,y,z+1) up in z
		
		if (a==5){
			ni=nx*L2*L3+ny*L3+nz-1;
			if  (nz==0) {
				ni=nx*L2*L3+ny*L3+nz+L3-1;
			}
		}//(x,y,z-1) down in z
		
		
		return ni;
		
	}
	
 	public double interactionEchange (int j)//function for interaction energy
 	{ 
		double Energy=0;
		double Energychange=0;
		if (R==0)
		{
			int b,k;
		    for(b=0; b<6;b++)
		    {
			k=Nneighber(b,j);
			Energy=Energy+J*spin[j]*spin[k];
		    }
		    Energychange=-2*Energy;
		    
		}
		if (R!=0)
		{
			int S=SumInRange(spin, j, R, shape);
			
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
				//totalintenergy+=IntEchange;
				
			}
			
			else
			{
				if(flip.nextDouble()<=Math.exp(-Echange/temperature))
						{
					            spin[j]=-spin[j];
								totalspin-=tempspin*2;
								magnetization=Magnetization();
								//totalintenergy+=IntEchange;	
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
	    
	    int rationalMCS= (int) (ratio*M);
	    for (int f=0; f< rationalMCS; f++)
	    {
		   j=(int) (flip.nextDouble()*M);
		   if(dynamics=="Metropolis")
			   MetropolisSpinflip(j, flip, T, H);
		   if(dynamics=="Glauber")
			   GlauberSpinflip(j, flip, T, H);
			   
	    }

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
 	
 	public double SpinIntEnergy(int j)
 	{
 		double SpinE=0;
 		if(spin[j]!=0)
 		{
 			if (R==0)
 			{
 
 				int b,k;
 			    for(b=0; b<6;b++)
 			    {
 				k=Nneighber(b,j);
 				SpinE+=J*spin[j]*spin[k];
 			    } 
 			}
 			
 			if (R!=0)
 			{
 				int S=SumInRange(spin, j, R, shape);

 				SpinE=J*spin[j]*S-J;
 				
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
 	
 	
    public int[] Clustergrowth(int spin[], int direction, double pb, int Sseed, int Rseed, boolean keeplargest)  // given the spin configuration and stable direction, generate the cluster size distribution data use cluster growth
    {
    	largestcluster=new int[M];
    	int Max=0;
    	int[] cltemp;   // store the location on the lattice of the current cluster
    	int[] clustersize= new int[M];
    	int[] ul=new int[M];      //usefullocation array store the spin's location on the lattice
    	int[] im=new int[M];      //indexmap stores the index of a spin in ul[] array
    	int ulpin=0;    // the index of the last non zero element in ul[] array, ulpin+1=total number of useful spins
    	
    	//now do the initalization of ul[] and im[] also clustersize[]
    	for(int i=0; i<M; i++)
    	{
    		ul[i]=-1;
    		im[i]=-1;
    		clustersize[i]=-1;
    	}
    	
    	//and take spin[] to make the real ul and im
    	for(int i=0;i<M;i++)
    	{
 
    		if(spin[i]==direction)
    		{
    			ul[ulpin]=i;
    			im[i]=ulpin;
    			ulpin++;   //now ulpin point to the first -1, ulpin= total
    		}
    	}
    	ulpin--;  //now ulpin point to the last spin location, so uplin+1=total
    	
    	
    	
    	////////////////////////////now end of the initialization/////////////////////////
    	Random srand=new Random(Sseed);
    	Random rrand=new Random(Rseed);
    	
    	int sizepin=0;
    	
    	while(ulpin>0)  // this loop is for the seed of each cluster, the seed site is randomly chosen from the left over useful sites
    	{
    		int totalleft=ulpin+1;  //total left over useful spin
    		cltemp=new int[totalleft];
    		int clpin=1;
    		for(int c=0; c<totalleft; c++)
    		{
    			cltemp[c]=-1;
    		}
    		
    		int j=(int)(srand.nextDouble()*(totalleft));    //pick the seed site from ul[]
    		//use j as the seed for the cluster, plant the seed
    		{
    			cltemp[0]=ul[j];    //write ul[j] into cltemp
    			im[ul[j]]=-1;       //delete j from indexmap
    			im[ul[ulpin]]=j;    //swap in indexmap
    			ul[j]=ul[ulpin];    //swap in ul[]
    			ulpin--;            //delete j from ul[]
    		}
    		
    		//now generate the cluster from this seed
            int k=0;
    		
    		if(ulpin>=0)
    		{
        		for(k=0; cltemp[k]>=0; k++)
        		{
        			int kx=getX(cltemp[k]);
        			int ky=getY(cltemp[k]);
        			int kz=getZ(cltemp[k]);
        			
        			if(R==0)       //for nearest neighbor case
        			{
        				for(int a=0; a<6; a++)
        				{
        					int s=Nneighber(a, cltemp[k]);
        					if(im[s]>=0)
        						if(rrand.nextDouble()<=pb)
        						{
        							cltemp[clpin]=s;  //add s into cltemp
        							clpin++;  //increase the length of cltemp, alway point to the first -1
        							im[ul[ulpin]]=im[s];   //swap in indexmap
        							ul[im[s]]=ul[ulpin];  //swap in ul[]			
        							im[s]=-1;  //delete from indexmap
        							ulpin--;   //delete from ul[]
        							}
        				}
        			}
        			
        			
        			if(R>0)
        			{
        				for(int m=-R; m<=R; m++)
        				for(int n=-R; n<=R; n++)
        				for(int l=-R; l<=R; l++)
        				{
        					int sx=X(kx+m);
        					int sy=Y(ky+n);
        					int sz=Z(kz+l);
  
        					int s=sx*L2*L3+sy*L3+sz;
        					if((im[s]>=0)&&(s!=cltemp[k]))
        						if(rrand.nextDouble()<=pb)
        						{
        							cltemp[clpin]=s;  //add s into cltemp
        							clpin++;  //increase the length of cltemp, alway point to the first -1
        							im[ul[ulpin]]=im[s];   //swap in indexmap
        							ul[im[s]]=ul[ulpin];  //swap in ul[]			
        							im[s]=-1;  //delete from indexmap
        							ulpin--;   //delete from ul[]
        							}
        				}
        		    }
        		}
    		}
    		
    		else
    			k=1;
    		
    		
    		if(k>Max)
    		{
    			Max=k;
    			if(keeplargest)
    				{
    				for(int l=0; l<M; l++)
    				{
    					largestcluster[l]=-1;
    				}
    				for(int l=0; l<totalleft; l++)
    				{
    					if(cltemp[l]>=0)
    						largestcluster[cltemp[l]]=2;
    				}
    				}
    					

    		}
    		clustersize[sizepin]=k;
    		sizepin++;
    		
    		
    	}
    	
    	
    	if(ulpin==0)
    	{
    		clustersize[sizepin]=1;      //the last site left over
    		sizepin++;
    	}
    	
    	
    	return clustersize;
    }
    
    
    
 	
}

	