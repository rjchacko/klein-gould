package kang.ising.BasicStructure;


import kang.ising.BasicStructure.StructureFactor;
import chris.util.Random;
//import chris.util.PrintUtil;

public class IsingStructure{
	
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
	
	
	
	public int L1, L2, M; //parameters for the lattice                                                                                                                                                                                                                                                                                                                                                                                                        
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
	public double biaspercent;
	
	//the function for this class IsingStructure
	
	public IsingStructure(int L1, int L2, int R, double NJ, double percent, double biaspercent, String shape)     //generating function
	{
		this.L1=L1;
		this.L2=L2;
		this.M=L1*L2;
		this.R=R;
		this.percent=percent;
		this.biaspercent=biaspercent;
		this.shape=shape;
		this.deadsites=0;

		this.NJ=NJ;
		if(R==0)
			this.J=NJ/4;
		if(R>0)
			{
			if(shape=="square")
				this.J=NJ/((2*R+1)*(2*R+1)-1);
			if(shape=="diamond")
				this.J=NJ/((R+1)*(R+1)+R*R-1);
			}
				
		this.spin=new int [M];
		this.initialcopy=new int [M];
		this.biaslabel=new int [M];
		this.dilutionmap=new double[M];
		this.damage=new double[M];
		
	}
	
	public IsingStructure clone()
	{
		IsingStructure copy= new IsingStructure(L1, L2, R, NJ, percent, biaspercent, shape);
		for(int t=0;t<M; t++)
		{
			copy.spin[t]=spin[t];
			copy.initialcopy[t]=initialcopy[t];
			copy.biaslabel[t]=biaslabel[t];
			copy.dilutionmap[t]=dilutionmap[t];
			copy.damage[t]=damage[t];

		}
		copy.biasA=biasA;
		copy.biasB=biasB;
		copy.totalspin=totalspin;
		copy.totalintenergy=totalintenergy;
		copy.deadsites=deadsites;
		
		
		return copy;
	}
	
	public IsingStructure Dperturbation(int seed)
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
				perturb.damage[j]=1;
				perturb.deadsites++;
				i--;
			}
		}
		perturb.Sinitialization(1, 1);
		
		return perturb;
	}
	
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
			dilutionmap[t]=0;
			damage[t]=0;
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
		
		/*or(int k=0; k<M; k++)
		{
			dilutionmap[k]=dilutionratio(R,k);
		} // draw the dilution map
	*/
	}
	
	public void dilutionmap(int range, String shape)
	{
		for(int k=0; k<M; k++)
		{
			if(shape=="Square")
				dilutionmap[k]=dilutionratioSquare(R,k);
			if(shape=="Circle")
				dilutionmap[k]=dilutionratioCircle(R,k);
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
				   {
				   spin[i]=1;
				   
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
	
	public void SpinSetup()
	{
		upspin=new double[M];
		downspin=new double[M];
		for(int j=0; j<M;j++)
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
	
	public double dilutionratioSquare(int r, int i)
	{
		double ratio;
		double dilutedsite;
		dilutedsite=0;
		double totalinrange;
		totalinrange=0;
		
		int nx=i/L2;
		int ny=i%L2;
		int kx, ky;
			
		for (int m=-R; m<=R; m++)
			for (int n=-R; n<=R; n++)
			    {
				   kx=X(nx+m);
				   ky=Y(ny+n);
				   totalinrange++;
				   if(spin[kx*L2+ky]==0)
					   dilutedsite++;			   	
			    }
	
		ratio=dilutedsite/totalinrange;
		return ratio;
	}
	
	public double dilutionratioCircle(int r, int i)    //calculate the dilution ratio within the circle with radius r
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
 			    for(b=0; b<4;b++)
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
 	
 	public double[] LocalEnergy(double field)
 	{
 		double localenergy[]=new double[M];
 		for(int le=0; le<M; le++)
 		{
 			if(spin[le]!=0)
 			{
 				localenergy[le]=SpinIntEnergy(le)-field*spin[le];
 			}
 			else 
 			    localenergy[le]=0;
 		}
 		
 		return localenergy;
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
 	    return totalintenergy+field*totalspin;	
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
  	
    public void SpinSF()
    {
    	SpinSetup();
    	SFup=new StructureFactor(L1, (double)L1);  //this only works if L1==L2
    	SFdown=new StructureFactor(L1, (double)L1);  //this only works if L1==L2
    	SFup.takeFT(upspin);
    	SFdown.takeFT(downspin);
    	SFup.SquareAverage();
    	SFdown.SquareAverage();
    }
    
    public void DilutionSF()
    {
    	
    	SFdilution=new StructureFactor(L1, (double)L1);  //this only works if L1==L2
    	SFdilution.takeFT(damage);
    	
    }
   
    public int[] ClusterInfo(int cluster[])    // the function to calculate the center position and the size of it and store as: clusterinfo[0]=size clusterinfo[1]=center
	{
		int clusterinfo[]=new int[2];
    	int sites[]= new int[cluster.length];
		int pin=0;
    	for(int ii=0; ii< cluster.length; ii++)
		{
			if(cluster[ii]==2)
				{
				sites[pin]=ii;
				pin++;
				}
		}
    	int size=pin;
    	
    	int fp,fx,fy;
		fp=sites[0];
		fx=fp/L2;
		fy=fp%L2;
		double sumx,sumy;
		sumx=0;
		sumy=0;
		int cx,cy;
		int clx[];
		int cly[];
		int tempx,tempy;
		clx= new int[size];
		cly= new int[size];
		clx[0]=fx;
		cly[0]=fy;
		sumx=fx;
		sumy=fy;
		
		for(int cl=1; cl<size; cl++)
		{
			tempx=sites[cl]/L2;
			tempy=sites[cl]%L2;
			clx[cl]=tempx;
			cly[cl]=tempy;
			if((tempx-fx)*(tempx-fx)>(tempx-L1-fx)*(tempx-L1-fx))
				clx[cl]=tempx-L1;
			if((tempx-fx)*(tempx-fx)>(tempx+L1-fx)*(tempx+L1-fx))
				clx[cl]=tempx+L1;
			if((tempy-fy)*(tempy-fy)>(tempy-L2-fy)*(tempy-L2-fy))
				cly[cl]=tempy-L2;
			if((tempy-fy)*(tempy-fy)>(tempy+L2-fy)*(tempy+L2-fy))
				cly[cl]=tempy+L2;
			sumx+=clx[cl];
			sumy+=cly[cl];
		}
		cx=(int)(sumx/size);
		cy=(int)(sumy/size);
		
		clusterinfo[0]=size;
		clusterinfo[1]=X(cx)*L2+Y(cy);
		
		return clusterinfo;
	}
    
    public int[] LargestCluster(int spin[], int direction, double pb, int Sseed, int Rseed)
    {
    	int Max=0;
    	int cluster[]=new int[spin.length];
    	int[] cltemp;   // store the location on the lattice of the current cluster
    	int[] clustersize= new int[M];
    	int[] ul=new int[M];      //usefullocation array store the spin's location on the lattice
    	int[] im=new int[M];      //indexmap stores the index of a spin in ul[] array
    	int ulpin=0;    // the index of the last non negative element in ul[] array, ulpin+1=total number of useful spins
    	
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
    		cltemp=new int[M];
    		int clpin=1;
    		for(int c=0; c<M; c++)
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
        			int kx=cltemp[k]/L2;
        			int ky=cltemp[k]%L2;
        			
        			if(ulpin>0)
        			
        			if(R==0)       //for nearest neighbor case
        			{
        				for(int a=0; a<4; a++)
        				{
        					int s=Nneighber(a, cltemp[k]);
        					if(im[s]>=0)
        						if(rrand.nextDouble()<=pb)
        						{
        							cltemp[clpin]=s;  //add s into cltemp
        							clpin++;  //increase the length of cltemp, always point to the first -1
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
        				{
        					int sx=X(kx+m);
        					int sy=Y(ky+n);
        					int s=sx*L2+sy;
        					if((im[s]>=0)&&(s!=cltemp[k]))
        						if(rrand.nextDouble()<=pb)
        						{
        							cltemp[clpin]=s;  //add s into cltemp
        							clpin++;  //increase the length of cltemp, always point to the first -1
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
    			for(int ii=0; ii<spin.length; ii++)
    			{
    				cluster[ii]=-1;
    			}
    			
    			for(int jj=0; jj<spin.length; jj++)
    			{
    				if(cltemp[jj]>=0)
    					cluster[cltemp[jj]]=2;
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
    	
    	
    	return cluster;
    }
    
    public int LargestClusterSize(int spin[], int direction, double pb, int Sseed, int Rseed)   //a shorter version of code just to return the size of the largest cluster
    {
    	int Max=0;
    	int[] cltemp;   // store the location on the lattice of the current cluster
    	int[] clustersize= new int[M];
    	int[] ul=new int[M];      //usefullocation array store the spin's location on the lattice
    	int[] im=new int[M];      //indexmap stores the index of a spin in ul[] array
    	int ulpin=0;    // the index of the last non negative element in ul[] array, ulpin+1=total number of useful spins
    	
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
    		cltemp=new int[M];
    		int clpin=1;
    		for(int c=0; c<M; c++)
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
        			int kx=cltemp[k]/L2;
        			int ky=cltemp[k]%L2;
        			
        			if(ulpin>0)
        			
        			if(R==0)       //for nearest neighbor case
        			{
        				for(int a=0; a<4; a++)
        				{
        					int s=Nneighber(a, cltemp[k]);
        					if(im[s]>=0)
        						if(rrand.nextDouble()<=pb)
        						{
        							cltemp[clpin]=s;  //add s into cltemp
        							clpin++;  //increase the length of cltemp, always point to the first -1
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
        				{
        					int sx=X(kx+m);
        					int sy=Y(ky+n);
        					int s=sx*L2+sy;
        					if((im[s]>=0)&&(s!=cltemp[k]))
        						if(rrand.nextDouble()<=pb)
        						{
        							cltemp[clpin]=s;  //add s into cltemp
        							clpin++;  //increase the length of cltemp, always point to the first -1
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
    		}
    		clustersize[sizepin]=k;
    		sizepin++;
    	}
    	
    	if(ulpin==0)
    	{
    		clustersize[sizepin]=1;      //the last site left over
    		sizepin++;
    	}
    	
    	
    	return Max;
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
        			int kx=cltemp[k]/L2;
        			int ky=cltemp[k]%L2;
        			
        			if(R==0)       //for nearest neighbor case
        			{
        				for(int a=0; a<4; a++)
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
        				{
        					int sx=X(kx+m);
        					int sy=Y(ky+n);
        					int s=sx*L2+sy;
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