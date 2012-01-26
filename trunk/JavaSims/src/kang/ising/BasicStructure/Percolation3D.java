package kang.ising.BasicStructure;

import chris.util.Random;
import kang.ising.BasicStructure.BasicTools;

/*the way to use this class is 
 * 
 * 1  construct Percolation object percolation
 * 2  input ising strucuture of the critical droplet configuration
 * 3  use  percolation.probability();
 * 4  use  percolation.Mapping();
 * 5  the set of clusters are stored in percolation.CS;
 * 6  the order in size of the clusters can be found in percolation.CS.order[];
 */


public class Percolation3D{
	
	public IsingStructure3D ISP;
	public ClusterSet3D CS;
	public Cluster3D cltemp;
	public double probability;
	
	
	public int clustersize[];
	public int clusterindex;
	public double meanclustersize;
	public double SDclustersize;
	public int totalclusters;   //the number of the total clusters
	public double totalsites;   // the number of the total sites for this percolation problem
	public double OP;  // order parameter
	
	public boolean span;
	public int spannumber; //the number of dimension in which the cluster spans
	
	
	public Percolation3D()
	{
		
	}

	public Percolation3D(IsingStructure3D IS, int number)
	{
		this.ISP= IS.clone();
		this.CS= new ClusterSet3D(number);
		clustersize=new int[IS.M];
		for(int jj=0; jj<IS.M; jj++)
		{
			clustersize[jj]=-1;   //preset the clustersize array
		}
	}

	public void SetProbability(double p)
	{
		probability=p;
	}
	
	public void probability(double T)     // the coniglio-klein pb near spinodal
	{
		double percolationM=ISP.magnetization;
		probability=1- Math.exp(-(-ISP.J*(1+(percolationM/(1-ISP.percent)))/T));
	}
	

	
	public void fastNNMapping(int Pseed)   // for nearest neighbor system with R=0
	{
		int stablespin[]=new int[ISP.M];
		int copymap[]=new int[ISP.M];    // backup copy of the spin configuration, -1 is the stable spin, 0 is the dilution, -3 is the label if the site has already been chosen in the previous cluster
		int a0[]=new int[ISP.M];    //a0[j]=1 (there is bond between the up neighbor for site j and site j), otherwise a0[j]=0
		int a1[]=new int[ISP.M];    //a1[j]=1 (there is bond between the right neighbor for site j and site j), otherwise a0[j]=0
		int a2[]=new int[ISP.M];    //a2[j]=1 (there is bond between the down neighbor for site j and site j), otherwise a0[j]=0
		int a3[]=new int[ISP.M];    //a3[j]=1 (there is bond between the left neighbor for site j and site j), otherwise a0[j]=0
		int a4[]=new int[ISP.M];    //a3[j]=1 (there is bond between the up in z neighbor for site j and site j), otherwise a0[j]=0
		int a5[]=new int[ISP.M];    //a3[j]=1 (there is bond between the down in z neighbor for site j and site j), otherwise a0[j]=0
		
		int SN=0;
		// first, find all the stable spins
		for(int q=0; q<ISP.M; q++)
		{
			a0[q]=0;
			a1[q]=0;   //initialize the bond configuration to be null
			a2[q]=0;
			a3[q]=0;
			a4[q]=0;
			a5[q]=0;
			
			
			copymap[q]=ISP.spin[q];    //back up the lattice site's visiting history
			if(ISP.spin[q]==-1)
			{
				stablespin[SN]=q;
				SN++;
			}
		}
		
		// second, determine all the possible bonds and restore to a0[ISP.M] and a1[ISP.M]

		totalsites=SN;
		Random bondsrand;
		int pin;
		bondsrand = new Random(Pseed);
		//boolean nearestneighbor;
		int totalbonds=0;
		
		for (int k=0; k<SN; k++)
			for(int a=0; a<6; a++)  //consider 6 neighbors with double counting
			{
				int j=ISP.Nneighber(a, stablespin[k]);
				if(ISP.spin[j]==-1)
				{
					if(bondsrand.nextDouble()<=probability)		// third, throw random numbers on all the possible bonds to determine the percolation bonds
				   {
						if(a==0)
							{
							   a0[stablespin[k]]=1;
							   totalbonds++;
							}
						if(a==1)
							{
							   a1[stablespin[k]]=1;
							   totalbonds++;
							}
						if(a==2)
					    {
						   a2[stablespin[k]]=1;
						   totalbonds++;
						}
					    if(a==3)
						{
						   a3[stablespin[k]]=1;
						   totalbonds++;
						}
					    if(a==4)
						{
						   a4[stablespin[k]]=1;
						   totalbonds++;
						}
					    if(a==5)
						{
						   a5[stablespin[k]]=1;
						   totalbonds++;
						}
				   }
				}
			}
		
		int temp[]= new int [SN];     //store the locations of the sites belonging to current cluster
		int neighbor=0;
		for(int tp=0; tp<SN; tp++)
		{
			temp[tp]=-2;
		}
		// finally, find the largest cluster and save it into the array largest[]
		
		if(totalbonds>0)
		{
		   clusterindex=0;
           for (int i2=0; i2<SN; i2++)     //the loop of finding different clusters
           {
        	   if(copymap[stablespin[i2]]==-1)    //the site stablespin[i2] has not been visited
        	   {
        		   for(int tpi=0; tpi<SN; tpi++)           // everytime reset temp[]
        			   temp[tpi]=-2;
        		   pin=0;
        		   temp[pin]=stablespin[i2];        //write stablespin[i2] into the current cluster as the starting site
        		   copymap[stablespin[i2]]=-3;     //delete this site from the label map since it has been visited  
        		   pin++;
        		   if(a0[stablespin[i2]]==1)     //now look at its up neighbor
        		   {
        	           neighbor=ISP.Nneighber(0,stablespin[i2]);
        	           if(copymap[neighbor]!=-3)
        		       {
        	        	   temp[pin]=neighbor;    //write the up neighbor into cluster
        	        	   copymap[neighbor]=-3;
            		       pin++;
        		       }       		  
        		   }
        		   if(a1[stablespin[i2]]==1)     //now look at its right neighbor
    			   {
    			       neighbor=ISP.Nneighber(1,stablespin[i2]); 
    			       if(copymap[neighbor]!=-3)
        		       {
        	        	   temp[pin]=neighbor;    //write the right neighbor into cluster
        	        	   copymap[neighbor]=-3;
            		       pin++;
        		       }
    			   }
        		   if(a2[stablespin[i2]]==1)     //now look at its down neighbor
        		   {
        	           neighbor=ISP.Nneighber(2,stablespin[i2]); 
        	           if(copymap[neighbor]!=-3)
        		       {
        	        	   temp[pin]=neighbor;    //write the down neighbor into cluster
        	        	   copymap[neighbor]=-3;
            		       pin++;
        		       }
        		   }
        		   if(a3[stablespin[i2]]==1)     //now look at its left neighbor
    			   {
    			       neighbor=ISP.Nneighber(3,stablespin[i2]); 
    			       if(copymap[neighbor]!=-3)
        		       {
        	        	   temp[pin]=neighbor;    //write the left neighbor into cluster
        	        	   copymap[neighbor]=-3;
            		       pin++;
        		       }
    			   }
        		   if(a4[stablespin[i2]]==1)     //now look at its up in z neighbor
    			   {
    			       neighbor=ISP.Nneighber(4,stablespin[i2]); 
    			       if(copymap[neighbor]!=-3)
        		       {
        	        	   temp[pin]=neighbor;    //write the up in z neighbor into cluster
        	        	   copymap[neighbor]=-3;
            		       pin++;
        		       }
    			   }
        		   if(a5[stablespin[i2]]==1)     //now look at its down in z neighbor
    			   {
    			       neighbor=ISP.Nneighber(5,stablespin[i2]); 
    			       if(copymap[neighbor]!=-3)
        		       {
        	        	   temp[pin]=neighbor;    //write the down in z neighbor into cluster
        	        	   copymap[neighbor]=-3;
            		       pin++;
        		       }
    			   }
        		   
        		   
        		   
        		   
        		   for(int i3=0; i3< pin; i3++)       //the loop of scanning all the spins in temp[]
        		   {
            		   if(a0[temp[i3]]==1)     //now look at its up neighbor
            		   {
            			   neighbor=ISP.Nneighber(0,temp[i3]); 
            			   if(copymap[neighbor]!=-3)         //if this site is not in the cluster
            		       {
            				   temp[pin]=neighbor;    //write the up neighbor into cluster
            				   copymap[neighbor]=-3;
                		       pin++;
            		       }
            		      
            		   }
            		   if(a1[temp[i3]]==1)     //now look at its right neighbor
        			   {
        			       neighbor=ISP.Nneighber(1,temp[i3]); 
        			       if(copymap[neighbor]!=-3)         //if this site is not in the cluster
        			       {
        			    	   temp[pin]=neighbor;    //write the right neighbor into cluster
        			    	   copymap[neighbor]=-3;
            			       pin++;
        			       }
        			       
        			   }
            		   if(a2[temp[i3]]==1)     //now look at its down neighbor
        			   {
        			       neighbor=ISP.Nneighber(2,temp[i3]); 
        			       if(copymap[neighbor]!=-3)         //if this site is not in the cluster
        			       {
        			    	   temp[pin]=neighbor;    //write the down neighbor into cluster
        			    	   copymap[neighbor]=-3;
            			       pin++;
        			       }
        			       
        			   }
            		   if(a3[temp[i3]]==1)     //now look at its left neighbor
        			   {
        			       neighbor=ISP.Nneighber(3,temp[i3]); 
        			       if(copymap[neighbor]!=-3)         //if this site is not in the cluster
        			       {
        			    	   temp[pin]=neighbor;    //write the left neighbor into cluster
        			    	   copymap[neighbor]=-3;
            			       pin++;
        			       }
        			       
        			   }
            		   
            		   if(a4[temp[i3]]==1)     //now look at its up in z neighbor
            		   {
            			   neighbor=ISP.Nneighber(4,temp[i3]); 
            			   if(copymap[neighbor]!=-3)         //if this site is not in the cluster
            		       {
            				   temp[pin]=neighbor;    //write the up in z neighbor into cluster
            				   copymap[neighbor]=-3;
                		       pin++;
            		       }
            		      
            		   }
            		   if(a5[temp[i3]]==1)     //now look at its down in z neighbor
        			   {
        			       neighbor=ISP.Nneighber(5,temp[i3]); 
        			       if(copymap[neighbor]!=-3)         //if this site is not in the cluster
        			       {
        			    	   temp[pin]=neighbor;    //write the down in z neighbor into cluster
        			    	   copymap[neighbor]=-3;
            			       pin++;
        			       }
        			       
        			   }
            		   
        		   }
     
        		   int copytemp[]= new int[ISP.M];   //we need to covert temp[SN]->copytemp[ISP.M] in order to match the cluster generating function
        		   for(int c=0;c<ISP.M;c++)
        		   {
        			   copytemp[c]=-2;        //initialize the cluster map
        		   }
        		   for(int ct=0;ct<SN;ct++)
        		   {
        			   copytemp[ct]=temp[ct];
        		   }
        		   cltemp=new Cluster3D(ISP.L1,ISP.L2,ISP.L3,copytemp);
        		   clustersize[clusterindex]=cltemp.size;
        		   clusterindex++;
        		   CS.AddCluster(cltemp);   
        	   }
           } // the end of i2 cycle
           CS.ordering();
           CS.findmaximum();
           CS.findminimum();
		}

	}
	
	public void SpanningCluster()       //if want to check the largest cluster then, map[]==CS.set[order[0]].lattice[]
	{
		//the logic of determine if the cluster spans the whole lattice:
		//1 the coordinate x,y or z covers 0 to L-1
		//2 the site of cluster on edge has neighbor on the other edge
		//if x covers 0 to L-1, then spanX=true, then check on the x=0 and x=L-1 plane, if there is point y,z on both plane to determine spanXpbc
		
		int L= ISP.L1;
	    span=false;
		boolean spanX=true;
		boolean spanY=true;
		boolean spanZ=true;
		boolean spanXpbc=false;     //the cluster wrap back to itself under periodic boundary condition              
		boolean spanYpbc=false;
		boolean spanZpbc=false;
		
		//int plane[]= new int[L*L];
		
		spannumber=0;
		
		int X[]=new int[L];
		int Y[]=new int[L];
		int Z[]=new int[L];
		int x, y, z;
		//int x1, y1, xL, yL;   //the coordinate j of the site in cluster on the edge
		
		for(int i=0; i<L; i++)
		{
			X[i]=1;
			Y[i]=1;
			Z[i]=1;
			
		}// preset the array to count the x and y coordinates
		for(int j=0; j<(L*L*L); j++)
		{
			if(CS.set[CS.maximumpin].lattice[j]==2)
			{
				
				x=(int)(j/(L*L));
				y=(int)((int)j%(L*L))/L;
				z=(int)((int)j%(L*L))%L;
				
				X[x]=0;
				Y[y]=0;
				Z[z]=0;
		
				
			}
		}// scan the whole cluster distribution to mark all the occupied x and y coordinates 
		for(int k=0; k<L; k++)
		{
			if(X[k]!=0)
				spanX=false;
			if(Y[k]!=0)
				spanY=false;
			if(Z[k]!=0)
				spanZ=false;
		}

		if(spanX)
		{
			for(int ty=0; ty<L; ty++)
				for(int tz=0; tz<L; tz++)
				{
					if((CS.set[CS.maximumpin].lattice[ty*L+tz]+CS.set[CS.maximumpin].lattice[(L-1)*L*L+ty*L+tz])==4)
						spanXpbc=true;
				}
			
			if(spanXpbc)
			{
				span=true;
				spannumber++;
			}
			
		}
		if(spanY)
		{
			for(int tx=0; tx<L; tx++)
				for(int tz=0; tz<L; tz++)
				{
					if((CS.set[CS.maximumpin].lattice[tx*L*L+tz]+CS.set[CS.maximumpin].lattice[tx*L*L+L*(L-1)+tz])==4)
						spanYpbc=true;
				}
			
			if(spanYpbc)
			{
				span=true;
				spannumber++;
			}
			
		}
		if(spanZ)
		{
			for(int ty=0; ty<L; ty++)
				for(int tx=0; tx<L; tx++)
				{
					if((CS.set[CS.maximumpin].lattice[tx*L*L+ty*L]+CS.set[CS.maximumpin].lattice[tx*L*L+ty*L+L-1])==4)
						spanZpbc=true;
				}
			
			if(spanZpbc)
			{
				span=true;
				spannumber++;
			}
			
		}
		
		 
	
	}
	
	public void ClusterSize()
	{
		double totalsize=0;
		totalclusters=0;
	    
		for(int i=0; clustersize[i]>=0; i++)
		{
			totalsize+=clustersize[i];
			totalclusters++;
		}
		if(span)   //get rid of the largest cluster
		{
			totalsize-=CS.maximumsize;
		}
		
		
		if(totalclusters>0)
		{
			if(span)
				meanclustersize=totalsize/(totalclusters-1);
			else
				meanclustersize=totalsize/totalclusters;
		}
			
		else
			meanclustersize=-1;    //this indicates that there is no cluster at all
	}
	
	public void SDClusterSize()
	{
		double totalSD=0;
		double SD=0;
		if(totalclusters>0)
		{
			for(int p=0; p<totalclusters;p++)
			{
				totalSD+=((clustersize[p]-meanclustersize)*(clustersize[p]-meanclustersize));
			}
			if(span)
			{
				totalSD-=((CS.maximumsize-meanclustersize)*(CS.maximumsize-meanclustersize));
			}
			
			if(span)
				SD=Math.sqrt(totalSD/(totalclusters-1));
			else
				SD=Math.sqrt(totalSD/totalclusters);
		}
		else
			SD=-1;  // this indicates that there is no cluster at all
		
		SDclustersize=SD;
	}
	
    public void OrderParameter()
    {
    	if(totalsites>0)
    		OP=((double)CS.maximumsize)/totalsites;
    	else
    		OP=-1;   // this indicates that there is no stable spin sites
    }
	
}