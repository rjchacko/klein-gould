package kang.ising.BasicStructure;

import chris.util.Random;


/*the way to use this class is 
 * 
 * 1  construct Percolation object percolation
 * 2  input ising strucuture of the critical droplet configuration
 * 3  use  percolation.probability();
 * 4  use  percolation.Mapping();
 * 5  the set of clusters are stored in percolation.CS;
 * 6  the order in size of the clusters can be found in percolation.CS.order[];
 */


public class Percolation{
	
	public IsingStructure ISP;
	public ClusterSet CS;
	public Cluster cltemp;
	public double probability;
	
	public int clustersize[];
	public int clusterindex;
	public double meanclustersize;
	public double SDclustersize;
	public int totalclusters;   //the number of the total clusters
	public double totalsites;   // the number of the total sites for this percolation problem
	public double OP;  // order parameter
	
	
	public Percolation()
	{
		
	}

	public Percolation(IsingStructure IS, int number)
	{
		this.ISP= IS.clone();
		this.CS= new ClusterSet(number);
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
	
	public void probability(double T)
	{
		double percolationM=ISP.magnetization;
		probability=1- Math.exp(-(-ISP.J*(1+(percolationM/(1-ISP.percent)))/T));
	}
	
	public void Mapping()    //for long range system with R as the interaction
	{
		int stablespin[]=new int[ISP.M];
		int SN=0;
		// first, find all the stable spins
		for(int q=0; q<ISP.M; q++)
		{
			if(ISP.spin[q]==-1)
			{
				stablespin[SN]=q;
				SN++;
			}	
		}
		
		// second, determine all the possible bonds
		int kx, ky, jx, jy;  // location for i,j
		int rx, ry;   // the difference between k and j
		int rx2,ry2;
		int PN=0;
		int PBsize=0;
		PBsize=(int)(SN*(2*ISP.R+1)*(2*ISP.R+1))/2;
		int RBH[]= new int [PBsize];
		int RBT[]= new int [PBsize];
		
		for(int pb=0; pb<PBsize; pb++)
		{
			RBH[pb]=-1;
			RBT[pb]=-1;
		}             //initialize all the possible bonds data structure
		
		Random bondsrand;
		int pin;
		bondsrand = new Random(1);
		
		
		for (int k=0; k<SN; k++)
			for(int j=0; j<k; j++)
			{
				kx= stablespin[k]/ISP.L2;
				ky= stablespin[k]%ISP.L2;
				jx= stablespin[j]/ISP.L2;
				jy= stablespin[j]%ISP.L2;
				rx= kx-jx;
				ry= ky-jy;
				rx2=rx*rx;
				if(rx2>=(rx-ISP.L1)*(rx-ISP.L1))
					rx2=(rx-ISP.L1)*(rx-ISP.L1);
				if(rx2>=(rx+ISP.L1)*(rx+ISP.L1))
					rx2=(rx+ISP.L1)*(rx+ISP.L1);
				
				ry2=ry*ry;
				if(ry2>=(ry-ISP.L2)*(ry-ISP.L2))
					ry2=(ry-ISP.L2)*(ry-ISP.L2);
				if(ry2>=(ry+ISP.L2)*(ry+ISP.L2))
					ry2=(ry+ISP.L2)*(ry+ISP.L2);
				
				
				if( (rx2<=ISP.R*ISP.R) & (ry2<=ISP.R*ISP.R) )
					if(bondsrand.nextDouble()<=probability)		// third, throw random numbers on all the possible bonds to determine the percolation bonds
				{
					RBH[PN]=stablespin[k];
					RBT[PN]=stablespin[j];
					PN++;
				}
				
			}
		
		int totalbonds=PN;
		int temp[]= new int [SN];
		int writeintotemp;
		for(int tp=0; tp<SN; tp++)
		{
			temp[tp]=-2;
		}
		// finally, find the largest cluster and save it into the array largest[]
		
		if(totalbonds>0)
		{
		  
           for (int i2=0; i2<totalbonds; i2++)     //the loop of finding different clusters
           {
        	   if((RBH[i2]!=-1)&(RBT[i2]!=-1))
        	   {
        		   for(int tpi=0; tpi<SN; tpi++)           // everytime reset temp[]
        			   temp[tpi]=-2;
        		   pin=0;
        		   temp[pin]=RBH[i2];
        		   RBH[i2]=-1;
        		   //PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/check.txt", temp[pin]);
        		   pin++;
        		   temp[pin]=RBT[i2];
        		   RBT[i2]=-1;
        		   //PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/check.txt", temp[pin]);
        		   pin++;
        		   for(int i3=0; i3< pin; i3++)       //the loop of scanning all the spins in temp[]
        		   {
        			   
        			   for(int i1=0; i1<totalbonds; i1++)   //the loop of scanning all the spins in RB[]
        			   {
        				   if(temp[i3]==RBH[i1])
        				   {
        					   writeintotemp=1;
        					   for(int scan1=0; scan1<pin; scan1++)
        						   {
        						   if(RBT[i1]==temp[scan1])
        							   writeintotemp=0;
        						   }
        					   if(writeintotemp==1)
        					   {temp[pin]=RBT[i1];
        					   //PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/check.txt", temp[pin]);
        					   RBT[i1]=-1;
        					   RBH[i1]=-1;
        					   pin++;}
        					   
        				   }
        				   if(temp[i3]==RBT[i1])
        				   {
        					   writeintotemp=1;
        					   for(int scan2=0; scan2<pin; scan2++)
        						   {
        						   if(RBH[i1]==temp[scan2])
        							   writeintotemp=0;
        						   }
        					   if(writeintotemp==1)
        					   {temp[pin]=RBH[i1];
        					   //PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/check.txt", temp[pin]);
        					   RBT[i1]=-1;
        					   RBH[i1]=-1;
        					   pin++;} 
        				   }   
        			   }
        		   }
        		   //PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/check.txt", "1111111111111111111");
        		   int copytemp[]= new int[ISP.M];   //we need to covert temp[SN]->copytemp[ISP.M] in order to match the cluster generating function
        		   for(int c=0;c<ISP.M;c++)
        		   {
        			   copytemp[c]=-2;
        		   }
        		   for(int ct=0;ct<SN;ct++)
        		   {
        			   copytemp[ct]=temp[ct];
        		   }
        		   cltemp=new Cluster(ISP.L1,ISP.L2,copytemp);
        		   CS.AddCluster(cltemp);   
        	   }
           } // the end of i2 cycle
           CS.ordering();
           CS.findmaximum();
           CS.findminimum();
		}

	}

	public void NNMapping(int Pseed)   // for nearest neighbor system with R=0
	{
		int stablespin[]=new int[ISP.M];
		int SN=0;
		// first, find all the stable spins
		for(int q=0; q<ISP.M; q++)
		{
			if(ISP.spin[q]==-1)
			{
				stablespin[SN]=q;
				SN++;
			}	
		}
		
		// second, determine all the possible bonds
		int PN=0;
		int PBsize=0;
		PBsize=(int)(SN*2);
		int RBH[]= new int [PBsize];
		int RBT[]= new int [PBsize];
		
		for(int pb=0; pb<PBsize; pb++)
		{
			RBH[pb]=-1;
			RBT[pb]=-1;
		}             //initialize all the possible bonds data structure
		
		Random bondsrand;
		int pin;
		bondsrand = new Random(Pseed);
		//boolean nearestneighbor;
		
		for (int k=0; k<SN; k++)
			for(int a=0; a<2; a++)  //only consider 2 neighbors in order to avoid double counting
			{
				int j=ISP.Nneighber(a, stablespin[k]);
				if(ISP.spin[j]==-1)
				{
					if(bondsrand.nextDouble()<=probability)		// third, throw random numbers on all the possible bonds to determine the percolation bonds
				   {
					RBH[PN]=stablespin[k];
					RBT[PN]=j;
					PN++;
				    }
				}
			}
		
		int totalbonds=PN;
		int temp[]= new int [SN];     //store the locations of the sites belonging to current cluster
		int writeintotemp;
		for(int tp=0; tp<SN; tp++)
		{
			temp[tp]=-2;
		}
		// finally, find the largest cluster and save it into the array largest[]
		
		if(totalbonds>0)
		{
		  
           for (int i2=0; i2<totalbonds; i2++)     //the loop of finding different clusters
           {
        	   if((RBH[i2]!=-1)&(RBT[i2]!=-1))
        	   {
        		   for(int tpi=0; tpi<SN; tpi++)           // everytime reset temp[]
        			   temp[tpi]=-2;
        		   pin=0;
        		   temp[pin]=RBH[i2];
        		   RBH[i2]=-1;
        		   //PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/check.txt", temp[pin]);
        		   pin++;
        		   temp[pin]=RBT[i2];
        		   RBT[i2]=-1;
        		   //PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/check.txt", temp[pin]);
        		   pin++;
        		   for(int i3=0; i3< pin; i3++)       //the loop of scanning all the spins in temp[]
        		   {
        			   
        			   for(int i1=0; i1<totalbonds; i1++)   //the loop of scanning all the spins in RB[]
        			   {
        				   if(temp[i3]==RBH[i1])
        				   {
        					   writeintotemp=1;
        					   for(int scan1=0; scan1<pin; scan1++)
        						   {
        						   if(RBT[i1]==temp[scan1])
        							   writeintotemp=0;
        						   }
        					   if(writeintotemp==1)
        					   {temp[pin]=RBT[i1];
        					   //PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/check.txt", temp[pin]);
        					   RBT[i1]=-1;
        					   RBH[i1]=-1;
        					   pin++;}
        					   
        				   }
        				   if(temp[i3]==RBT[i1])
        				   {
        					   writeintotemp=1;
        					   for(int scan2=0; scan2<pin; scan2++)
        						   {
        						   if(RBH[i1]==temp[scan2])
        							   writeintotemp=0;
        						   }
        					   if(writeintotemp==1)
        					   {temp[pin]=RBH[i1];
        					   //PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/check.txt", temp[pin]);
        					   RBT[i1]=-1;
        					   RBH[i1]=-1;
        					   pin++;} 
        				   }   
        			   }
        		   }
        		   //PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/check.txt", "1111111111111111111");
        		   int copytemp[]= new int[ISP.M];   //we need to covert temp[SN]->copytemp[ISP.M] in order to match the cluster generating function
        		   for(int c=0;c<ISP.M;c++)
        		   {
        			   copytemp[c]=-2;
        		   }
        		   for(int ct=0;ct<SN;ct++)
        		   {
        			   copytemp[ct]=temp[ct];
        		   }
        		   cltemp=new Cluster(ISP.L1,ISP.L2,copytemp);
        		   CS.AddCluster(cltemp);   
        	   }
           } // the end of i2 cycle
           CS.ordering();
           CS.findmaximum();
           CS.findminimum();
		}

	}
	
	public void fastNNMapping(int Pseed)   // for nearest neighbor system with R=0
	{
		int stablespin[]=new int[ISP.M];
		int copymap[]=new int[ISP.M];    // backup copy of the spin configuration, -1 is the stable spin, 0 is the dilution, -3 is the label if the site has already been chosen in the previous cluster
		int a0[]=new int[ISP.M];    //a0[j]=1 (there is bond between the up neighbor for site j and site j), otherwise a0[j]=0
		int a1[]=new int[ISP.M];    //a1[j]=1 (there is bond between the right neighbor for site j and site j), otherwise a0[j]=0
		int a2[]=new int[ISP.M];    //a2[j]=1 (there is bond between the down neighbor for site j and site j), otherwise a0[j]=0
		int a3[]=new int[ISP.M];    //a3[j]=1 (there is bond between the left neighbor for site j and site j), otherwise a0[j]=0
		
		
		int SN=0;
		// first, find all the stable spins
		for(int q=0; q<ISP.M; q++)
		{
			a0[q]=0;
			a1[q]=0;   //initialize the bond configuration to be null
			a2[q]=0;
			a3[q]=0;
			
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
			for(int a=0; a<4; a++)  //consider 4 neighbors with double counting
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
        		   cltemp=new Cluster(ISP.L1,ISP.L2,copytemp);
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
	
	public void ClusterSize()
	{
		double totalsize=0;
		totalclusters=0;
	    
		for(int i=0; clustersize[i]>=0; i++)
		{
			totalsize+=clustersize[i];
			totalclusters++;
		}
		if(totalclusters>0)
			meanclustersize=totalsize/totalclusters;
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