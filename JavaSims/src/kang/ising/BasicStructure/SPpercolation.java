package kang.ising.BasicStructure;

import chris.util.Random;
import kang.ising.BasicStructure.Cluster;
import kang.ising.BasicStructure.ClusterSet;

/*the way to use this class is 
 * 
 * 1  construct Percolation object percolation
 * 2  input ising strucuture of the critical droplet configuration(default setting +1 is the stable direction)
 * 3  input the dilution map for the calculation of local bond probability--use SetDensity()
 * 3  use  SetProbability(T);
 * 4  use  Mapping();
 * 5  the set of clusters are stored in SPpercolation.CS;
 * 6  the order in size of the clusters can be found in SPpercolation.CS.order[];
 */


public class SPpercolation{
	
	public IsingStructure ISP;
	public ClusterSet CS;
	public Cluster cltemp;
	public double rho[];      //density of the occupied sites
	public double probability[];    //the distribution of bond probability
	public double Pb;
	public int stabledirection;
	
	public SPpercolation()
	{
		
	}

	public SPpercolation(IsingStructure IS, int number)
	{
		this.ISP= IS.clone();
		this.CS= new ClusterSet(number);
		rho= new double[IS.M];
		probability=new double[IS.M];
		stabledirection=1;
	}

	
	public void SetDensity()
	{
		ISP.dilutionmap(ISP.R, "Square");
	
		for(int j=0; j<ISP.M; j++)
		{
			rho[j]=ISP.dilutionmap[j];
		}
	}
	
	
	public void SetProbability(double T)
	{
		double percolationM=ISP.magnetization;
		for(int ii=0; ii<ISP.M; ii++)
		{
			if(stabledirection==1)
			probability[ii]=1- Math.exp(-(-2*ISP.J*(rho[ii]-percolationM))/T);
			}
		
	}
	
	public void Mapping()    //for long range system with R as the interaction
	{
		int stablespin[]=new int[ISP.M];
		int SN=0;
		// first, find all the stable spins
		for(int q=0; q<ISP.M; q++)
		{
			if(ISP.spin[q]==stabledirection)
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
				
				Pb=(probability[stablespin[k]]+probability[stablespin[j]])/2;
				
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
					if(bondsrand.nextDouble()<=Pb)		// third, throw random numbers on all the possible bonds to determine the percolation bonds
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
	
}
