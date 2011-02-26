package kang.oldising;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import chris.util.PrintUtil;
import chris.util.Random;


import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DoubleValue;

public class percolationtest extends Simulation{
	
	Grid grid2 = new Grid("The largest cluster");
	Grid grid1 = new Grid("ising configuration snapshot");
	
	public int L1, L2, M; //parameters for the lattice
    public int i, x, y;  //parameters for the index
	public double percent;  //diluted percent
	public Random spinrand;
	public Random dilutionrand;
	
	
	
	
	
	public int R;   //interaction range R=0 is NN interaction
	public double magnetization;
	public int deadsites;
	public int isingspin[];     //the array of the data
	//public int initialcopy[];
	public int criticaldroplet[];
	
	// the new parameters for the percolation mapping20
	
	public double probability;//
	public int stablespin[];  //the array to store the location of all stable spins(-1), the range is from 0 to M
	public int SN;  //the number of the stable spins in the configuration
	public int PBH[];  //possible bonds head
	public int PBT[];  //possible bonds tail
	public int PBS[];  //possible bonds status
	public int PN;  //the number of the possible bonds in range R
	public int RBH[];
	public int RBT[];
	public int RN;
	public int totalbonds;  // the number of the real bonds after throwing the random number
	public int temp[]; //the array to save a cluster for finding the largest one
	public int largest[];
	
	
	
	public void percolation(int spin[], double percolationP){
		stablespin=new int[M];
		SN=0;
		// first, find all the stable spins
		for(int q=0; q<M; q++)
		{
			if(spin[q]==-1)
			{
				stablespin[SN]=q;
				SN++;
			}
				
		}
		
		// second, determine all the possible bonds
		int kx, ky, jx, jy;  // location for i,j
		int rx, ry;   // the difference between k and j
		PN=0;
		int PBsize=0;
		PBsize=(int)(SN*(2*R+1)*(2*R+1))/2;
		RBH= new int [PBsize];
		RBT= new int [PBsize];
		
		Random bondsrand;
		int pin;
		bondsrand = new Random(1);
		
		
		for (int k=0; k<SN; k++)
			for(int j=0; j<k; j++)
			{
				kx= stablespin[k]/L2;
				ky= stablespin[k]%L2;
				jx= stablespin[j]/L2;
				jy= stablespin[j]%L2;
				rx= kx-jx;
				ry= ky-jy;
				if( (rx*rx<=R*R) & (ry*ry<=R*R) )
					if(bondsrand.nextDouble()<=percolationP)		// third, throw random numbers on all the possible bonds to determine the percolation bonds
				{
					RBH[PN]=stablespin[k];
					RBT[PN]=stablespin[j];
					PN++;
				}
				
			}
		
		totalbonds=PN;
		
		temp= new int [SN];
		largest = new int [SN];
		int clustersizetemp;
		int writeintotemp;
		
		for(int tp=0; tp<SN; tp++)
		{
			temp[tp]=-2;
			largest[tp]=-2;
		}
		
		// finally, find the largest cluster and save it into the array largest[]
		
		if(totalbonds>0)
		{
           int largestsize=0;
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
        						   {if(RBT[i1]==temp[scan1])
        							   writeintotemp=0;}
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
        						   {if(RBH[i1]==temp[scan2])
        							   writeintotemp=0;}
        					   if(writeintotemp==1)
        					   {temp[pin]=RBH[i1];
        					   //PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/check.txt", temp[pin]);
        					   RBT[i1]=-1;
        					   RBH[i1]=-1;
        					   pin++;}
        					 
        				   }
        				   
        			   }
        			   
        			  
        		   }
        		   clustersizetemp= pin;
        		   //PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/check.txt", "1111111111111111111");
        		   if(clustersizetemp>largestsize)
        			   {
        			   largestsize=clustersizetemp;
        			   for(int tp=0; tp< SN; tp++)
        				   largest[tp]=temp[tp];
        			   }
        		   
        	   }
        	   
        	   
           }
		}

	}
	
	public void cluster(int largestcluster[])
	{
		for(int in=0; in<M; in++)
			{
			if(isingspin[in]==0)
		       criticaldroplet[in]=0;
			else
				criticaldroplet[in]=-1;
			}
		for(int cl=0; cl<SN; cl++)
			if(largestcluster[cl]!=-2)
				criticaldroplet[largestcluster[cl]]=-2;
	}
	
	public void initialize(int spin[], int type, double p, int spinseed, int dilutionseed)
	{
		spinrand= new Random(spinseed);
		dilutionrand= new Random(dilutionseed);
		deadsites=0;
		if(type==0)
			for (i=0; i<M; i++)
		    {
			   spin[i]=-1;
			   if (spinrand.nextDouble()> 0.5)
				   spin[i]=1;
				
		    }
		if(type==1)
			for (i=0; i<M; i++)
		    {
			   spin[i]=1;
		    }
		if(type==2)
			for (i=0; i<M; i++)
		    {
			   spin[i]=-1;
		    }
		
		for (i=0; i<M; i++){
			if(dilutionrand.nextDouble()<p)
			{
				spin[i]=0;
				deadsites++;
			}
			//initialcopy[i]=spin[i];  // here, make the copy of the system
		}    // here, the sequence of dilution and initialization is not important
		
	}

	public void animate(){
		
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);
		ising.setColor(-1, Color.WHITE);
		ising.setColor(0, Color.RED);
		ising.setColor(-2, Color.GREEN);
		
		grid1.setColors(ising);
		grid1.registerData(L1, L2, isingspin);
		grid2.setColors(ising);
		grid2.registerData(L1, L2, criticaldroplet);

		
	}
	
	public void clear(){
		grid1.clear();
		grid2.clear();
	}
	
	
	public static void main (String[] kangliu){
		new Control(new percolationtest(), "Kang Liu's percolationtest" );
	}
	
	
	public void load(Control percolationtest){
		percolationtest.frame (grid1);
		percolationtest.frame (grid2);
	
		params.add("lattice's width", 20);
		params.add("lattice's length", 20);
		params.add("Interaction range", 2);
		params.add("Diluted Percentage", new DoubleValue(0.30,0,1).withSlider());
		params.add("Percolation P", new DoubleValue(0.30,0,1).withSlider());
	}
	
	
	public void run(){
		
		
		R = (int)params.fget("Interaction range");
		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");
		M = L1 * L2;
		
		percent = params.fget("Diluted Percentage");
		probability = params.fget("Percolation P");
		
		criticaldroplet = new int[M];
		isingspin = new int[M];
		
		initialize(isingspin, 2, percent, 1, 1);
		
		Job.animate();
		
		percolation(isingspin, probability);
		
		cluster(largest);
		
		Job.animate();
		
		
	}
	
	
	
	
	
	
}