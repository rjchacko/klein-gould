package kang.ising;

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

public class biasintervention extends Simulation{
	
	private DecimalFormat fmt = new DecimalFormat("0000000");
	
	Grid grid1 = new Grid("Bias Diluted Ising lattice 2d for nucleation");
	Grid grid2 = new Grid("Bias Diluted Ising initial copy");
	Grid grid3 = new Grid("Bias Diluted Ising bias critical droplet");
	Grid grid4 = new Grid("Bias Diluted Ising bias intervention");
	Grid grid5 = new Grid("Percolation cluster");

	

	public int L1, L2, M; //parameters for the lattice
    public int i, x, y;  //parameters for the index
	public double T;     //temperature
	public double QuenchT;  //the temperature after the quench
	public double H;     //field
	public double field;
	public double J;     //interaction constant after normalization
	public double NJ;    //interaction constant before normalization
	public double percent;  //diluted percent


	public double magnetization;
	public double Ms;  //saturated magnetization
	public int R;   //interaction range R=0 is NN interaction
	
	//bias variables
	public int biasrange;
	public double biaspercent;
	

	//public int numberofcopies;
	public int lifetime;   // the time starts at the moment when quench of the fields happens
	public int quenchstart; // the time of the quench

	public int deadsites;
	
	public int prestep;  // step for preparation of intial configuration
	public int step; //Monte Carlo step
	public int rationalMCS; // use this we can run half of the MCS
	public int steplimit; //upperlimit of MC step
	public int meanlifetime;
	public int totallifetime;
	
	public int growthpercentage;
	public int decaynumber;
	public int grownumber;
	public int interventionsteplimit;
	public int interventioncopies;
	public double interventionstart;    // this double format is used for pinpoint the critical droplet
	public double interventionP;
	
	

	public int isingspin[];     //the array of the data
	public int initialcopy[];   //the array of the initial copy of the system
	public int isingcopy[];     //the array for reset
	public int interventioncopy[];

	public Random spinfliprand;
    public Random spinrand;
	public Random dilutionrand;
	public Random newpositionrand;
	public Random newspinfliprand;
	public Random biasrand;
	

	public int spinflipseed;
	public int spinseed;
	public int dilutionseed;
	public int biasseed;
	public int type;  // initialization type(0-random spin    1-all spins up    2-all spins down)
	
	
	public double probability;//
	public double percolationM; //the magnetization for the percolation mapping
	public int largestsize;
	
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
	public int percolationcluster[];

	public int min, max;
	public double avgmag[];
	
 	public int Nneighber(int a,int i )
 	{// function for the index of nearest neighbor
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

 	public double interactionEchange (int range, double constJ, int spin[], int j)
 	{ //function for interaction energy
		double Energy=0;
		double Energychange=0;
		if (range==0)
		{
			J= constJ/4;                     // normalize the interaction constant
			int b,k;
		    for(b=0; b<4;b++)
		    {
			k=Nneighber(b,j);
			Energy=Energy+J*spin[j]*spin[k];
		    }
		    Energychange=-2*Energy;
		    
		}
		if (range!=0)
		{
			int S=0;
			int nx=j/L2;
			int ny=j%L2;
			int kx, ky;
			J= constJ/((2*range+1)*(2*range+1)-1);
			for (int m=-range; m<=range; m++)
				for (int n=-range; n<=range; n++)
				{
					kx=nx+m;
					ky=ny+n;
					if(nx+m<0)
						kx=nx+m+L1;
					if(nx+m>L1-1)
						kx=nx+m-L1;
					if(ny+n<0)
					    ky=ny+n+L2;
					if(ny+n>L2-1)
						ky=ny+n-L2;
					S+=spin[kx*L2+ky];	
				}
			Energy=J*spin[j]*S-J;
			Energychange=-2*Energy;
			
		}
		
		return Energychange;	
    }
	
 	public void spinflip(int range, double constJ, int spin[], int j, Random flip, double temperature, double field)
	{
		double ZeemanE=2*field*spin[j]; //the change in Zeeman's energy if we flip the spin[j]
		double Echange=ZeemanE+interactionEchange(range, constJ, spin,j);
		
		if(Echange<0)
		{
			spin[j]=-spin[j];
		}
		
		else
		{
			if(flip.nextDouble()<=Math.exp(-Echange/temperature))
					spin[j]=-spin[j];
		}
		
	}

 	public double magnetization(int spin[])   // function to calculate the magnetization for spin[]
	{
		double totalmagnetization=0;	
		for(int s=0; s<M; s++)
		{
			totalmagnetization+=spin[s];
		}
		double m=totalmagnetization/M;
		return m;
	}
	
	public void MCS(int spin[], Random flip, double ratio)
	{
	    int j=0;
		H = params.fget("Field");
	    T = params.fget("Temperature");
	    NJ = params.fget("Interaction Constant");
	    rationalMCS= (int) (ratio*M);
	    for (int f=0; f< rationalMCS; f++)
	    {
		   j=(int) (flip.nextDouble()*M);
		   spinflip(R, NJ, spin, j, flip, T, H);
		  
	    }
	    magnetization=magnetization(spin);
	    params.set("magnetization", magnetization);
	    
	}

	public void temperaturequench(double finaltemp)
	{
		params.set("Temperature", finaltemp);
	}
	
	public void fieldquench()
	{
	    double finalfield=-H;
	    params.set("Field", finalfield);
	    
	}
	
	public void animate()
	{
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);
		ising.setColor(-1, Color.WHITE);
		ising.setColor(0, Color.RED);
		ising.setColor(2, Color.BLUE);
		ising.setColor(-2, Color.GREEN);
		
		grid1.setColors(ising);
		grid1.registerData(L1, L2, isingspin);
		grid2.setColors(ising);
		grid2.registerData(L1, L2, initialcopy);
		grid3.setColors(ising);
		grid3.registerData(L1, L2, interventioncopy);
		grid4.setColors(ising);
		grid4.registerData(L1, L2, isingcopy);
		grid5.setColors(ising);
		grid5.registerData(L1, L2, percolationcluster);

		
	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
		grid3.clear();
		grid4.clear();
		grid5.clear();

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
			initialcopy[i]=spin[i];  // here, make the copy of the system
			//isingcopy[i]=spin[i];
		}    // here, the sequence of dilution and initialization is not important
		
		
	}
	
	public void biasdilution (int spin[], int biasR, double biasP, int biasSeed)  // function for dilution bias
	{
		int cx, cy; //x,y index for the center
		cx=(int)L1/2;
		cy=(int)L2/2;
		int bx, by; //x,y for the bias point
		biasrand= new Random(biasSeed);
		
		
		for(bx=cx-biasR; bx<cx+biasR; bx++)
			for(by=cy-biasR; by<cy+biasR; by++)
			{
				if(spin[bx*L2+by]!=0)
					if(biasrand.nextDouble()<biasP)
						{
						spin[bx*L2+by]=0;
						isingcopy[bx*L2+by]=0;
						deadsites++;
						initialcopy[bx*L2+by]=2;
						}
					
			}
		
		
	}
	
	public void movie(Grid grid, int number, int copynumber)   //function to capture the grid
	{
		
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/criticaldroplet/homogeneous/q=0.111/evolution1880-1890/pic_"+fmt.format(copynumber)+"_"+fmt.format(number)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
	}
	
	public void movie2(Grid grid, int size, double etime)   //function to capture the grid
	{
		
		PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/criticaldroplet/homogeneous/q=0.111/evolution1880-1890/evolution.txt", etime, size);
			
		    String SaveAs = "/Users/liukang2002507/Desktop/simulation/criticaldroplet/homogeneous/q=0.111/evolution1880-1890/pic_"+fmt.format(etime)+"_"+fmt.format(size)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
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
	
	public void printdilutionratio(int spin[], int estep, int dropletcenter, int clustercenter)
	{
	    int distance;
	    distance=distance(dropletcenter, clustercenter);
	    int drx,dry;
	    int clx,cly;
	    drx=dropletcenter/L2;
	    dry=dropletcenter%L2;
	    clx=clustercenter/L2;
	    cly=clustercenter%L2;
	    double dropletratio;
	    double clusterratio;
	    dropletratio= dilutionratio(spin, R, dropletcenter, M);
	    clusterratio= dilutionratio(spin, R, clustercenter, M);
	    
		PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/criticaldroplet/homogeneous/q=0.111/evolution1880-1890/dilutionratio.txt", estep, drx, dry, clx, cly, distance, dropletratio, clusterratio );
	}
	
    
	
	public void intervention(int copies, int steplimit, double percentage)
	{
		//growthpercentage=0;
		decaynumber=0;
		grownumber=0;
		for(int k=0; k<copies; k++)
		{
			//newpositionrand=new Random(k);
			newspinfliprand=new Random(k+spinflipseed);

			for(int l=0; l<M; l++)
				isingcopy[l]=interventioncopy[l];
			for(int s=0; s<steplimit; s++)
			{
				MCS(isingcopy, newspinfliprand, 1);
				Job.animate();
				
			}
			movie(grid4, steplimit, k);
			if(magnetization>(percentage*Ms))
				decaynumber++;
			else
				grownumber++;
			
			params.set("u#copies", k);
			params.set("grow",grownumber);
			params.set("decay",decaynumber);
			
			
		}
		
		
	}
	
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
		   largestsize=0;
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
	
	public void cluster(int largestcluster[], int spin[])
	{
		for(int in=0; in<M; in++)
			{
		       percolationcluster[in]=initialcopy[in];
			if(initialcopy[in]!=2)
				percolationcluster[in]=spin[in];
			}
		for(int cl=0; cl<SN; cl++)
			if(largestcluster[cl]!=-2)
				percolationcluster[largestcluster[cl]]=-2;
	}
	
	public void clusterevolution(int spin[], double steplength, int estart)
	{
		int dropletC, clusterC;   // centers of the droplet and cluster
		for(int estep=0; estep<(int)(10/steplength); estep++)
		{
		percolationM=magnetization;
		probability=1- Math.exp(-(-J*(1+percolationM)/T));
		params.set("percolationP", probability);
		
		percolation(spin, probability);
		cluster(largest, spin);
		
		Job.animate();
		
		dropletcenter(spin, R, M);
		dropletC= min;
		clusterC= clustercenter(largest);
		
		movie2(grid5, largestsize, estep);  // the percolation cluster
		printdilutionratio(spin, estep, dropletC, clusterC);
		params.set("percolationP", probability);
		params.set("clustersize", largestsize);
		
		MCS(spin, spinfliprand, steplength);
		}
		
		
		
	}
	
	public double dilutionratio(int spin[], int r, int i, int total)
	{
		double ratio;
		double dilutedsite;
		dilutedsite=0;
		double totalinrange;
		totalinrange=0;
		int j;
		int disij;

		for(j=0; j<total;j++)
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
	
	public void dropletcenter(int spin[], int r, int total)
	{
		int i,ix,iy;
		double sum;
		int mx,my,j;
		int jx,jy;
		
		for(i=0;i<total;i++)
			{
			ix=i/L2;
			iy=i%L2;
			sum=0;
			for(mx=-r;mx<=r;mx++)
				for(my=-r;my<=r;my++)
				{
					jx=ix+mx;
					jy=iy+my;
					if(ix+mx<0)
						jx=ix+mx+L1;
					if(ix+mx>L1-1)
						jx=ix+mx-L1;
					if(iy+my<0)
					    jy=iy+my+L2;
					if(iy+my>L2-1)
						jy=iy+my-L2;
					j=jx*L2+jy;
					sum+=spin[j];
				}
			avgmag[i]=sum/((2*r+1)*(2*r+1));
			}
		max=0;
		min=0;
		for(i=1;i<total;i++)
		{
			if(avgmag[i]>avgmag[max])
				max=i;
			if(avgmag[i]<avgmag[min])
				min=i;
		}
		
	}
	
	public int clustercenter(int largestcluster[])
	{
		int fp,fx,fy;
		fp=largestcluster[0];
		fx=fp/L2;
		fy=fp%L2;
		double sumx,sumy;
		sumx=0;
		sumy=0;
		int center,cx,cy;
		int clx[];
		int cly[];
		int tempx,tempy;
		clx= new int[SN];
		cly= new int[SN];
		clx[0]=fx;
		cly[0]=fy;
		sumx=fx;
		sumy=fy;
		
		for(int cl=1; cl<SN; cl++)
		{
			tempx=largestcluster[cl]/L2;
			tempy=largestcluster[cl]%L2;
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
		cx=(int)(sumx/SN);
		cy=(int)(sumy/SN);
		
		if(cx<0)
			cx=cx+L1;
		if(cx>L1-1)
			cx=cx-L1;
		if(cy<0)
			cy=cy+L2;
		if(cy>L2-1)
			cy=cy-L2;
		
		
		center=cx*L2+cy;
		return center;
	}
	
	
	
	
	public static void main (String[] biasintervention){
		new Control(new biasintervention(), "Kang Liu's bias diulted ising model" );
	}
	
	public void load(Control biasintervention){
		biasintervention.frame (grid1);
		biasintervention.frame (grid2);
		biasintervention.frame (grid3);
		biasintervention.frame (grid4);
		biasintervention.frame (grid5);
		
		
		params.add("lattice's width", 200);
		params.add("lattice's length", 200);
		params.add("Diluted Percentage", new DoubleValue(0.111,0,1).withSlider());
		params.add("Bias percent", new DoubleValue(0, 0, 1).withSlider());
		
		params.add("Quench starts at", 100);
		
		params.addm("Interaction Constant", -4);
		params.add("Interaction range", 10);
		params.add("Bias range", 10);
		params.add("Monte Carlo step's limit", 1000000);

		//params.add("Spin seed", 1);
		//params.add("Dilution seed", 1);
		//params.add("Bias seed", 1);
		//params.add("Initialization type", 0);
		params.addm("Quench temperature", new DoubleValue(1.591, 0, 10).withSlider());
		params.addm("Temperature", new DoubleValue(9, 0, 10).withSlider());
		params.addm("Field", new DoubleValue(1.022, -5, 5).withSlider());
		
		
		params.add("intervention start", new DoubleValue(2880, 0, 99999999).withSlider());
		params.add("intervention steplimit", 50);
		params.add("intervention copies", 50);
		params.add("intervention percentage", new DoubleValue(0.85,0,1).withSlider());
		
		///below is the parameters that would be displayed on the panel
		
		params.add("MC time");
		params.add("Metastable state time");
		params.add("magnetization");
		params.add("Saturated magnetization");
		params.add("lifetime");

		
		
        //params.add("spinfliprand");
        //params.add("spinrand");
        //params.add("dilutionrand");
        params.add("u#copies");
        params.add("grow");
        params.add("decay");
        params.add("percolationP");
        params.add("clustersize");
 

		
	}
	

    public void run(){
		
    	int estart;   //start of the cluster evolution
    	
		
		R = (int)params.fget("Interaction range");
		biasrange= (int)params.fget("Bias range");
	    steplimit = (int)params.fget("Monte Carlo step's limit");


		quenchstart=(int)params.fget("Quench starts at");
		QuenchT=params.fget("Quench temperature");

		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");

		M = L1 * L2 ;

		
		isingspin = new int[M];
		initialcopy = new int[M];
		isingcopy = new int[M];
		interventioncopy = new int[M];
		percolationcluster = new int [M];
		avgmag = new double[M];
		

		//spinseed = (int)params.fget("Spin seed");
		//dilutionseed = (int)params.fget("Dilution seed");
		//biasseed =(int)params.fget("Bias seed");
		spinseed=1;
		dilutionseed=1;
		biasseed=1;
	
		percent = params.fget("Diluted Percentage");
		biaspercent = params.fget("Bias percent");
		//type = (int)params.fget("Initialization type");
		type=0;
		field = params.fget("Field");
		
		interventionstart=params.fget("intervention start");
		interventioncopies=(int)params.fget("intervention copies");
		interventionsteplimit=(int)params.fget("intervention steplimit");
		interventionP=params.fget("intervention percentage");
		
		initialize(isingspin, type, percent, spinseed, dilutionseed);
		biasdilution(isingspin, biasrange, biaspercent, biasseed);
		Job.animate();
		int biasint= (int)(biaspercent*100);
		movie(grid2, biasint, 6666);// save the initial configuration with the bias dilution
		
		spinflipseed = (int)((field+1)*QuenchT*10000);
		spinfliprand= new Random (spinflipseed);
		
		for (prestep=0; prestep < 50; prestep++)
		{
		MCS(isingspin, spinfliprand ,1);
		params.set("MC time", prestep-50);
		Job.animate();
		}
		
		//params.set("spinfliprand", spinflipseed);
		//params.set("dilutionrand",dilutionrand.nextDouble());
		temperaturequench(QuenchT);
		Ms=magnetization;
		
		for(step=0; steplimit>0; step++)
		{
			if(step==quenchstart)
			{
			fieldquench();	
			}
			if(step==quenchstart+1)
			{
				Ms=magnetization;
				params.set("Saturated magnetization", Ms);
			}
			

			
			if(step-quenchstart==(int)interventionstart)
			{
				double residue=0;
				residue=interventionstart-step+quenchstart;
				if(residue!=0)
					MCS(isingspin, spinfliprand, residue);
				
				percolationM=magnetization;
				probability=1- Math.exp(-(-J*(1+percolationM)/T));
				params.set("percolationP", probability);
				for (int h=0; h<M; h++)
				{
				interventioncopy[h]=isingspin[h];
				}
				
				
				intervention(interventioncopies, interventionsteplimit, interventionP);
				steplimit=-1;
				
				
				
				//now do the percolation mapping
				
				percolation(interventioncopy, probability);
				cluster(largest, interventioncopy);
				
				Job.animate();
				movie(grid3, 8888, 8888);  // the critical droplet configuration
				movie(grid5, 1234, step-quenchstart);  // the percolation cluster
				params.set("percolationP", probability);
				params.set("clustersize", largestsize);
			}
			
			
			estart=1880;
			if(step-quenchstart==estart)
			{
				clusterevolution(isingspin, 0.1, estart);
				steplimit=-2;
			}
			
			
			else{
			MCS(isingspin, spinfliprand, 1);
			params.set("MC time", step);
			Job.animate();
			params.set("Metastable state time", step-quenchstart);
			if(interventionstart>9999)
				movie(grid1, step, 0000);
			}
			//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/biasdata/nucleation/data50.txt",step-quenchstart, (magnetization/Ms));
		}
		
	
	
    }// the end of run()
	
}