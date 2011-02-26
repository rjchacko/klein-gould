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

public class nucleation extends Simulation{
	
	private DecimalFormat fmt = new DecimalFormat("0000000");
	
	Grid grid1 = new Grid("Diluted Ising lattice 2d for nucleation");
	Grid grid2 = new Grid("Diluted Ising copy");
	Grid grid3 = new Grid("Critical Droplet snapshot");
	Grid grid4 = new Grid("The percolation cluster");
	

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

	public int numberofcopies;
	public int lifetime;   // the time starts at the moment when quench of the fields happens
	public int quenchstart; // the time of the quench
	public int interventionstart;
	public int deadsites;
	public int prestep;  // step for preparation of intial configuration
	public int step; //Monte Carlo step
	public int steplimit; //upperlimit of MC step
	public int interventionsteplimit;
	public double interventionthreshold;

	public int isingspin[];     //the array of the data
	public int initialcopy[];   //the array of the initial copy of the system
	public int interventioncopy[];   //the array of the isingspin at the beginning of intervention
	public int isingcopy[];     //the array of the copy
	public int criticaldroplet[];

	public double Mag[];   //the array for magnetization used for the susceptibiliy
	
	//public Random positionrand;
	public Random spinfliprand;
    public Random spinrand;
	public Random dilutionrand;
	public Random newpositionrand;
	public Random newspinfliprand;
	
	//public int positionseed;
	public int spinflipseed;
	public int spinseed;
	public int dilutionseed;
	public int type;  // initialization type(0-random spin    1-all spins up    2-all spins down)
	
	public int decaynumber;
	public int grownumber;
	public double growthpercentage;
	
	public int Xstart;
	public int Xend;
	public int Mend;  // the end of recording magnetization
	public double X;
	
	
	// the new parameters for the percolation mapping
	
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
	public int percolationcluster[];
	
	
	
	
	
	
	
	
	
////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
 	public int Nneighber(int a,int i ){// function for the index of nearest neighbor
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
	
	public double interactionEchange (int range, double constJ, int spin[], int j){ //function for interaction energy
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

	public void movie(Grid grid, int number, int copynumber)   //function to capture the grid
	{
		
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/spinodaldata/p=0.2/pic_"+fmt.format(copynumber)+"_"+fmt.format(number)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
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
		}    // here, the sequence of dilution and initialization is not important
		
	}
	
	public void MCS(int spin[], Random flip)
	{
	    int j=0;
		H = params.fget("Field");
	    T = params.fget("Temperature");
	    NJ = params.fget("Interaction Constant before normalization");
	    for (int f=0; f< M; f++)
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
	
	public void intervention(int copies, int steplimit, double percentage)
	{
		growthpercentage=0;
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
				MCS(isingcopy, newspinfliprand);
				Job.animate();
				
			}
			movie(grid2, steplimit, k);
			if(magnetization>(percentage*Ms))
				decaynumber++;
			else
				grownumber++;
			
			params.set("#copies", k);
			//params.set("Growthpercent",growthpercentage);
			params.set("grow",grownumber);
			params.set("decay",decaynumber);
			
			
		}
		
		
	}
	
	public void susceptibility(double M[]){
		double m=0;
		double m2=0;
		for(int r=Xstart; r<Xend; r++)
		{
			m+=M[r];
			m2+=(M[r]*M[r]);
		}
		X=(m2/(Xend-Xstart))-(m*m)/((Xend-Xstart)*(Xend-Xstart));
		params.set("susceptibility", X);
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
		       percolationcluster[in]=0;
			else
				percolationcluster[in]=criticaldroplet[in];
			}
		for(int cl=0; cl<SN; cl++)
			if(largestcluster[cl]!=-2)
				percolationcluster[largestcluster[cl]]=-2;
	}
	
	
	public void animate(){
		
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);
		ising.setColor(-1, Color.WHITE);
		ising.setColor(0, Color.RED);
		ising.setColor(-2, Color.GREEN);    // the color to represent the largest percolation cluster
		
		grid1.setColors(ising);
		grid1.registerData(L1, L2, isingspin);
		grid2.setColors(ising);
		grid2.registerData(L1, L2, isingcopy);
		grid3.setColors(ising);
		grid3.registerData(L1, L2, criticaldroplet);
		grid4.setColors(ising);
		grid4.registerData(L1, L2, percolationcluster);
		
	}
	
	public void clear(){
		grid1.clear();
		grid2.clear();
		grid3.clear();
		grid4.clear();
	}

	
	public static void main (String[] kangliu){
		new Control(new nucleation(), "Kang Liu's diulted ising model for nucleation" );
	}
	
	
	public void load(Control nucleation){
		nucleation.frame (grid1);
		nucleation.frame (grid2);
		nucleation.frame (grid3);
		nucleation.frame (grid4);
		
		params.add("lattice's width", 400);
		params.add("lattice's length", 400);
		params.add("Diluted Percentage", new DoubleValue(0.30,0,1).withSlider());
		params.add("Intervention copies", 20);
		params.add("Quench starts at", 100);
		params.add("Intervention starts at",1034000000);
		params.add("Mag ends", 500);
		params.add("Xstart",100);
		params.add("Xend",400);
		params.add("Intervention threshold",new DoubleValue(0.7,-1,1).withSlider());
		params.add("Intervention steplimit", 20);
		
		params.addm("Interaction Constant before normalization", -4);
		params.add("Interaction range", 20);
		
		params.add("Monte Carlo step's limit", 1000000);
		//params.add("Position seed", 1001);
		params.add("Spinflip seed", 1001);
		params.add("Spin seed", 1);
		params.add("Dilution seed", 1);
		params.add("Initialization type", 0);
		
		params.addm("Quench temperature", new DoubleValue(2.835, 0, 10).withSlider());
		params.addm("Temperature", new DoubleValue(9, 0, 10).withSlider());
		params.addm("Field", new DoubleValue(0, -5, 5).withSlider());
		params.add("MC time");
		params.add("Metastable state time");
		params.add("magnetization");
		params.add("Saturated magnetization");
		params.add("interventionM");
		params.add("#copies");
 
        params.add("grow");
        params.add("decay");
       // params.add("positionrand");  
        params.add("spinfliprand");
        params.add("spinrand");
        params.add("dilutionrand");
        params.add("susceptibility");

		
	}
	
	public void run(){
		
	
		R = (int)params.fget("Interaction range");

		steplimit = (int)params.fget("Monte Carlo step's limit");
		numberofcopies=(int)params.fget("Intervention copies");
		interventionstart=(int)params.fget("Intervention starts at");
		quenchstart=(int)params.fget("Quench starts at");
		QuenchT=params.fget("Quench temperature");
		interventionsteplimit=(int)params.fget("Intervention steplimit");
	    interventionthreshold=params.fget("Intervention threshold");

		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");
		M = L1 * L2;
		Mend=(int)params.fget("Mag ends");
		
		isingspin = new int[M];
		initialcopy = new int[M];
		isingcopy = new int[M];
		interventioncopy =new int[M];
		criticaldroplet = new int[M];
		percolationcluster = new int[M];
		
		Mag = new double [Mend];
		Xstart= (int)params.fget("Xstart");
		Xend= (int)params.fget("Xend");
	//	positionseed = (int)params.fget("Position seed");
		//spinflipseed = (int)params.fget("Spinflip seed");
		spinseed = (int)params.fget("Spin seed");
		dilutionseed = (int)params.fget("Dilution seed");
		
		percent = params.fget("Diluted Percentage");
		type = (int)params.fget("Initialization type");
		
	//	positionrand= new Random (positionseed);
		field = params.fget("Field");

		


		
		initialize(isingspin, type, percent, spinseed, dilutionseed);
		
		Job.animate();
		//movie(grid1, 99999, 1);
		//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/spinodaldata/errorbar200HX.txt", "         ");
		for(int u=0; u<16; u++)
		{
			params.set("Temperature",9);	//reset the high temperature everytime
			params.set("Field",field);
			params.set("Quench temperature",2.835-0.005*u);
			QuenchT=params.fget("Quench temperature");
			spinflipseed = (int)((field+1)*QuenchT*10000+u);
			spinfliprand= new Random (spinflipseed);
	
		for (prestep=0; prestep < 50; prestep++)
			{
			MCS(isingspin, spinfliprand);
			params.set("MC time", prestep-50);
			Job.animate();
			}
		
		params.set("spinfliprand", spinflipseed);
		params.set("spinrand", u);
		params.set("dilutionrand",dilutionrand.nextDouble());
		temperaturequench(QuenchT);
		
		//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/spinodaldata/errorbar200.txt","field = ", H , u);
		//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/spinodaldata/errorbar200.txt","             ");
		
		for (step=0; step-quenchstart<Mend+1; step++)
		{

			if(step==quenchstart)
				{
				fieldquench();
				}
			if(step==quenchstart+100)
			{
				Ms=magnetization;
				params.set("Saturated magnetization", Ms);
			}
			MCS(isingspin, spinfliprand);
			params.set("MC time", step);
			Job.animate();
			//movie(grid1, step, 1001);
			params.set("Metastable state time", step-quenchstart);
			if(((step-quenchstart)%1==0)&((step-quenchstart)>=0)&(interventionstart>100000)&(step-quenchstart<Mend))
				
				{
				//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/nucleationmovie/data200.txt",step-quenchstart, (magnetization/Ms));
				Mag[step-quenchstart]=magnetization;
				//movie(grid1,step-quenchstart, 1001);
				}
			if((step-quenchstart)==Mend)
				{
				/*if(magnetization<0)
					{
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/spinodaldata/errorbar200.txt","lifetime too short", H , u );
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/spinodaldata/errorbar200HX.txt","lifetime too short", H , u );
					}*/
				susceptibility(Mag);
				PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/spinodaldata/400large300Tc.txt",H , QuenchT, u, Mend, Xend-Xstart, X);
				//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/spinodaldata/errorbar200HX.txt",H , QuenchT, u, Mend, Xend-Xstart, X);
				}

		    	
		}
		
	}
		for (int h=0; h<M; h++)
			{
			interventioncopy[h]=isingspin[h];
			criticaldroplet[h]=interventioncopy[h];
			}
		params.set("interventionM", magnetization(interventioncopy));
		movie(grid1, interventionstart, 99999);

		
		//the beginning of intervention
		intervention(numberofcopies, interventionsteplimit, interventionthreshold);
	
		Job.animate();
		movie(grid3, interventionstart, 88888);
		
		//here, start the percolation and save the data onto grid4

        
		//end of all MCS steps
		
		
	}


	
	
	
		
}