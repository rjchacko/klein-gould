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
	public int steplimit; //upperlimit of MC step
	public int meanlifetime;
	public int totallifetime;
	
	public int growthpercentage;
	public int decaynumber;
	public int grownumber;
	public int interventionsteplimit;
	public int interventioncopies;
	public int interventionstart;
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
	
	public void animate()
	{
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);
		ising.setColor(-1, Color.WHITE);
		ising.setColor(0, Color.RED);
		ising.setColor(2, Color.BLUE);
		
		grid1.setColors(ising);
		grid1.registerData(L1, L2, isingspin);
		grid2.setColors(ising);
		grid2.registerData(L1, L2, initialcopy);
		grid3.setColors(ising);
		grid3.registerData(L1, L2, interventioncopy);
		grid4.setColors(ising);
		grid4.registerData(L1, L2, isingcopy);

		
	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
		grid3.clear();
		grid4.clear();

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
		
			String SaveAs = "/Users/liukang2002507/Desktop/biasdata/nucleation/p0bias50/pic_"+fmt.format(copynumber)+"_"+fmt.format(number)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
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
				MCS(isingcopy, newspinfliprand);
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
	
	public static void main (String[] biasintervention){
		new Control(new biasintervention(), "Kang Liu's bias diulted ising model" );
	}
	
	public void load(Control biasintervention){
		biasintervention.frame (grid1);
		biasintervention.frame (grid2);
		biasintervention.frame (grid3);
		biasintervention.frame (grid4);
		
		
		params.add("lattice's width", 200);
		params.add("lattice's length", 200);
		params.add("Diluted Percentage", new DoubleValue(0,0,1).withSlider());
		params.add("Bias percent", new DoubleValue(0.5, 0, 1).withSlider());
		
		params.add("Quench starts at", 100);
		
		params.addm("Interaction Constant before normalization", -4);
		params.add("Interaction range", 10);
		params.add("Bias range", 10);
		params.add("Monte Carlo step's limit", 1000000);

		params.add("Spin seed", 1);
		params.add("Dilution seed", 1);
		params.add("Bias seed", 1);
		params.add("Initialization type", 0);
		
		params.addm("Quench temperature", new DoubleValue(1.778, 0, 10).withSlider());
		params.addm("Temperature", new DoubleValue(9, 0, 10).withSlider());
		params.addm("Field", new DoubleValue(0.95, -5, 5).withSlider());
		
		
		params.add("intervention start", 1912);
		params.add("intervention steplimit", 50);
		params.add("intervention copies", 20);
		params.add("intervention percentage", new DoubleValue(0.85,0,1).withSlider());
		
		///below is the parameters that would be displayed on the panel
		
		params.add("MC time");
		params.add("Metastable state time");
		params.add("magnetization");
		params.add("Saturated magnetization");
		params.add("lifetime");

		
		
        params.add("spinfliprand");
        params.add("spinrand");
        params.add("dilutionrand");
        params.add("u#copies");
        params.add("grow");
        params.add("decay");
 

		
	}
	

    public void run(){
		
		
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


		spinseed = (int)params.fget("Spin seed");
		dilutionseed = (int)params.fget("Dilution seed");
		biasseed =(int)params.fget("Bias seed");
		percent = params.fget("Diluted Percentage");
		biaspercent = params.fget("Bias percent");
		type = (int)params.fget("Initialization type");
		field = params.fget("Field");
		
		interventionstart=(int)params.fget("intervention start");
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
		MCS(isingspin, spinfliprand);
		params.set("MC time", prestep-50);
		Job.animate();
		}
		
		params.set("spinfliprand", spinflipseed);
		params.set("dilutionrand",dilutionrand.nextDouble());
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

			
			if(step-quenchstart==interventionstart)
			{
				for (int h=0; h<M; h++)
				{
				interventioncopy[h]=isingspin[h];
				}
				intervention(interventioncopies, interventionsteplimit, interventionP);
				steplimit=-1;
				movie(grid3, 8888, 8888);
			}
			MCS(isingspin, spinfliprand);
			params.set("MC time", step);
			Job.animate();
			params.set("Metastable state time", step-quenchstart);
			//movie(grid1, step, 0000);
			//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/biasdata/nucleation/data50.txt",step-quenchstart, (magnetization/Ms));
		}
		
		
	
	
	
	
    }// the end of run()
	
}