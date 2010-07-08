package kang.ising;

import java.awt.Color;

//import java.text.Decima



import chris.util.PrintUtil;
import chris.util.Random;


import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DoubleValue;

public class ising3d extends Simulation{
	
//	private DecimalFormat fmt = new DecimalFormat("0000000");
	
	Grid grid1 = new Grid("Diluted Ising lattice 3d for nucleation");
	Grid grid2 = new Grid("Diluted Ising 3d initial copy");


	public int L1, L2, L3, M; //parameters for the lattice
	public int layer; //index for layer (nx)
    public int i, x, y, z;  //parameters for the index
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
//	public int lifetime;   // the time starts at the moment when quench of the fields happens
	public int quenchstart; // the time of the quench

	public int deadsites;
	
	public int prestep;  // step for preparation of intial configuration
	public int step; //Monte Carlo step
	public int steplimit; //upperlimit of MC step

	public int isingspin[];     //the array of the data
	public int initialcopy[];   //the array of the initial copy of the system
	public int isinglayer[];
	public int copylayer[];

	public int initialtype;     //the type of initialization 0-random 1-all up 2-all down
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
	

	public int Xstart;
	public int Xend;
	public int Mend;  // the end of recording magnetization
	public double X;
	
	public int totallifetime;
	public int meanlifetime;
	
	
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
	
	public double interactionEchange (int range, double constJ, int spin[], int j){ //function for interaction energy
		double Energy=0;
		double Energychange=0;
		if (range==0)
		{
			J= constJ/6;                     // normalize the interaction constant
			int b,k;
		    for(b=0; b<6;b++)
		    {
			k=Nneighber(b,j);
			Energy=Energy+J*spin[j]*spin[k];
		    }
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
	
	
	public void animate(){
		
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);
		ising.setColor(-1, Color.WHITE);
		ising.setColor(0, Color.RED);
		
		for(int w=0; w<(L2*L3); w++)
		{
			isinglayer[w]=isingspin[(layer*L2*L3)+w];
			copylayer[w]=initialcopy[(layer*L2*L3)+w];
		}
		
		grid1.setColors(ising);
		grid1.registerData(L2, L3, isinglayer);
		grid2.setColors(ising);
		grid2.registerData(L2, L3, copylayer);

		
	}
	public void clear(){
		grid1.clear();
		grid2.clear();

	}
	
	public static void main (String[] ising3d){
		new Control(new ising3d(), "Kang Liu's diulted ising model for nucleation" );
	}
	
	public void load(Control ising3d){
		ising3d.frame (grid1);
		ising3d.frame (grid2);

		
		params.add("lattice's width", 50);
		params.add("lattice's length", 50);
		params.add("lattice's height", 50);
		params.add("Diluted Percentage", new DoubleValue(0.40,0,1).withSlider());

		params.add("Quench starts at", 200);

		params.add("Mag ends", 500);
		params.add("Xstart",100);
		params.add("Xend",400);


		params.addm("Interaction Constant before normalization", -6);
		params.add("Interaction range", 0);
		
		params.add("Monte Carlo step's limit", 1000);
		//params.add("Position seed", 1001);
		params.add("Spinflip seed", 1001);
		params.add("Spin seed", 1);
		params.add("Dilution seed", 1);
		params.add("Initialization type", 0);
		
		params.addm("Quench temperature", new DoubleValue(0.8073, 0, 10).withSlider());
		params.addm("Temperature", new DoubleValue(9, 0, 10).withSlider());
		params.addm("Field", new DoubleValue(0.65, -5, 5).withSlider());
		params.add("MC time");
		params.add("Metastable state time");
		params.add("mean life time");
		params.add("magnetization");
		params.add("Saturated magnetization");
		

		params.add("u#copies");
	    params.add("v#fieldruns");

        params.add("spinfliprand");
    
        params.add("dilutionrand");
        params.add("susceptibility");

		
	}
	
	public void run(){
		
		
		R = (int)params.fget("Interaction range");

		steplimit = (int)params.fget("Monte Carlo step's limit");


		quenchstart=(int)params.fget("Quench starts at");
		QuenchT=params.fget("Quench temperature");

		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");
		L3 =(int)params.fget("lattice's height");
		M = L1 * L2 * L3;
		Mend=(int)params.fget("Mag ends");
		
		isingspin = new int[M];
		initialcopy = new int[M];
		isinglayer = new int[L2*L3];
		copylayer =new int[L2*L3];
		layer = (int)L2/2;

		//Mag = new double [Mend];
		//Xstart= (int)params.fget("Xstart");
		//Xend= (int)params.fget("Xend");

		spinseed = (int)params.fget("Spin seed");
		dilutionseed = (int)params.fget("Dilution seed");
		
		percent = params.fget("Diluted Percentage");
		type = (int)params.fget("Initialization type");
		field = params.fget("Field");
		//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/3disingdata/lowtemp/q=0.40.txt","quench temperature=   ", QuenchT);

		for(int v=0; v<60; v++)
		{
		
		initialize(isingspin, type, percent, spinseed, dilutionseed);
		
		Job.animate();
		
		totallifetime=0;
		params.set("v#fieldruns", v);

		for(int u=0; u<30; u++)
		{
			for(int r=0; r<M; r++)
			{
				isingspin[r]=initialcopy[r];
			}
			params.set("Temperature",9);	//reset the high temperature everytime
			params.set("Field",field-0.01*v);
			params.set("Quench temperature",0.8073);
			QuenchT=params.fget("Quench temperature");
			spinflipseed = (int)((field+1+v)*QuenchT*10000+u);
			spinfliprand= new Random (spinflipseed);
	
		for (prestep=0; prestep < 50; prestep++)
			{
			MCS(isingspin, spinfliprand);
			params.set("MC time", prestep-50);
			Job.animate();
			}
		
		params.set("spinfliprand", spinflipseed);
		params.set("u#copies", u);
		params.set("dilutionrand",dilutionrand.nextDouble());
		temperaturequench(QuenchT);
		steplimit=1000;
	
		
		for (step=0; steplimit>0; step++)
		{

			if(step==quenchstart)
				{
				Ms=magnetization;
				params.set("Saturated magnetization", Ms);
				fieldquench();
				}
			if((magnetization<(0.9*Ms))&(step>quenchstart))
			{
				totallifetime+=(step-quenchstart);
				steplimit=-1;// set steplimit so that the MCS ends
				
			}
			MCS(isingspin, spinfliprand);
			params.set("MC time", step);
			Job.animate();
			params.set("Metastable state time", step-quenchstart);

		}
		
		meanlifetime=(int)(totallifetime/(u+1));
		params.set("mean life time",meanlifetime);
		
	}// end of u circle
		PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/3disingdata/lowtemp/q=0.40.txt",-H, meanlifetime);
		
	}// end of v circle(the sweep of field)
		
	}// end of run()
	
	
}