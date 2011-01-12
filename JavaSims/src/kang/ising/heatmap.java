package kang.ising;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import chris.util.PrintUtil;
import chris.util.Random;


import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DoubleValue;

public class heatmap extends Simulation{
	
	private DecimalFormat fmt = new DecimalFormat("0000000");
	
	Grid grid1 = new Grid("Ising lattice 2d for nucleation");
	Grid grid2 = new Grid("Heat map for dilution");
	Grid grid3 = new Grid("Nucleation events distribution");
	Grid grid4 = new Grid("Nucleating droplet");
	

	

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
	
	public double threshold;
	public int copies;
	

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
	public double dilutionmap[];     //the array for reset
	public double nucleationevents[];
	public double nodeadsites[];

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
	
	public void fieldset(double finalfield)
	{
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
		
		ColorGradient heatmap = new ColorGradient();
		
		
		grid1.setColors(ising);
		grid1.registerData(L1, L2, isingspin);
		grid2.setColors(heatmap);
		grid2.registerData(L1, L2, dilutionmap);
		grid3.setColors(heatmap);
		grid3.registerData(L1, L2, nucleationevents);
		grid4.setColors(heatmap);
		grid4.registerData(L1, L2, nodeadsites);
		
	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
		grid3.clear();
		grid4.clear();

	}
	
	public int X(int bx)
	{
		int realx=bx;
		if (bx>L1)
			realx=bx-L1;
		if (bx<0)
			realx=bx+L1;
		return realx;
	}
	
	public int Y(int by)
	{
		int realy=by;
		if (by>L2)
			realy=by-L2;
		if (by<0)
			realy=by+L2;
		return realy;
	}
	
	public void square(int label[], int range, int cx, int cy)  //the function to draw a square
	{
		int bx, by;
		int x, y;
		for(int t=0; t<M; t++)
			label[t]=0;
		for(bx=cx-range; bx<cx+range; bx++)
			for(by=cy-range; by<cy+range; by++)
			{
				x=X(bx);
				y=Y(by);
				label[x*L2+y]=1;
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
	
	public void initialize(int spin[], int type, double p, int spinseed, int biasR, int dilutionseed, double biasp, int biasseed)
	{
		spinrand= new Random(spinseed);
		dilutionrand= new Random(dilutionseed);
		biasrand= new Random(biasseed);
		deadsites=0;
		
		int biasrangelabel[];
		biasrangelabel= new int [M];    // the array to label the location of bias diluted site(1 for bias diluted site)
		
		for (int b=0;b<M;b++)
			biasrangelabel[b]=0;
			
		
		int cx, cy; //x,y index for the center
		cx=(int)L1/2;
		cy=(int)L2/2;
		if(biasp!=p)
			square(biasrangelabel, biasR, cx, cy);
		
		for(int t=0; t<M; t++)
			spin[t]=1;
		
		
		for(int j=0; j<M; j++)
		{
			if (biasrangelabel[j]==1)
				if(biasrand.nextDouble()<biasp)
					spin[j]=0;
			if (biasrangelabel[j]==0)
				if(dilutionrand.nextDouble()<p)
					spin[j]=0;
		}
		
		if(type==0)
			for (i=0; i<M; i++)
		    {
			   if(spin[i]!=0)
			   {
				   spin[i]=-1;
			   if (spinrand.nextDouble()> 0.5)
				   spin[i]=1;
			   }
				
		    }
		if(type==1)
			for (i=0; i<M; i++)
		    {
			   if(spin[i]!=0)
				   spin[i]=1;
		    }
		if(type==2)
			for (i=0; i<M; i++)
		    {
			   if(spin[i]!=0)
				   spin[i]=-1;
		    }
		
		for (i=0; i<M; i++){

			initialcopy[i]=spin[i];  // here, make the copy of the system
			
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
	
	public void dilutionheatmap(int spin[], int r)
	{
		double map[];
		map= new double[M];
		for(int s=0; s<M; s++)
		{
			map[s]=dilutionratio(spin, r, s, M);
			dilutionmap[s]=map[s];
		}
		
	}
	
	public void nucleationmap(int spin[], double events[], int copies, double threshold, double dilutionp,int dilutionseed, double biasp, int biasseed, int biasR)
	{
		int totalspin[]= new int[M];
		for(int y=0;y<M;y++)
			totalspin[y]=0;
		for(int ir=0; ir<copies; ir++)
		{
			params.set("replica",ir+1);
			initialize(spin, 0, dilutionp, ir+1, biasR, 1, biasp, 1);
			spinfliprand= new Random (ir+1);
			temperaturequench(9);
			params.set("Field", field);
			
			for (prestep=0; prestep < 50; prestep++)
			{
			MCS(isingspin, spinfliprand ,1);
			params.set("MC time", prestep-50);
			Job.animate();
			}
			temperaturequench(QuenchT);
			for(int bquench=0; bquench<100; bquench++)
			{
			MCS(isingspin, spinfliprand ,1);
			params.set("MC time", bquench);
			Job.animate();
			}
			Ms=magnetization;
			
			fieldquench();
			int nucleationend=-1;
			for(int j=0; nucleationend<0; j++)
			{
				MCS(isingspin, spinfliprand, 1);
				params.set("MC time", j);
				Job.animate();
				if(magnetization<(Ms*threshold))
					{
					nucleationend=1;
					params.set("Metastable state time",j);
					}
			}
			for(int x=0;x<M;x++)
			{
				totalspin[x]+=isingspin[x];
			}
		
			double scale=1;
			for (int xx=0; xx<M; xx++)
			{
				if(totalspin[xx]<scale)
					scale=totalspin[xx];
			}
			
			for (int yy=0; yy<M;yy++)
			{
				events[yy]=totalspin[yy]/scale;
			}
			
			Job.animate();
            movie(grid3, ir+1, ir+1);

		
		
		
		
		}
		

		
		
		
	}
	
	
	public void criticalquenchmap(int spin[], double events[], int copies, double threshold, double dilutionp,int dilutionseed, double biasp, int biasseed, int biasR)
	{
		int totalspin[]= new int[M];
		//int skipthisrun=0;
		int direction=1;
		for(int y=0;y<M;y++)
			totalspin[y]=0;
		for(int ir=0; ir<copies; ir++)
		{
			//skipthisrun=0;
			params.set("replica",ir+1);
			initialize(spin, 0, dilutionp, ir+1, biasR, 1, biasp, 1);
			spinfliprand= new Random (ir+1);
			temperaturequench(9);
			params.set("Field", 0);
			
			for (prestep=0; prestep < 20; prestep++)
			{
			MCS(isingspin, spinfliprand ,1);
			params.set("MC time", prestep-20);
			Job.animate();
			}
			temperaturequench(QuenchT);

			//int nucleationend=-1;
			for(int j=0; j<2; j++)
			{
				MCS(isingspin, spinfliprand, 1);
				params.set("MC time", j);
				Job.animate();

			}
			
			
			if(magnetization>0)
			{
			direction=-1;
			}
			if(magnetization<0)
			{
			direction=1;
			}
			
			for(int x=0;x<M;x++)
			{
				totalspin[x]+=(isingspin[x]*direction);
			}
		
			double scale=1;
			for (int xx=0; xx<M; xx++)
			{
				if(totalspin[xx]<scale)
					scale=totalspin[xx];
			}
			
			for (int yy=0; yy<M;yy++)
			{
				events[yy]=totalspin[yy]/scale;
			}
			
			Job.animate();
            movie(grid3, ir+1, ir+1);
		
		
		
		}
		

		
		
		
	}

	public void offcriticalquenchmap(int spin[], double events[], int copies, double threshold, double dilutionp,int dilutionseed, double biasp, int biasseed, int biasR, double hi, double hf)
	{
		int totalspin[]= new int[M];

		for(int y=0;y<M;y++)
			totalspin[y]=0;
		for(int ir=0; ir<copies; ir++)
		{
			//skipthisrun=0;
			params.set("replica",ir+1);
			initialize(spin, 0, dilutionp, ir+1, biasR, 1, biasp, 1);
			spinfliprand= new Random (ir+1);
			temperaturequench(9);
			params.set("Field",hi);
			int direction=1;
			if(hf>0)
				direction=-1;
			
			for (prestep=0; prestep < 20; prestep++)
			{
			MCS(isingspin, spinfliprand ,1);
			params.set("MC time", prestep-20);
			Job.animate();
			}
			temperaturequench(QuenchT);
			params.set("Field",hf);

			int nucleationend=-1;
			for(int j=0; nucleationend<0; j++)
			{
				MCS(isingspin, spinfliprand, 1);
				params.set("MC time", j);
				Job.animate();
				if(-magnetization*direction>(threshold*(1-dilutionp)))
					nucleationend=1;

			}
			
			
			for(int x=0;x<M;x++)
			{
				totalspin[x]+=(isingspin[x]*direction);
			}
		
			double scale=1;
			for (int xx=0; xx<M; xx++)
			{
				if(totalspin[xx]<scale)
					scale=totalspin[xx];
			}
			
			for (int yy=0; yy<M;yy++)
			{
				events[yy]=totalspin[yy]/scale;
			}
			
			Job.animate();
            movie(grid3, ir+1, ir+1);
		
		
		
		}
		
	}
	
	public void takeoutzeros(double events[])
	{
		double totalevents=0;
		for(int ss=0; ss<M; ss++)
		{
			totalevents+=events[ss];
		}
		double avgevents=totalevents/(M);
		for(int dd=0; dd<M; dd++)
		{
			nodeadsites[dd]=events[dd];
			if(events[dd]==0)
				nodeadsites[dd]=avgevents;
			
		}
	}
	
	public void movie(Grid grid, int number, int copynumber)   //function to capture the grid
	{
		
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/heatmap/homogeneous/q=0.111/pic_"+fmt.format(copynumber)+"_"+fmt.format(number)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
	}
	
	public static void main (String[] heatmap){
		new Control(new heatmap(), "Kang Liu's heatmap for diulted ising model" );
	}
	
	public void load(Control heatmap){
		heatmap.frame (grid1);
		heatmap.frame (grid2);
		heatmap.frame (grid3);
		heatmap.frame (grid4);
		
		
		
		params.add("lattice's width", 200);
		params.add("lattice's length", 200);
		params.add("Diluted Percentage", new DoubleValue(0.111,0,1).withSlider());
		params.add("Bias percent", new DoubleValue(0.111, 0, 1).withSlider());
		
		params.addm("Interaction Constant", -4);
		params.add("Interaction range", 10);
		params.add("Bias range", 10);
		params.add("Monte Carlo step's limit", 1000000);
		params.add("copies",100);
		params.add("Threshold", new DoubleValue(0.6, 0, 1).withSlider());

		//params.add("Spin seed", 1);
		//params.add("Dilution seed", 1);
		//params.add("Bias seed", 1);
		//params.add("Initialization type", 0);
		params.addm("Quench temperature", new DoubleValue(1.591, 0, 10).withSlider());
		params.addm("Temperature", new DoubleValue(9, 0, 10).withSlider());
		params.addm("Field", new DoubleValue(-0.1, -5, 5).withSlider());
		
		
		
		///below is the parameters that would be displayed on the panel
		
		params.add("MC time");
		params.add("Metastable state time");
		params.add("magnetization");
		params.add("Saturated magnetization");
		params.add("lifetime");

		
        params.add("replica");
        
		
	}
	
    public void run(){
		
    	
		
		R = (int)params.fget("Interaction range");
		biasrange= (int)params.fget("Bias range");
	    steplimit = (int)params.fget("Monte Carlo step's limit");
		QuenchT=params.fget("Quench temperature");

		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");

		M = L1 * L2 ;

		
		isingspin = new int[M];
		initialcopy = new int[M];
		dilutionmap = new double [M];
		nucleationevents = new double[M];
		nodeadsites = new double [M];
		

		//spinseed = (int)params.fget("Spin seed");
		//dilutionseed = (int)params.fget("Dilution seed");
		//biasseed =(int)params.fget("Bias seed");
		dilutionseed=1;
		biasseed=1;
	
		percent = params.fget("Diluted Percentage");
		biaspercent = params.fget("Bias percent");
		//type = (int)params.fget("Initialization type");
		type=0;
		field = params.fget("Field");
		
		initialize(isingspin, 0, percent, 1, R, 1, biaspercent, 1);
		dilutionheatmap(initialcopy, R);
		Job.animate();
		movie(grid2, 9999, 9999);
		
		copies= (int)params.fget("copies");
		threshold= params.fget("Threshold"); 
		
		
		//nucleationmap(isingspin, nucleationevents, copies, threshold, percent,1, biaspercent, 1, R);
		//criticalquenchmap(isingspin, nucleationevents, copies, threshold, percent,1, biaspercent, 1, R);
		offcriticalquenchmap(isingspin, nucleationevents, copies, threshold, percent,1, biaspercent, 1, R, field, -field);
		
		takeoutzeros(nucleationevents);
		Job.animate();
		movie(grid4,6666,6666);
		
	
    }
	
	
}