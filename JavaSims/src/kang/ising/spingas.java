package kang.ising;

import java.awt.Color;

import chris.util.PrintUtil;

import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

public class spingas extends Simulation{
	Grid grid1 = new Grid("Spin gas lattice 2d");
	public int L1, L2,M; //parameters for the lattice
	public int i, x, y;  //parameters for the index
	public double T;     //temperature 
	public double chp;     //chemical potential
	public double K;     //interaction constant after normalization
	public double NK;    //interaction constant before normalization
	
	public int R;   //interaction range
	public int step; //Monte Carlo step
	public int steplimit; //upperlimit of MC step
	public int metricstart; //the start of metric calculation
	public int N; //how many points in the metric plot
	public double P; //dilution percentage
	
	
	
	public int isingspin[];     //the array of the data
	public int particleposition[]; //tracking the particles
	
	
	public double timetotalE[];  
	public double timeaverageE[]; //time average energy for metric
	public double Metric[];  //array for plotting metric
	
	
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
	
	public double interactionE (int j){ //function for interaction energy
		double Energy=0;
		int b,k;
		for(b=0; b<4;b++){
			k=Nneighber(b,j);
			if (isingspin[k]!=0)
				Energy=Energy+K*(isingspin[j]+1)*(isingspin[k]+1)/4;
	
		}
		return Energy;	
	}
	
	public double longrangeE (int i){
		double Energy=0;
		int totalN=0;
		int nx=i/L2;
		int ny=i%L2;
		int kx, ky;
		
		for (int m=-R; m<=R; m++)
			for (int n=-R; n<=R; n++)
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
				if(isingspin[kx*L2+ky]!=0)
					totalN+=(isingspin[kx*L2+ky]+1)/2;
			}
		Energy=K*(totalN-(isingspin[i]+1)/2)*(isingspin[i]+1)/2;
		return Energy;
	}
	
	
	public static void main (String[] kangliu){
		new Control(new spingas(), "Kang Liu's spin gas model" );
	}
	
	
	public void load(Control kang){
		kang.frame (grid1);
		params.add("lattice's width", 100);
		params.add("lattice's length", 100);
		params.addm("Temperature", new DoubleValue(1, 0, 100).withSlider());
		params.addm("Chemical potential", new DoubleValue(0, -2, 2).withSlider());
		params.addm("Interaction Constant", 1);
		params.add("Interaction range", 10);
		params.add("Monte Carlo step's limit", 1000000);
		params.add("Metric Start at",5000);
		params.add("Metric points", 2000);
		params.add("Diluted Percentage", new DoubleValue(0,0,1).withSlider());
		
		//params.add("Model", new ChoiceValue("noninteracting", "interacting"));
		//params.add("Mechanics", new ChoiceValue("Metropolis", "Kawasaki"));
		params.add("MC time");
		//params.add("magnetization");
		
	}
	
	
    public void animate(){
		
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);   //particles
		ising.setColor(-1, Color.WHITE);    //holes
		ising.setColor(0, Color.RED);
		
		grid1.setColors(ising);
		grid1.registerData(L1, L2, isingspin);
		
		
	}
    
    
    public void clear(){
		grid1.clear();
		
	}
	
    public void run(){
    	int i,j;
		R = (int)params.fget("Interaction range");
		steplimit = (int)params.fget("Monte Carlo step's limit");
		metricstart = (int)params.fget("Metric Start at");
		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");
		M = L1 * L2;
		isingspin = new int[M];
		timetotalE = new double[M];
		timeaverageE = new double[M];
		N = (int)params.fget("Metric points");
		Metric = new double[N];
		
		P = params.fget("Diluted Percentage");
		
		double NMetric=0;
		double totalE=0;
		double averageE=0;// definition for metric variables
		
		int X=0;  //the number of the particles in the spin gas
		
		//randomize the initial state
		for (i=0; i<M; i++){
			if (Math.random()< P)
				isingspin[i]=0;               //dilute the lattice first
			if (isingspin[i]!=0)              //if the site is not dead
			isingspin[i]=-1;
			if (Math.random()> 0.5)
				{
				isingspin[i]=1;
				X+=1;                         //count the number of the particles
				}
				
		}
		
		
		// track all the particles
		
		particleposition = new int[X];
		j=0;                     //initialize the index for particleposition array
		for (i=0; i<M; i++){
			if (isingspin[i]==1)
			{
				particleposition[j]=i;         //record the position of the jth particle
				j++;
			}
			
		}
    	
    	//enter the MCS
		
		for (step=0; step< steplimit; step++){
		
			
			
			
			
			
			
		}
		
		
		
    	
    }
	
}