package kang.ising;

import java.awt.Color;

import chris.util.PrintUtil;
import chris.util.Random;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DoubleValue;

public class wanglandau extends Simulation{
	Grid grid1 = new Grid("Ising lattice 2d");

	public int L1, L2, M; //parameters for the lattice
	public int i, x, y;  //parameters for the index
	public int R; // interaction range
	public int steplimit;  //limit of mcs steps
	public int step;
	public int outputstep;
	public int wanglandaustep;
	public int N; // number of the beams in the energy space
	public int isingspin[];
	public int histogram[]; // array for the histogram
	public double densityG[]; // array for the density of states
	public double logdensityG[];  // array for the log of G
	public double T;     //temperature 
	public double H;     //field
	public double J;     //interaction constant after normalization
 	public double f;     //normalization factor
	public double P;     //dilution percentage
	
	public double NJ;    //interaction constant before normalization
	
	public double totalmagnetization;
	public double magnetization;

	public int Nneighber(int a,int i, int length1, int length2)// function for the index of nearest neighbor
	{
		int nx,ny; //index for neighbor
		int ni=0;
		nx=(int)i/length2;
		ny=(int)i%length2;
		
		if (a==0) {
			ni=nx*length2+ny-1;
			if  (ny==0) {
				ni=nx*length2+ny+length2-1;
			}
			
		}//(x,y-1) up
		
		if (a==1){
			ni=(nx+1)*length2+ny;
			if  (nx==length1-1) {
				ni=(nx-length1+1)*length2+ny;
			}
			
		}//(x+1,y) right
		
		if (a==2){
			ni=nx*length2+ny+1;
			if  (ny==length2-1) {
				ni=nx*length2+ny-length2+1;
			}
			
		}//(x,y+1) down
		
		if (a==3){
			ni=(nx-1)*length2+ny;
			if  (nx==0) {
				ni=(nx+length1-1)*length2+ny;
			}
		}//(x-1,y) left
		
		return ni;
		
	}
	
	public double interactionE (int j, int spin[], double J, int length1, int length2){ //function for interaction energy
		double Energy=0;
		int b,k;
		for(b=0; b<4;b++){
			k=Nneighber(b,j,length1, length2);
			Energy=Energy+J*spin[j]*spin[k];
		}
		return Energy;
			
	}
	
	public double totalE (int spin[], double J, int length1, int length2)
	{
		double energy=0;
		int i, spinnumber;
		spinnumber=length1*length2;
		for(i=0; i<spinnumber; i++)
			energy+=interactionE(i, spin, J, length1, length2);
		return energy/2;
		
	}
	
	public int beamposition(double energy)
	{
		int p;// the position in the energy space
		p=(int)((energy+2*M*NJ/4)/(M*NJ/N));
		return p;
	}
	
	
    public void wanglandauspinflip(int i, int spin[], int R, double T, double H, double f)
    {
		double ZeemanE1=-H*spin[i];// initial field energy
		double InterE1=0;
        double InterE2=0;
        double logf=Math.log(f);
		
		if (R==0) {
			J=NJ/4;
			InterE1=totalE(spin, J, L1, L2); // initial interaction energy
			InterE2=InterE1-2*interactionE(i,spin, J, L1, L2); //interaction energy after flip
		}
		
			
		double ZeemanE2=-H*(-spin[i]); //field energy after flipping
		
		double E1=ZeemanE1+InterE1;
		double E2=ZeemanE2+InterE2;
		
		int p1=beamposition(E1);
		int p2=beamposition(E2);
		
		if (logdensityG[p1]>=logdensityG[p2]){

			spin[i] = -spin[i];
			histogram[p2]++;
			//densityG[p2]=f*densityG[p2];
			logdensityG[p2]=logdensityG[p2]+logf;
			
			
		}// flip the spin
		
		
		else{
			if (Math.random()<(Math.exp(logdensityG[p1]-logdensityG[p2])))
			{
				spin[i]=-spin[i];
				histogram[p2]++;
				//densityG[p2]=f*densityG[p2];
				logdensityG[p2]=logdensityG[p2]+logf;
			
			}
			else
			{
				histogram[p1]++;
				//densityG[p1]=f*densityG[p1];
				logdensityG[p1]=logdensityG[p1]+logf;
			}
		}
    }
	
	public static void main (String[] wanglaudautest){
		new Control(new wanglandau(), "Kang Liu's wang-landau algorithm for NN ising model" );
	}
	
	public void load(Control liu){
		liu.frame (grid1);
		params.add("lattice's width", 100);
		params.add("lattice's length", 100);
		params.add("Number of the energy beams", 100);
		params.addm("Temperature", new DoubleValue(1, 0, 10).withSlider());
		params.addm("Field", new DoubleValue(0, -2, 2).withSlider());
		params.addm("Interaction Constant before normalization", -4);
		params.add("Interaction range", 0);
		params.add("Monte Carlo step's limit", 1000000);
		params.add("Output step",2000);
		params.add("Wanglandau step", 200);

		params.add("Diluted Percentage", new DoubleValue(0,0,1).withSlider());
		
		params.add("MC time");
		params.add("magnetization");
		
	}
	
	public void animate(){
		
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);
		ising.setColor(-1, Color.WHITE);
		ising.setColor(0, Color.RED);
		
		grid1.setColors(ising);
		grid1.registerData(L1, L2, isingspin);
		
		
	}
	
	public void clear(){
		grid1.clear();
	}
	
	public void run()
	{
		int i,j,k;
		R = (int)params.fget("Interaction range");
		steplimit = (int)params.fget("Monte Carlo step's limit");
		outputstep = (int)params.fget("Output step");
		wanglandaustep = (int)params.fget("Wanglandau step");
		P = params.fget("Diluted Percentage");
		N = (int) params.fget("Number of the energy beams");
		f = Math.exp(1);

		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");
		M = L1 * L2;
		
		isingspin= new int[M];
		histogram= new int[N+1];
		densityG= new double[N+1];
		logdensityG= new double[N+1];
		
		//randomize the initial state
		for (i=0; i<M; i++){
			isingspin[i]=-1;
			if (Math.random()> 0.5)
				isingspin[i]=1;
		}
		
		for (i=0; i<M; i++){
			if(Math.random()< P)
				isingspin[i]=0;
		}// here, the sequence of dilution and initialization is not important
		
		for(k=0; k<N+1; k++){
			histogram[k]=0;
			densityG[k]=1;
			logdensityG[k]=0;
		
		}
		
		
		
		for (step=0; step< steplimit; step++){
			
			for (int z=0; z< M; z++){
			
			    H = params.fget("Field");
			    T = params.fget("Temperature");
			    NJ = params.fget("Interaction Constant before normalization");
			
				j=(int) (Math.random()*M); //choose a spin randomly
				
				wanglandauspinflip(j, isingspin, R, T, H, f);
			
				

				Job.animate();
				
			}
			
			if(step % 1 == 0) {
				params.set("MC time", step);
			}
			
			totalmagnetization=0;	
			for(int s=0; s<M; s++)
			{
				totalmagnetization+=isingspin[s];
			}
			magnetization=totalmagnetization/M;
			params.set("magnetization", magnetization);
			
			if(step % wanglandaustep ==0)
			{
				for(int d=0; d<N+1; d++)
					histogram[d]=0;
				f=Math.sqrt(f);
			}
			
			
			
			if(step==outputstep-1){
				for(int b=0; b<N+1; b++)
					{
					densityG[b]=Math.exp(logdensityG[b]-logdensityG[N/2]);
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/histogram.txt",b, histogram[b]);
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/logdensity.txt",b, logdensityG[b]);
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/density.txt",b, densityG[b]);
					}
			}
			
			
			
			
		}

		
		
	}
	
	
}