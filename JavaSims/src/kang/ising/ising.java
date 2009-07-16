package kang.ising;// what is this used for?
/**
 * Kang, this is used to organize things. So, for example, if you open 
 * one of my files, for example ergodicDiagramApp then you will see that the package
 * basically tells JAVA / Eclipse what directory and subdirectory to store the file in. 
 * 
 * Also, you can have public, private, and protected variables and classed e.g.
 * private double x = 0
 * public double x = 0
 * protected double = 0
 * 
 * and protected numbers can be accessed from any class in the same package as the 
 * protected variable / class.
 * 
 *  -- Chris
 * 
 */

import java.awt.Color;

import chris.util.PrintUtil;

import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

public class ising extends Simulation{
	Grid grid1 = new Grid("Ising lattice 2d");
	public int L1, L2,M; //parameters for the lattice
	public int i, x, y;  //parameters for the index
	public double T;     //temperature 
	public double H;     //field
	public double J;     //interaction constant after normalization
	public double NJ;    //interaction constant before normalization
	public double totalmagnetization;
	public double magnetization;
	public int R;   //interaction range
	public int step; //Monte Carlo step
	public int steplimit; //upperlimit of MC step
	public int metricstart; //the start of metric calculation
	public int N; //how many points in the metric plot
	public double P;
	
	
	
	public int isingspin[];     //the array of the data
	
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
//			System.out.print("neighber=  ");
//			System.out.println(k);
			Energy=Energy+J*isingspin[j]*isingspin[k];
		}
		return Energy;
			
	}
	

	
	
	public double longrangeE (int i){
		double Energy=0;
		int S=0;
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
				S+=isingspin[kx*L2+ky];	
			}
		Energy=J*isingspin[i]*S-J;
		return Energy;
	}
	
	
	public static void main (String[] kangliu){
		new Control(new ising(), "Kang Liu's ising model" );
	}
	
	public void load(Control liu){
		liu.frame (grid1);
		params.add("lattice's width", 100);
		params.add("lattice's length", 100);
		params.addm("Temperature", new DoubleValue(1, 0, 100).withSlider());
		params.addm("Field", new DoubleValue(0, -2, 2).withSlider());
		params.addm("Interaction Constant", 1);
		params.add("Interaction range", 10);
		params.add("Monte Carlo step's limit", 1000000);
		params.add("Metric Start at",5000);
		params.add("Metric points", 2000);
		params.add("Diluted Percentage", new DoubleValue(0,0,1).withSlider());
		
		//params.add("Model", new ChoiceValue("noninteracting", "interacting"));
		//params.add("Mechanics", new ChoiceValue("Metropolis", "Kawasaki"));
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
		
		
		
		//randomize the initial state
		for (i=0; i<M; i++){
			isingspin[i]=-1;
			if (Math.random()> 0.5)
				isingspin[i]=1;
				//System.out.print("spin = ");
				//System.out.println(isingspin[i]);
				//PrintUtil.printlnToFile("/Users/cserino/Desktop/foo.txt","spin = ", isingspin[i]);
		}
		
		for (i=0; i<M; i++){
			if(Math.random()< P)
				isingspin[i]=0;
		}// here, the sequence of dilution and initialization is not important
		
		
		for (step=0; step< steplimit; step++){
			
			for (int f=0; f< M; f++){
			
			    H = params.fget("Field");
			    T = params.fget("Temperature");
			    NJ = params.fget("Interaction Constant");
			
				j=(int) (Math.random()*M); //choose a spin randomly
				
//				System.out.print("Step=");
//				System.out.println(step);
//				System.out.print("j=");
//				System.out.println(j);
				
			
				double ZeemanE1=-H*isingspin[j];// initial field energy
				double InterE1=0;
                double InterE2=0;
				
				if (R==0) {
					J=NJ/4;
					InterE1=interactionE(j); // initial interaction energy
					InterE2=-interactionE(j);
				}
				
				if (R!=0) {
					J=NJ/((2*R+1)*(2*R+1)-1);
					InterE1=longrangeE(j);
					InterE2=-longrangeE(j)-2*J;
				}
				
				
				
				
				double ZeemanE2=-H*(-isingspin[j]); //field energy after flipping
				
				
				
				double E1=ZeemanE1+InterE1;
				double E2=ZeemanE2+InterE2;
				
				if (E1>E2){
//					int tempS=isingspin[j];
//					isingspin[j]=-tempS;
					isingspin[j] = -isingspin[j];
				}// flip the spin
				
				if (E1<E2){
					if (Math.random()<Math.exp((E1-E2)/T)){
						//int tempS=isingspin[j];
						isingspin[j]=-isingspin[j];
					
					}
				}
			
				if(step % 10 == 0) params.set("MC time", step);
				Job.animate();
			}
			
			if(step > metricstart){
			
			totalE=0;
			
			for (int y=0; y<M; y++)
			{
				if(R==0)
					timetotalE[y]+=interactionE(y);
				if(R!=0)
					timetotalE[y]+=longrangeE(y);
				
				timeaverageE[y]=timetotalE[y]/(step-metricstart);
			 
			}
			
			for (int x=0; x< M; x++)
			{
					totalE+=timeaverageE[x];
			}
			
			averageE=totalE/M;

			NMetric=0;
			for (int z=0; z< M; z++)
			{
				NMetric+=(timeaverageE[z]-averageE)*(timeaverageE[z]-averageE);
			}
			
		
			if(step<N+metricstart+1)
				{
				Metric[step-metricstart-1]=NMetric/M;
				PrintUtil.printlnToFile("F:/data/dmetric6.txt",step-metricstart, NMetric/M);
				//PrintUtil.printlnToFile("/Users/cserino/Desktop/metric2.txt",step-metricstart, NMetric/M);
				}
			
			if(step==N+metricstart+1)
			{
				PrintUtil.printlnToFile("F:/data/dmetric1.txt","");
				PrintUtil.printlnToFile("F:/data/dmetric1.txt","Lattice length=", L1);
				PrintUtil.printlnToFile("F:/data/dmetric1.txt","Lattice width=", L2);
				PrintUtil.printlnToFile("F:/data/dmetric1.txt","J=",NJ);
				PrintUtil.printlnToFile("F:/data/dmetric1.txt","Interaction Range=",R);
				PrintUtil.printlnToFile("F:/data/dmetric1.txt","Mectric starts at", metricstart);
				PrintUtil.printlnToFile("F:/data/dmetric1.txt","Final Tempearture=",T);
				PrintUtil.printlnToFile("F:/data/dmetric1.txt","Final Field=",H);
				PrintUtil.printlnToFile("F:/data/dmetric1.txt","Diluted Percentage=",P);
				
			}
			
			}
			
		totalmagnetization=0;	
		for(int s=0; s<M; s++)
		{
			totalmagnetization+=isingspin[s];
		}
		magnetization=totalmagnetization/M;
		params.set("magnetization", magnetization);
		}
		
			
				
	}
}
	
	
