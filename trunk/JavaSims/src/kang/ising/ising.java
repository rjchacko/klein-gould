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
	public double J;     //interaction constant
	public int R;   //interaction range
	public int step; //Monte Carlo step
	public int steplimit; //upperlimit of MC step
	
	
	public int isingspin[];     //the array of the data
	
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
			System.out.print("neighber=  ");
			System.out.println(k);
			Energy=Energy+J*isingspin[j]*isingspin[k];
		}
		return Energy;
			
	}
	
	
	/*public double longrangeE (int j){
		double Energy=0;
		int a;
	}*/
	
	
	public static void main (String[] kangliu){
		new Control(new ising(), "Kang Liu's ising model" );
	}
	
	public void load(Control liu){
		liu.frame (grid1);
		params.add("lattice's width", 12);
		params.add("lattice's length", 24);
		params.add("Temperature", new DoubleValue(5, 0, 100).withSlider());
		params.addm("Field", new DoubleValue(5, 0, 500).withSlider());
		params.add("Interaction Constant", 5);
		params.add("interaction range", 0);
		params.add("Monte Carlo step's limit", 100);
		
		
		params.add("Model type", new ChoiceValue("noninteracting", "interacting"));
		params.add("Mechanics", new ChoiceValue("Metropolis", "Kawasaki"));
		
	}
	
	public void animate(){
		
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);
		ising.setColor(-1, Color.WHITE);
		
		grid1.setColors(ising);
		grid1.registerData(L1, L2, isingspin);
		
		
	}
	
	public void clear(){
		grid1.clear();
		
	}
	
	public void run(){
		
		int i,j;
		T = params.fget("Temperature");
		H = params.fget("Field");
		J = params.fget("Interaction Constant");
		steplimit = (int)params.fget("Monte Carlo step's limit");
		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");
		M = L1 * L2;
		isingspin = new int[M];
		
		
		//randomize the initial state
		for (i=0; i<M; i++){
			isingspin[i]=-1;
			if (Math.random()> 0.5)
				isingspin[i]=1;
				//System.out.print("spin = ");
				//System.out.println(isingspin[i]);
				//PrintUtil.printlnToFile("/Users/cserino/Desktop/foo.txt","spin = ", isingspin[i]);
		}
		
		
		for (step=0; step< steplimit; step++){
				j=(int) (Math.random()*M); //choose a spin randomly
				
				System.out.println("j=");
				System.out.println(j);
			
				double ZeemanE1=-H*isingspin[j];// initial field energy
				double InterE1=interactionE(j); // initial interaction energy
				
				
				double ZeemanE2=-H*(-isingspin[j]); //field energy after flipping
				double InterE2=-interactionE(j);
				
				
				double E1=ZeemanE1+InterE1;
				double E2=ZeemanE2+InterE2;
				
				if (E1>E2){
					int tempS=isingspin[j];
					isingspin[j]=-tempS;
				}// flip the spin
				
				if (E1<E2){
					if (Math.random()<Math.exp((E1-E2)/T)){
						int tempS=isingspin[j];
						isingspin[j]=-tempS;
					
					}
				}	
				
				Job.animate();
				
			}
		
		
			
				
	}
}
	
	
