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

public class randomnumbertest extends Simulation{
	
	Grid grid1 = new Grid("Diluted Ising lattice 2d");
	Grid grid2 = new Grid("copy");
	
	public int L1, L2, M; //parameters for the lattice
	
	public int i,j,k,l,m;
	
	public int step;
	
	public int clonestep;
	
	public int isingspin[];     //the array of the data
	public int isingcopy[];
	
	private DecimalFormat fmt = new DecimalFormat("0000000");

	private Random rand;
	
	private int SEED1, SEED2;
	
	public static void main (String[] randomnumber){
		new Control(new randomnumbertest(), "Kang Liu's random number test program" );
	}
	
	public void load(Control kangliurandom){

		kangliurandom.frame (grid1);
		kangliurandom.frame (grid2);
		
		params.add("lattice's width", 100);
		params.add("lattice's length", 100);
		
		params.add("Seed 1", 1);
		params.add("Seed 2", 2);
		params.add("random number 1");
		params.add("random number 2");
		params.add("random number copy");
		params.add("step");
		params.add("clone step",10);
		
		
	}
	
    public void animate(){
		
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);
		ising.setColor(-1, Color.WHITE);
		ising.setColor(0, Color.RED);
		
		grid1.setColors(ising);
		grid2.setColors(ising);
		grid1.registerData(L1, L2, isingspin);
		grid2.registerData(L1, L2, isingcopy);

		
		
	}
	
	public void clear(){
		grid1.clear();
		grid2.clear();
	
	}
	
	public void run(){
		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");
		SEED1 =(int)params.fget("Seed 1");
		SEED2 =(int)params.fget("Seed 2");
		clonestep = (int)params.fget("clone step");
		
		Random prand;
		Random srand;
		
		M=L1 * L2;
		isingspin = new int[M];
		isingcopy = new int[M];
		
		for (i=0; i<M; i++){
			isingspin[i]=-1;
			isingcopy[i]=-1;
			//Job.animate();
		}
		

		
		prand = new Random(SEED2); // random number for probability
		srand = new Random(SEED1); // random number for spin position
		
		
		for (j=0; j<M; j++){
			k=(int) (srand.nextDouble()*M); 
			l=(int) (prand.nextDouble()*M);
			if (j==clonestep)
				rand=srand.clone();
			if (j>clonestep)
				m=(int) (rand.nextDouble()*M);
				
			isingspin[k]=-isingspin[k];
			isingcopy[l]=-isingcopy[l];
			Job.animate();
			params.set("random number 1", k);
			params.set("random number 2", l);
			params.set("random number copy", m);
			params.set("step", j);
		}
//		for(int jj = 0 ; jj < 100 ; jj++)
//			rand.nextDouble();
		//Random newr = rand.clone();
		
		
	}
	
	
}