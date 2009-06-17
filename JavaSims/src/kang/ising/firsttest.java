package kang.ising;// what is this used for?

import java.awt.Color;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

public class firsttest extends Simulation{
	Grid grid1 = new Grid("Ising lattice 2d");
	Grid grid2 = new Grid("Heisenberg lattice 2d");
	public int L1, L2,M;
	public double T;
	public double H;
	public int isingspin[];
	public double heisenbergspin[];
	
	//items below this line are the functions
	
	public static void main (String[] kangliu){
		new Control(new firsttest(), "Kang Liu's first try" );
	}
	
	public void load(Control liu){
		liu.frame (grid1);
		liu.frame (grid2);
		params.add("lattice's width", 12);
		params.add("lattice's length", 24);
		params.add("Temperature", new DoubleValue(300, 0, 800).withSlider());
		params.addm("Field", new DoubleValue(90, 0, 10000).withSlider());
		
		params.add("Model type", new ChoiceValue("noninteracting", "interacting"));
		
		
	}
	
	public void animate(){
		int i;
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);
		ising.setColor(-1, Color.WHITE);
		
		ColorGradient heisenberg = new ColorGradient();
		
		for (i=0; i < M; i++){
			heisenberg.getColor(heisenbergspin[i], -1., 1.);
		}
		grid1.setColors(ising);
		grid2.setColors(heisenberg);
		
		grid1.registerData(L1, L2, isingspin);
		grid2.registerData(L1, L2, heisenbergspin);
		
		
		
	}
	
	public void clear(){
		grid1.clear();
		grid2.clear();
	}
	
	public void run (){
		int i,j;
		T = params.fget("Temperature");
		H = params.fget("Field");
		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");
		M = L1 * L2;
		isingspin = new int[M];
		heisenbergspin = new double[M];
		
		for (i=0; i<M; i++){
			isingspin[i]=-1;
			if (Math.random()> 0.5)
				isingspin[i]=1;
			
		}
		
		for (j=0; j<M; j++){
			heisenbergspin[j]= Math.random();
		}
		
		Job.animate();
	}
}
