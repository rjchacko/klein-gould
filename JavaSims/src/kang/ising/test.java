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

public class test extends Simulation{
	
	Grid grid1 = new Grid("Grid 1");
	Grid grid2 = new Grid("Grid 2");
	
	public int L1, L2, M; //parameters for the lattice
	public int int1, int2;
    public double double1, double2;
    public double Dres1,Dres2;     //results with double format
    public int Ires1,Ires2;       //results with int format
    
	public int isingspin[];     //the array of the data
	public int initialcopy[];   //the array of the initial copy of the system
	public double dilutionmap[];
    
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

	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
	
	}
	
	public static void main (String[] test){
		new Control(new test(), "Kang Liu's test program" );
	}
	
	public void load(Control test){
		test.frame (grid1);
		test.frame (grid2);

		
		
		params.add("lattice's width", 200);
		params.add("lattice's length", 200);
		params.add("double1", new DoubleValue(100.5,0,1000).withSlider());
		params.add("double2", new DoubleValue(0, 0, 1).withSlider());
		
		params.add("int1", 100);
		params.add("int2", -4);
		params.add("Ires1");
		params.add("Ires2");
		params.add("Dres1");
		params.add("Dres2");
	
	}
	
    public void run(){
		
		
		int1 = (int)params.fget("int1");
		int2 = (int)params.fget("int2");

		double1=params.fget("double1");
		double2=params.fget("double2");

		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");

		M = L1 * L2 ;

		isingspin = new int[M];
		initialcopy = new int[M];
		dilutionmap = new double[M];
		double tt=0;
		
		for (int t=0; t<M; t++)
			{
			tt=t;
			dilutionmap[t]=tt;
			isingspin[t]=t%5-2;
			}
		
		Dres1=double1-int1;
		
		params.set("Ires1", Ires1);
		params.set("Ires2", Ires2);
		params.set("Dres1", Dres1);
		params.set("Dres2", Dres2);
		
		Job.animate();
		
    }
	
}