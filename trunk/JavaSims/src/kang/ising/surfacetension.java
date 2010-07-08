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

public class surfacetension extends Simulation{
	
	Grid grid1 = new Grid("Diluted Ising lattice 2d surface");
	public int L1, L2, M; //parameters for the lattice
	public int i, x, y;  //parameters for the index
	public double percent;  //diluted percent
	public int deadsites;
	public int numberofcopies;
	public int isingspin[]; 
	public int bonds;
	public int totalbonds;
	public double bondsratio;
	public double totalbondsratio;
	public int averagebonds;
	
	public Random dilutionrand;
	public int dilutionseed;
	
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
	
	public void initialize(int spin[], double p, Random rand, int totalspins)
	{
		deadsites=0;
		for(int j=0; j<totalspins; j++)
			spin[j]=1;
		for(int k=0; k<totalspins; k++)
		{
			if(rand.nextDouble()<=p)
				{
				spin[k]=0;
				deadsites++;
				}
		}
	}
	
	public int numberofbonds(int spin[], int totalspins)
	{
		int bondsnumber=0;
		for(int l=0;l<totalspins;l++)
			for(int m=0; m<4; m++)
			{
				bondsnumber+=spin[l]*spin[Nneighber(m,l)];
			}
		return bondsnumber;
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
	
	public static void main (String[] kangliu){
		new Control(new surfacetension(), "Kang Liu's surface tension for classical nucleation" );
	}
	
	public void load(Control surfacetension){
		surfacetension.frame (grid1);

		
		params.add("lattice's width", 100);
		params.add("lattice's length", 100);
		params.add("Diluted Percentage", new DoubleValue(0.30,0,1).withSlider());
		params.add("copies", 1000);
		params.add("runs");
		params.add("bonds");
		params.add("averagebonds");
		
		
	}
	
	public void run(){
	
		int u=0;

		numberofcopies=(int)params.fget("copies");


		L1 =(int)params.fget("lattice's width");
		L2 =(int)params.fget("lattice's length");
		M = L1 * L2;
		isingspin = new int[M];
		//percent = params.fget("Diluted Percentage");
		
		for(int v=0; v<1000; v++)
			{
			
			totalbonds=0;
			averagebonds=0;
			
			percent=v*0.001;
			for(u=0; u<numberofcopies; u++)
			{
			dilutionseed=u;
			dilutionrand= new Random(dilutionseed);
			initialize(isingspin, percent, dilutionrand,M);
			Job.animate();
			bonds=numberofbonds(isingspin,M);
			bondsratio=bonds/(4*M);
			totalbonds+=bonds;
			averagebonds=totalbonds/(u+1);
			params.set("bonds", bonds);
			params.set("averagebonds", averagebonds);
			params.set("runs", u);
			}
			PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/surfacetension30.txt", percent, averagebonds);
		}
		
		
		
	}
	
	

	
	
	
	
}