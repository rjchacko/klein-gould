package kang.AEM;

import java.awt.Color;
import java.text.DecimalFormat;


import kang.AEM.BasicStructure.AEMStructure;
import kang.util.PrintUtil;
import kang.ising.BasicStructure.BasicTools;
import kang.ising.BasicStructure.IsingStructure;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import chris.util.Random;


public class Globaltrading extends Simulation{
	
	Grid grid1=new Grid("wealth distribution");     // the map to display the wealth distribution	
	
	
	public AEMStructure AH;
	public AEMStructure AHtemp;


	
	public BasicTools Tools;
	
	private DecimalFormat fmt = new DecimalFormat("000");
	private DecimalFormat bmt = new DecimalFormat("0000");
	
	
	
	//initialization parameters
	public int L,R;
	public double M;


	public double percent;
	public double tax;
	public double growth;

	public int time;
	public double order;
	

	
	public void animate()
	{
	
	    ColorGradient heatmap = new ColorGradient();
		
		
		grid1.setColors(heatmap);
		grid1.registerData(L, L, AHtemp.wealth);
	
		
	}
	
	public void clear()
	{
		grid1.clear();
	
	}
	
	public static void main (String[] Globaltrading){
		new Control(new Globaltrading(), "Kang Liu's global trading AEM" );
	}
	
	public void printdata(String path, int[] data)
	{
		int j=0;
		while(data[j]>=0)
		{
			PrintUtil.printlnToFile(path, j+1, data[j]);
			j++;
		}
	}
	
	
	public void load(Control Globaltrading)
	{

		Globaltrading.frameTogether("Display", grid1);

		params.add("L", 800);
		params.add("R", 800);
		
	    params.add("tax");
		params.add("percent", 0.00);
		params.add("growth", 0.00);
		
		params.addm("time", 0);
		params.addm("order", 1.26);

	}
	
	public void run(){
		
		
		L = (int)params.fget("L");
		R =(int)params.fget("R");
		M = L * L;


		percent=params.fget("percent");
		tax=params.fget("tax");
		growth=params.fget("growth");
		
		
	    AH=new AEMStructure(L,L,R,percent,tax, 0, growth);   
	    AHtemp=new AEMStructure(L,L,R,percent,tax, 0,growth);

	    
	    Tools=new BasicTools();

	    
	    {//initialization
	    	
	    	AH.Uinitialization(100);

	    	AHtemp=AH.clone();
	    	
	    }
	    
	    Job.animate();
	   
	    singlerun(AHtemp, steplimit, 2000, true, 1);

	    
	    Job.animate();

	}
	
}