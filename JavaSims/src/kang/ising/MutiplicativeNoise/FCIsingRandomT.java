package kang.ising.MutiplicativeNoise;

import java.text.DecimalFormat;

import kang.ising.BasicStructure.FCIsing;
import kang.ising.BasicStructure.MeanFieldIsingStructure;

import java.awt.Color;
import java.text.DecimalFormat;

import chris.util.PrintUtil;
import chris.util.Random;

import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.Control;
import scikit.jobs.params.DoubleValue;


public class FCIsingRandomT extends Simulation
{
	public double NJ;
	public double T,H;
	public FCIsing ising;
	
	private DecimalFormat fmt = new DecimalFormat("000");
	
 	public void Heatup(FCIsing ising, int steplimit)
 	{
	    Random heatrand=new Random(999);
 	    for(int heat=0; heat<steplimit; heat++)
	    {
	    	ising.MCS(heatrand,heatrand,999,0,1);
	    	Job.animate();
			params.set("MCS", heat-steplimit);
	    }
 	}
 	
 	
 	
	
}