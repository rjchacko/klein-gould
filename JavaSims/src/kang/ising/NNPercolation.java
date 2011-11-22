package kang.ising;

import java.awt.Color;
import java.text.DecimalFormat;

import kang.ising.BasicStructure.IsingStructure;
import kang.ising.BasicStructure.Cluster;
import kang.ising.BasicStructure.ClusterSet;
import kang.ising.BasicStructure.Percolation;

import chris.util.PrintUtil;
import chris.util.Random;

import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.params.DoubleValue;

public class NNPercolation extends Simulation
{
	//Grid grid1=new Grid("grid1");     // the map to display the dilution configuration
	Grid grid2=new Grid("grid2");     // the map to display the largest cluster
	
	public IsingStructure IS;
	public Percolation NNP;   //nearest neighbor percolation object
	
	private DecimalFormat fmt = new DecimalFormat("000");
	private DecimalFormat qmt = new DecimalFormat("00000");
	
	public int L;
	public double M;
	public double pb;
	public double pmin,pmax,increment;
	public double percent;
	public int deadsite;
	public double ratio;
	public int Dseed, Pseed;
	public boolean span;
	
	public int size;    // size of the largest cluster
	
	public int clustermap[];
	
	
	public void animate()
	{
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);      //up spin
		ising.setColor(-1, Color.WHITE);     //down spin
		ising.setColor(0, Color.RED);        //normal dilution
		ising.setColor(2, Color.BLUE);       //clusters
		ising.setColor(-2, Color.GREEN);     //
		ising.setColor(3, Color.darkGray);    // the centers of the clusters
		
		
		//grid1.setColors(ising);
		//grid1.registerData(IS.L1, IS.L2, IS.spin);
		grid2.setColors(ising);
		grid2.registerData(IS.L1, IS.L2, clustermap);
		params.set("span", span);
		params.set("ratio", ratio);
	}
	
	public void clear()
	{
		//grid1.clear();
		grid2.clear();
	}
	
	public static void main (String[] NNPercolation){
		new Control(new NNPercolation(), "Kang Liu's nearest neighbor site-diluted ising model's percolation problem" );
	}
	
	public void load(Control NNPercolation)
	{
		//NNPercolation.frame (grid1);
		NNPercolation.frame (grid2);

		params.add("L", 100);
		params.add("M");
	    params.add("deadsites");
	    params.add("livingsites");
		params.add("percent", 0.0);
		
		params.add("pb",1.0);     //bond probability
		params.add("pmin",0.01); 
		params.add("pmax",1.0); 
		params.add("increment",0.01); 
		
        params.add("Np");  //size of the largest cluster
        params.add("ratio",0.00);  //the raito of largest cluster/the total occupied sites
		
        params.add("span",false);
		params.add("Dseed",1);    //seed for dilution configuration
		params.add("Pseed",1);    //seed for percolation
		
		params.add("check",99.0);
		

	}
	
	public boolean SpanningCluster(int map[], int L)
	{
		//the logic of determine if the cluster spans the whole lattice:
		//1 the coordinate x covers 0 to L-1
		//2 the coordinate y covers 0 to L-1
		//3 the site of cluster on edge has neighbor on the other edge
		boolean span=false;
		boolean spanX=true;
		boolean spanY=true;
		
		
		int X[]=new int[L];
		int Y[]=new int[L];
		int x, y;
		//int x1, y1, xL, yL;   //the coordinate j of the site in cluster on the edge
		
		for(int i=0; i<L; i++)
		{
			X[i]=1;
			Y[i]=1;
			
		}// preset the array to count the x and y coordinates
		for(int j=0; j<(L*L); j++)
		{
			if(map[j]==2)
			{
				x=j/L;
				y=j%L;
				X[x]=0;
				Y[y]=0;
				
			}
		}// scan the whole cluster distribution to mark all the occupied x and y coordinates 
		for(int k=0; k<L; k++)
		{
			if(X[k]!=0)
				spanX=false;
			if(Y[k]!=0)
				spanY=false;
		}

		if(spanX)
			span=true;
		if(spanY)
			span=true;
		 
		return span;
	}
	
	public int ClusterSize(int map[], int L)
	{
		int size=0;
		for(int j=0; j<(L*L); j++)
		{
			if(map[j]==2)
				size++;
		}
		return size;
	}
	
	public void singlerun(IsingStructure Ising, int Pseed, double probability)
	{
		String data="/Users/liukang2002507/Desktop/simulation/NNP/data/L="+fmt.format(L)+"-q=0."+fmt.format(percent*1000)+".txt";
		//String image="/Users/liukang2002507/Desktop/simulation/NNP/image/L="+fmt.format(L)+"-q=0."+fmt.format(percent*1000)+"-Pb=0."+qmt.format(probability*10000)+".txt";
		
		for(int j=0; j<(L*L); j++)
		{
			clustermap[j]=Ising.spin[j];
		}
		
		Job.animate();
		
		NNP=new Percolation(Ising, 1);// keep track of 2 largest clusters
		NNP.SetProbability(probability);
		NNP.fastNNMapping(Pseed);
		
		for(int j=0; j<(L*L); j++)
		{
			if(Ising.spin[j]==-1)
			{
				if(NNP.CS.set[NNP.CS.maximumpin].lattice[j]==2)
					{
					clustermap[j]=2;
					}
			}
		}
		
		span=SpanningCluster(clustermap,L);
		size=ClusterSize(clustermap,L);
		ratio=size/(M-deadsite);
		Job.animate();
		if(span)
			PrintUtil.printlnToFile(data, probability, ratio, 1);
		else
			PrintUtil.printlnToFile(data, probability, ratio, 0);
		//save image
	}
	
	public void multiruns(IsingStructure Ising, int Pseed, double pmin, double pmax, double increment)
	{
		int ii=1;
		for(double pp=pmin; pp<=pmax; pp+=increment)
		{
			ii++;
			
			params.set("pb", pp);
			singlerun(Ising, ii, pp);
		}
	}
	
	public void run()
	{
		
		
		L = (int)params.fget("L");
		M = L*L;
		
		size=0;
		ratio=0;
		percent = params.fget("percent");
		Dseed = (int)params.fget("Dseed");
		Pseed = (int)params.fget("Pseed");
		clustermap= new int[L*L];
		pb= params.fget("pb");
		pmin= params.fget("pmin");
		pmax= params.fget("pmax");
		increment= params.fget("increment");
		
		IS=new IsingStructure(L,L,0,-4,percent,percent,"square");
		IS.Dinitialization(Dseed, Dseed, 10, 10);
		IS.Sinitialization(2, Dseed);    //set all the occupied sites to have down spins(white in display)
		params.set("deadsites", IS.deadsites);
		params.set("livingsites", M-IS.deadsites);
		
		for(int j=0; j<M; j++)
		{
			clustermap[j]=IS.spin[j];
		}
		
		Job.animate();

		//singlerun(IS, Pseed, pb);
        multiruns(IS, Pseed, pmin, pmax, increment);
		
	}
	
	
	
}


