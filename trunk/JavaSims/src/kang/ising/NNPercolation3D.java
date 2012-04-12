package kang.ising;

import java.awt.Color;
import java.text.DecimalFormat;

import kang.ising.BasicStructure.IsingStructure3D;
import kang.ising.BasicStructure.Cluster3D;
import kang.ising.BasicStructure.ClusterSet3D;
import kang.ising.BasicStructure.Percolation3D;

import chris.util.PrintUtil;
import chris.util.Random;

import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.params.DoubleValue;

public class NNPercolation3D extends Simulation
{
	Grid grid1=new Grid("grid1");     // the map to display the dilute ising configuration
	//Grid grid2=new Grid("grid2");     // the map to display the largest cluster
	Grid gridX=new Grid("X plane");
	Grid gridY=new Grid("Y plane");
	Grid gridZ=new Grid("Z plane");     // the map to display the largest cluster(choose 3 characteristic planes)
	
	
	
	public IsingStructure3D IS;
	public Percolation3D NNP;   //nearest neighbor percolation object
	
	private DecimalFormat fmt = new DecimalFormat("000");
	private DecimalFormat qmt = new DecimalFormat("000000");
	
	public int L;
	public double M;
	public double pb;
	public double pmin,pmax,increment;
	public double percent;
	public int deadsite;
	public double ratio;
	public int Dseed, Pseed;
	
	public double T;
	
	public int size;    // size of the largest cluster
	
	public int clustermap[];
	public int displayX[];
	public int displayY[];
	public int displayZ[];
	public int isingdisplay[];
	public int dx,dy,dz;   //x, y, z indices for side-show display
	
	
	public void animate()
	{
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);      //up spin
		ising.setColor(-1, Color.WHITE);     //down spin
		ising.setColor(0, Color.RED);        //normal dilution
		ising.setColor(2, Color.BLUE);       //clusters
		ising.setColor(-2, Color.GREEN);     //
		ising.setColor(3, Color.darkGray);    // the centers of the clusters
		
		
		grid1.setColors(ising);
		grid1.registerData(L, L, isingdisplay);
		gridX.setColors(ising);
		gridX.registerData(L, L, displayX);
		gridY.setColors(ising);
		gridY.registerData(L, L, displayY);
		gridZ.setColors(ising);
		gridZ.registerData(L, L, displayZ);
		
		
		
		//params.set("span", NNP.span);
		params.set("ratio", ratio);
	}
	
	public void clear()
	{
		grid1.clear();
		//grid2.clear();
		gridX.clear();
		gridY.clear();
		gridZ.clear();
	
	}
	
	public static void main (String[] NNPercolation3D){
		new Control(new NNPercolation3D(), "Kang Liu's nearest neighbor site-diluted ising model's 3D percolation problem" );
	}
	
	public void load(Control NNPercolation3D)
	{
		
		NNPercolation3D.frameTogether("Display", grid1, gridX, gridY, gridZ);

		params.add("L", 50);
		params.add("M");
	    params.add("deadsites");
	    params.add("livingsites");
		params.add("percent", 0.20);
		
		params.add("pb",0.358);     //bond probability
		params.add("pmin",0.35); 
		params.add("pmax",0.50); 
		params.add("increment",0.01); 
		
        params.add("Np");  //size of the largest cluster
        params.add("ratio",0.00);  //the raito of largest cluster/the total occupied sites
		
        params.add("run", 0);
        params.add("span",false);
		params.add("Dseed",1);    //seed for dilution configuration
		params.add("Pseed",1);    //seed for percolation
		
		params.add("prestep",20);
		params.add("steplimit",10000);
		params.add("MCS");
		params.add("T",3.50590);
		
	}
	
	public int ClusterSize(int map[], int L)
	{
		int size=0;
		for(int j=0; j<(L*L*L); j++)
		{
			if(map[j]==2)
				size++;
		}
		return size;
	}
	
	public void singlerun(IsingStructure3D Ising, int Pseed, double probability)
	{
		String data="/Users/liukang2002507/Desktop/simulation/NNP3D/data/L="+fmt.format(L)+"-q=0."+fmt.format(percent*1000)+"-T="+qmt.format(T*100000)+".txt";
		//String image="/Users/liukang2002507/Desktop/simulation/NNP/image/L="+fmt.format(L)+"-q=0."+fmt.format(percent*1000)+"-Pb=0."+qmt.format(probability*10000)+".txt";
		
		for(int j=0; j<(L*L*L); j++)
		{
			clustermap[j]=Ising.spin[j];
		}
		
		NNP=new Percolation3D(Ising, 1);// keep track of 2 largest clusters
		NNP.SetProbability(probability);
		NNP.fastNNMapping(Pseed);
		NNP.SpanningCluster();
		
		NNP.ClusterSize();
		NNP.SDClusterSize();
		NNP.OrderParameter();
		
		NNP.CS.set[NNP.CS.maximumpin].Center();
		dx=NNP.CS.set[NNP.CS.maximumpin].cx;
		dy=NNP.CS.set[NNP.CS.maximumpin].cy;
		dz=NNP.CS.set[NNP.CS.maximumpin].cz;

		Job.animate();
		params.set("span", NNP.span);
		for(int j=0; j<(L*L*L); j++)
		{
			if(Ising.spin[j]==-1)
			{
				if(NNP.CS.set[NNP.CS.maximumpin].lattice[j]==2)
					{
					clustermap[j]=2;
					}
			}
		}
		for(int di=0; di<L; di++)
			for(int dj=0; dj<L; dj++)
			{
				displayX[di*L+dj]=clustermap[dx*L*L+di*L+dj];
				displayY[di*L+dj]=clustermap[dj*L*L+dy*L+di];
				displayZ[di*L+dj]=clustermap[di*L*L+dj*L+dz];	
			}
		
		//size=ClusterSize(clustermap,L);
		//ratio=size/(M-deadsite);
		Job.animate();
		if(NNP.span)
			{
			PrintUtil.printlnToFile(data, probability, NNP.OP, NNP.meanclustersize, NNP.SDclustersize, NNP.totalclusters, 1, NNP.spannumber);
			}
		else
			{
			PrintUtil.printlnToFile(data, probability, NNP.OP, NNP.meanclustersize, NNP.SDclustersize, NNP.totalclusters, 0, NNP.spannumber);
			}
		//save image
	}
	
	public void multiruns(IsingStructure3D Ising, int Pseed, double pmin, double pmax, double increment)
	{
		int ii=1;
		for(double pp=pmin; pp<=pmax; pp+=increment)
		{
			ii++;
			
			params.set("pb", pp);
			singlerun(Ising, ii, pp);
		}
	}
	
	public void multiaverageruns(IsingStructure3D Ising, int totalruns, double pmin, double pmax, double increment)
	{
		int ii=1;
		for(double pp=pmin; pp<=pmax; pp+=increment)
		{
			ii++;
			
			params.set("pb", pp);
			averagerun(Ising, totalruns, pp);
		}
	}
	
	public void averagerun(IsingStructure3D Ising, int totalruns, double probability)
	{
		String averagedata="/Users/liukang2002507/Desktop/simulation/NNP3D/averagedata/L="+fmt.format(L)+"-q=0."+fmt.format(percent*1000)+"-T="+qmt.format(T*100000)+".txt";
		//String image="/Users/liukang2002507/Desktop/simulation/NNP/image/L="+fmt.format(L)+"-q=0."+fmt.format(percent*1000)+"-Pb=0."+qmt.format(probability*10000)+".txt";
		
		for(int j=0; j<(L*L*L); j++)
		{
			clustermap[j]=Ising.spin[j];
		}
		
		double totalspan=0;
		for(int r=0; r<totalruns; r++)
		{
			params.set("run", r+1);
			NNP=new Percolation3D(Ising, 1);// keep track of the largest cluster
			NNP.SetProbability(probability);
			NNP.fastNNMapping(r+1);
			NNP.SpanningCluster();
			
			NNP.ClusterSize();
			NNP.SDClusterSize();
			NNP.OrderParameter();
			
			NNP.CS.set[NNP.CS.maximumpin].Center();
			dx=NNP.CS.set[NNP.CS.maximumpin].cx;
			dy=NNP.CS.set[NNP.CS.maximumpin].cy;
			dz=NNP.CS.set[NNP.CS.maximumpin].cz;

			Job.animate();
			params.set("span", NNP.span);
			if(NNP.span)
			{
				totalspan++;
			}
			
			
			for(int j=0; j<(L*L*L); j++)
			{
				if(Ising.spin[j]==-1)
				{
					if(NNP.CS.set[NNP.CS.maximumpin].lattice[j]==2)
						{
						clustermap[j]=2;
						}
				}
			}
			for(int di=0; di<L; di++)
				for(int dj=0; dj<L; dj++)
				{
					displayX[di*L+dj]=clustermap[dx*L*L+di*L+dj];
					displayY[di*L+dj]=clustermap[dj*L*L+dy*L+di];
					displayZ[di*L+dj]=clustermap[di*L*L+dj*L+dz];	
				}
			
			//size=ClusterSize(clustermap,L);
			//ratio=size/(M-deadsite);
			Job.animate();
		}
		

		PrintUtil.printlnToFile(averagedata, probability, totalspan/totalruns);

	}
	
	public void run()
	{
		
		
		L = (int)params.fget("L");
		M = L*L*L;
		
		size=0;
		ratio=0;
		percent = params.fget("percent");
		Dseed = (int)params.fget("Dseed");
		Pseed = (int)params.fget("Pseed");
		clustermap= new int[L*L*L];
		displayX=new int[L*L];
		displayY=new int[L*L];
		displayZ=new int[L*L];
		isingdisplay=new int[L*L];
		
		pb= params.fget("pb");
		pmin= params.fget("pmin");
		pmax= params.fget("pmax");
		increment= params.fget("increment");
		int prelimit=(int)params.fget("prestep");
		int steplimit=(int)params.fget("steplimit");;
		Random mcflip=new Random(47);
		
		IS=new IsingStructure3D(L,L,L,0,-6,percent,percent,"square");
		IS.Dinitialization(Dseed, Dseed,10, 10, 10);
		IS.Sinitialization(0, Dseed);    //set all the occupied sites to have down spins(white in display)
		params.set("deadsites", IS.deadsites);
		params.set("livingsites", M-IS.deadsites);
		for(int prestep=0; prestep<prelimit; prestep++)
		{
			IS.MCS(99, 0, mcflip, 1, "Glauber");
			for(int dis=0; dis<L*L; dis++)
			{
				isingdisplay[dis]=IS.spin[L*L/2+dis];
			}
			Job.animate();
			params.set("MCS",prestep-prelimit);
		}
		for(int step=0; step<steplimit; step++)
		{
			T=params.fget("T");
			IS.MCS(T, 0, mcflip, 1, "Glauber");
			for(int dis2=0; dis2<L*L; dis2++)
			{
				isingdisplay[dis2]=IS.spin[L*L/2+dis2];
			}
			Job.animate();
			params.set("MCS",step);
		}
		if(IS.Magnetization()>0)
			for(int jj=0; jj<(L*L*L);jj++)
			{
				IS.spin[jj]=-IS.spin[jj];  //here, choose the stable direction,but didn't change the energy and magnetization, be careful
			}
		
		for(int j=0; j<M; j++)
		{
			clustermap[j]=IS.spin[j];
		}
		
		Job.animate();

		//singlerun(IS, Pseed, pb);
        //multiruns(IS, Pseed, pmin, pmax, increment);
		//averagerun(IS, 10, pb);
		multiaverageruns(IS, 10, pmin, pmax, increment);
		
	}
	
	
	
}


