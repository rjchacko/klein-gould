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
	
	public String dynamics="Glauber";
	
	public IsingStructure3D IS;
	public IsingStructure3D Istemp;
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
	
	public int totalup, totaldown;
	public int absupdown;  //abs(totalup-totaldown)
	public double totalabsM;  //total absolute magnetization used to determine if the system is ordered
	
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

		params.add("L", 60);
		params.add("M");
	    params.add("deadsites");
	    params.add("livingsites");
		params.add("percent", 0.00);
		
		params.add("pb",0.358);     //bond probability
		params.add("pmin",0.35); 
		params.add("pmax",0.45); 
		params.add("increment",0.005); 
		
        params.add("Np");  //size of the largest cluster
        params.add("ratio",0.00);  //the raito of largest cluster/the total occupied sites
		
        params.add("run", 0);
        params.add("span",false);
		params.add("Dseed",1);    //seed for dilution configuration
		params.add("Pseed",1);    //seed for percolation
		
		params.add("prestep",20);
		params.add("steplimit",10000);
		params.add("MCS");
		params.add("T",4.525);
		
		params.add("magnetization", 0.0);
		params.add("totalup", 0.0);
		params.add("totaldown", 0.0);

		
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
		
		
		
		double totalspan=0;
		for(int r=0; r<totalruns; r++)
		{
			for(int j=0; j<(L*L*L); j++)
			{
				clustermap[j]=Ising.spin[j];
			}
			
			
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
	
	public void prepare(IsingStructure3D Ising,int prelimit, int steplimit, Random mcflip)
	{
		for(int prestep=0; prestep<prelimit; prestep++)
		{
			Ising.MCS(99, 0, mcflip, 1, "Glauber");
			for(int dis=0; dis<L*L; dis++)
			{
				isingdisplay[dis]=Ising.spin[L*L/2+dis];
			}
			Job.animate();
			params.set("MCS",prestep-prelimit);
		}
		for(int step=0; step<steplimit; step++)
		{
			T=params.fget("T");
			Ising.MCS(T, 0, mcflip, 1, "Glauber");
			for(int dis2=0; dis2<L*L; dis2++)
			{
				isingdisplay[dis2]=Ising.spin[L*L/2+dis2];
			}
			Job.animate();
			params.set("MCS",step);
		}
		if(Ising.Magnetization()>0)
			for(int jj=0; jj<(L*L*L);jj++)
			{
				Ising.spin[jj]=-Ising.spin[jj];  //here, choose the stable direction,but didn't change the energy and magnetization, be careful
			}
		
		for(int j=0; j<M; j++)
		{
			clustermap[j]=Ising.spin[j];
		}
		
		Job.animate();
	}
	
	public double susceptibility(IsingStructure3D Ising, double T, double H, int equstep, int length, int runs)//measure the susceptibility
	{
		double chi=0;
		
		double tempM[];
	    tempM= new double [length];
	    double totalX=0;
	    absupdown=0;
	    totalabsM=0;
		
		for(int rr=0; rr<runs; rr++)
		{
			params.set("run", rr+1);
			Istemp=Ising.clone();
			Random cflip=new Random(rr+1);
			totalup=0;
		    totaldown=0;
		    params.set("totalup", totalup);
			params.set("totaldown", totaldown);
			
			
			for(int heat=0; heat<50; heat++)
			{
				Istemp.MCS(9, H, cflip, 1, dynamics);
				params.set("magnetization", Istemp.Magnetization());
				for(int dis=0; dis<L*L; dis++)
				{
					isingdisplay[dis]=Istemp.spin[L*L/2+dis];
				}
				Job.animate();
				params.set("MCS", -9999);

			}
			for(int prestep=0; prestep< equstep; prestep++)
			{
				
				Istemp.MCS(T, H, cflip, 1, dynamics);
				params.set("magnetization", Istemp.Magnetization());
				if(prestep%20==0){
					for(int dis=0; dis<L*L; dis++)
					{
						isingdisplay[dis]=Istemp.spin[L*L/2+dis];
					}
					Job.animate();
				}
				params.set("MCS", prestep-equstep);
			}
			for(int step=0; step<length; step++)
			{
				
				//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/NNP3D/progress.txt", T, rr, step);
				Istemp.MCS(T, H, cflip, 1, dynamics);
				tempM[step]=Istemp.TotalSpin();
				totalabsM+=Math.abs(tempM[step]/Istemp.M);
				if(tempM[step]>0)
					totalup++;
				else
					totaldown++;
				params.set("magnetization", Istemp.Magnetization());
				params.set("totalup", totalup);
				params.set("totaldown", totaldown);
				Job.animate();
				if(step%20==0){
					for(int dis=0; dis<L*L; dis++)
					{
						isingdisplay[dis]=Istemp.spin[L*L/2+dis];
					}
					Job.animate();
				}
				
				params.set("MCS", step);
			}
			totalX+=IS.Fluctuation(tempM, length);
			absupdown+=Math.abs(totalup-totaldown);
		}
		
		chi=totalX/runs;
		return chi;
	}
	
	public void chiscanforTc(IsingStructure3D Ising, double minT, double maxT, double dT, int equstep, int length, int runs)//measure the susceptibility for different temperature to determine the critical temperature
	{
		double Xtemp;
		String path="/Users/liukang2002507/Desktop/simulation/NNP3D/chiforTc <L="+fmt.format(L)+"-q=0."+fmt.format(percent*1000)+"> Dseed="+fmt.format(Dseed)+".txt";
		for(double t=minT; t<maxT; t+=dT)
		{
			params.set("T", t);
			Xtemp=susceptibility(Ising, t, 0, equstep, length, runs);
			PrintUtil.printlnToFile(path, t, Xtemp, Xtemp/(Ising.M-Ising.deadsites), totalabsM, absupdown);
		}
		
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
		Istemp=new IsingStructure3D(L,L,L,0,-6,percent,percent,"square");
		IS.Dinitialization(Dseed, Dseed,10, 10, 10);
		IS.Sinitialization(0, Dseed);    //set all the occupied sites to have down spins(white in display)
		params.set("deadsites", IS.deadsites);
		params.set("livingsites", M-IS.deadsites);
		
		
		prepare(IS, prelimit, steplimit, mcflip);

		//singlerun(IS, Pseed, pb);
        //multiruns(IS, Pseed, pmin, pmax, increment);
		//averagerun(IS, 10, pb);
		
		multiaverageruns(IS, 10, pmin, pmax, increment);
		
		//chiscanforTc(IS, 3.39, 3.59, 0.005, 5000, 1000, 10);
		
	}
	
	
	
}


