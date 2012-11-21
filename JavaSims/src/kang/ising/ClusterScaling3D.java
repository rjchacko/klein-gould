package kang.ising;

import java.awt.Color;
import java.text.DecimalFormat;

import kang.ising.BasicStructure.BasicTools;
import kang.ising.BasicStructure.IsingStructure3D;
import kang.util.PrintUtil;
import kang.ising.BasicStructure.BasicTools;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import chris.util.Random;


public class ClusterScaling3D extends Simulation{
	
	Grid grid1=new Grid("x ising");     // the map to display the dilution configuration
	Grid grid2=new Grid("y ising");     // the map to display the simulation
	Grid grid3=new Grid("z ising");
	Grid gridx=new Grid("x cluster");
	Grid gridy=new Grid("y cluster");
	Grid gridz=new Grid("z cluster");
	
	public IsingStructure3D IS;
	public IsingStructure3D Istemp;

	public Random Erand;
	
	public BasicTools Tools;
	
	private DecimalFormat fmt = new DecimalFormat("000");
	private DecimalFormat bmt = new DecimalFormat("0000");
	
	
	
	//initialization parameters
	public int L,la,lb,lc,R, Lp;
	public double M;
	public double NJ;
	public int Dseed, Bseed, Sseed;
	public double percent;
	public double biaspercent;
	public int deadsite;
	public String dynamics;

	//dynamic parameters
	public double T, H;
	
	//cluster parameters
	public int[] clustersize;
	
	//arrays to display in 3D
	public int isingx[];
	public int isingy[];
	public int isingz[];
	
	public int clusterx[];
	public int clustery[];
	public int clusterz[];
	
	public int clcx, clcy, clcz;  // location of the cluster center in 3D
	
	
	public void animate()
	{
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);      //up spin
		ising.setColor(-1, Color.WHITE);     //down spin
		ising.setColor(0, Color.RED);        //normal dilution
		ising.setColor(2, Color.BLUE);       //clusters
		ising.setColor(-2, Color.GREEN);     //
		ising.setColor(3, Color.darkGray);    // the centers of the clusters

		for(int i=0;i<L; i++)
			for(int j=0; j<L; j++)
			{
				isingx[i*L+j]=Istemp.spin[L/2*L*L+i*L+j];
				isingy[i*L+j]=Istemp.spin[i*L*L+L/2*L+j];
				isingz[i*L+j]=Istemp.spin[i*L*L+j*L+L/2];
			}
		
		grid1.setColors(ising);
		grid1.registerData(L, L, isingx);
		grid2.setColors(ising);
		grid2.registerData(L, L, isingy);
		grid3.setColors(ising);
		grid3.registerData(L, L, isingz);
		
		for(int i=0;i<L; i++)
			for(int j=0; j<L; j++)
			{
				clusterx[i*L+j]=Istemp.largestcluster[clcx*L*L+i*L+j];
				clustery[i*L+j]=Istemp.largestcluster[i*L*L+clcy*L+j];
				clusterz[i*L+j]=Istemp.largestcluster[i*L*L+j*L+clcz];
			}
		
		
		gridx.setColors(ising);
		gridx.registerData(L, L, clusterx);
		gridy.setColors(ising);
		gridy.registerData(L, L, clustery);
		gridz.setColors(ising);
		gridz.registerData(L, L, clusterz);
		
	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
		grid3.clear();

		gridx.clear();
		gridy.clear();
		gridz.clear();
	}
	
	public static void main (String[] ClusterScaling3D){
		new Control(new ClusterScaling3D(), "Kang Liu's cluster scaling in 3D" );
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
	
	public int adddata(String path, int[] data, int start)
	{
		int j=0;
		int end=0;
		while(data[j]>=0)
		{
			PrintUtil.printlnToFile(path, j+start, data[j]);
			j++;
		}
		
		end=j+start;
		return end;
	}
	
	public void singlerunTc(IsingStructure3D ising, double T, int steplimit, boolean keeplargest,int seed)
	{
		String singlerun="<L="+fmt.format(L)+", R="+fmt.format(R)+", T="+fmt.format(T*100)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String singlepath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling3D/"+dynamics+"/"+singlerun+"positive.txt";
		String secondpath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling3D/"+dynamics+"/"+singlerun+"negative.txt";
		
		Random rand= new Random(seed);
		for(int prestep=0; prestep<50; prestep++)
		{
			params.set("T", 999);
			params.set("H", 0);
			ising.MCS(999, 0, rand, 1, dynamics);
			Job.animate();
			params.set("Emcs", prestep-50);
			params.set("magnetization", ising.magnetization);
		}
		
		for(int step=0; step<steplimit; step++)
		{
			params.set("T", T);
			params.set("H", 0);
			ising.MCS(T, 0, rand, 1, dynamics);
			
			Job.animate();
			params.set("Emcs", step);
			params.set("magnetization", ising.magnetization);
			
		}
		double pb=1-Math.exp(2*ising.J/T);
		int direction=1;
		if(ising.magnetization<0)
			direction=-1;
			
		clustersize=ising.Clustergrowth(ising.spin, direction, pb, seed, seed, keeplargest);
		int center=ising.ClusterInfo(ising.largestcluster)[1];
		clcx=ising.getX(center);
		clcy=ising.getY(center);
		clcz=ising.getZ(center);
		
		Job.animate();
		printdata(singlepath, clustersize);
		
		clustersize=ising.Clustergrowth(ising.spin, -direction, pb, seed, seed, keeplargest);
	    center=ising.ClusterInfo(ising.largestcluster)[1];
		clcx=ising.getX(center);
		clcy=ising.getY(center);
		clcz=ising.getZ(center);
		
		Job.animate();
		printdata(secondpath, clustersize);
	}
	
	public void multiplerunsTc(IsingStructure3D ising, double T, int steplimit, boolean keeplargest,int seed, int copies)// average over multiple spin configurations, positve represent all up spin clusters
	{
		String multirun="multi <L="+fmt.format(L)+", R="+fmt.format(R)+", T="+fmt.format(T*100)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String positivepath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling3D/"+dynamics+"/"+multirun+"positive.txt";
		String negativepath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling3D/"+dynamics+"/"+multirun+"negative.txt";
		
		int spincopies[][]=new int [copies][ising.M];
		
		
		Random rand= new Random(seed);
		for(int prestep=0; prestep<50; prestep++)
		{
			params.set("T", 999);
			params.set("H", 0);
			ising.MCS(999, 0, rand, 1, dynamics);
			Job.animate();
			params.set("Emcs", prestep-50);
			params.set("magnetization", ising.magnetization);
		}
		
		for(int step=0; step<steplimit; step++)
		{
			params.set("T", T);
			params.set("H", 0);
			ising.MCS(T, 0, rand, 1, dynamics);
			
			Job.animate();
			params.set("Emcs", step);
			params.set("magnetization", ising.magnetization);
			
		}
		
		double pb=1-Math.exp(2*ising.J/T);
		
		for(int astep=0; astep<(10*copies); astep++)
		{
			params.set("T", T);
			params.set("H", 0);
			ising.MCS(T, 0, rand, 1, dynamics);
			
			Job.animate();
			params.set("Emcs", steplimit+astep);
			params.set("magnetization", ising.magnetization);
			
			if(astep%10==0)
			{
				for(int jj=0; jj<ising.M; jj++)
				{
					spincopies[astep/10][jj]=ising.spin[jj];
				}
			}
		}
		
		int Pstart=0;
		int Nstart=0;
        for(int cc=0; cc<copies; cc++)
        {
        	
        	clustersize=ising.Clustergrowth(spincopies[cc], 1, pb, seed, seed, keeplargest);
    		int center=ising.ClusterInfo(ising.largestcluster)[1];
    		clcx=ising.getX(center);
    		clcy=ising.getY(center);
    		clcz=ising.getZ(center);
    		
    		Job.animate();
    		Pstart=adddata(positivepath, clustersize, Pstart);
    		
    		clustersize=ising.Clustergrowth(spincopies[cc], -1, pb, seed, seed, keeplargest);
    	    center=ising.ClusterInfo(ising.largestcluster)[1];
    		clcx=ising.getX(center);
    		clcy=ising.getY(center);
    		clcz=ising.getZ(center);
    		
    		Job.animate();
    		Nstart=adddata(negativepath, clustersize, Nstart);
        	
        }
	
			
		
	}
	
	public void wrongrunsTc(IsingStructure3D ising, double T, int steplimit, boolean keeplargest,int seed, int copies)// the runs with pb=1, just to compare with the right bond probability
	{
		String wrongrun="wrong <L="+fmt.format(L)+", R="+fmt.format(R)+", T="+fmt.format(T*100)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String positivepath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling3D/"+dynamics+"/"+wrongrun+"positive.txt";
		String negativepath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling3D/"+dynamics+"/"+wrongrun+"negative.txt";
		
		int spincopies[][]=new int [copies][ising.M];
		
		
		Random rand= new Random(seed);
		for(int prestep=0; prestep<50; prestep++)
		{
			params.set("T", 999);
			params.set("H", 0);
			ising.MCS(999, 0, rand, 1, dynamics);
			Job.animate();
			params.set("Emcs", prestep-50);
			params.set("magnetization", ising.magnetization);
		}
		
		for(int step=0; step<steplimit; step++)
		{
			params.set("T", T);
			params.set("H", 0);
			ising.MCS(T, 0, rand, 1, dynamics);
			
			Job.animate();
			params.set("Emcs", step);
			params.set("magnetization", ising.magnetization);
			
		}
		
		double pb=1;
		
		for(int astep=0; astep<(10*copies); astep++)
		{
			params.set("T", T);
			params.set("H", 0);
			ising.MCS(T, 0, rand, 1, dynamics);
			
			Job.animate();
			params.set("Emcs", steplimit+astep);
			params.set("magnetization", ising.magnetization);
			
			if(astep%10==0)
			{
				for(int jj=0; jj<ising.M; jj++)
				{
					spincopies[astep/10][jj]=ising.spin[jj];
				}
			}
		}
		
		int Pstart=0;
		int Nstart=0;
        for(int cc=0; cc<copies; cc++)
        {
        	
        	clustersize=ising.Clustergrowth(spincopies[cc], 1, pb, seed, seed, keeplargest);
    		int center=ising.ClusterInfo(ising.largestcluster)[1];
    		clcx=ising.getX(center);
    		clcy=ising.getY(center);
    		clcz=ising.getZ(center);
    		
    		Job.animate();
    		Pstart=adddata(positivepath, clustersize, Pstart);
    		
    		clustersize=ising.Clustergrowth(spincopies[cc], -1, pb, seed, seed, keeplargest);
    	    center=ising.ClusterInfo(ising.largestcluster)[1];
    		clcx=ising.getX(center);
    		clcy=ising.getY(center);
    		clcz=ising.getZ(center);
    		
    		Job.animate();
    		Nstart=adddata(negativepath, clustersize, Nstart);
        	
        }
	
			
		
	}
	
	public void singlerunHs(IsingStructure3D ising, double T, double H, int steplimit, boolean keeplargest,int seed)
	{
		String singlerun="Hs <L="+fmt.format(L)+", R="+fmt.format(R)+", T="+fmt.format(T*100)+", H="+bmt.format(H*1000)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String singlepath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling3D/"+dynamics+"/"+singlerun+".txt";
		
		Random rand= new Random(seed);
        double totalm=0;
        int mnumber=0;
		
		for(int step=0; step<steplimit; step++)
		{
			params.set("T", T);
			params.set("H", -H);
			ising.MCS(T, -H, rand, 1, dynamics);
			
			Job.animate();
			params.set("Emcs", step);
			params.set("magnetization", ising.magnetization);
			if(step>steplimit-100)
			{
				totalm+=ising.magnetization;
				mnumber++;
			}
			
			
		}
		
		double mag=totalm/mnumber;
	    double pb=1-Math.exp((1+mag/(1-ising.percent))*ising.J/T);
		int direction=-1;
			
		clustersize=ising.Clustergrowth(ising.spin, direction, pb, seed, seed, keeplargest);
		printdata(singlepath, clustersize);
		
	}
	
	public void approachHs(IsingStructure3D ising, double T, double Hmin, double Hmax, double dH, int steplimit, int seed)
	{
		String SPCrun="SPC data <L="+fmt.format(L)+", R="+fmt.format(R)+", T="+fmt.format(T*100)+"p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String SPCpath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling3D/"+dynamics+"/"+SPCrun+".txt";
		double pb=1-Math.exp((1+0.745/(1-ising.percent))*ising.J/T);;     //need to use the bond probability at spinodal
		
		int spantemp=0;
		for(double h=Hmin; h<Hmax; h+=dH)
		{
		  
			spantemp=findlargestcluster(ising, T, h, steplimit, pb, seed);
			PrintUtil.printlnToFile(SPCpath, h, spantemp);
			
		}
		
	}
	
	public int findlargestcluster(IsingStructure3D ising, double T, double H, int steplimit, double pb, int seed)
	{
		
		int largestsize=0;
		
		String FLCrun="FLC <L="+fmt.format(L)+", R="+fmt.format(R)+", T="+fmt.format(T*100)+", H="+bmt.format(H*1000)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String FLCpath = "/Users/liukang2002507/Desktop/simulation/ClusterScaling3D/"+dynamics+"/"+FLCrun+".txt";
		String FLCpic = "/Users/liukang2002507/Desktop/simulation/ClusterScaling3D/"+dynamics+"/FLCpic/cluster";
		
		ising.Sinitialization(1, Sseed);
		Random rand= new Random(seed);
		for(int step=0; step<steplimit; step++)
		{
			params.set("T", T);
			params.set("H", -H);
			ising.MCS(T, -H, rand, 1, dynamics);
			
			Job.animate();
			params.set("Emcs", step);
			params.set("magnetization", ising.magnetization);
		}
		clustersize=ising.Clustergrowth(ising.spin, -1, pb, seed, seed, true);
		Job.animate();
		printdata(FLCpath, clustersize);
		Tools.Picture(grid3, (int)(H*1000), (int)(percent*1000), FLCpic);
		
		int[] clinfotemp=new int [2];
		clinfotemp=ising.ClusterInfo(ising.largestcluster);
		largestsize=clinfotemp[0];
		
		return largestsize;
	}
	
	public void load(Control ClusterScaling)
	{

		ClusterScaling.frameTogether("Display", grid1 ,grid2, grid3, gridx, gridy, gridz);

		params.add("L", 100);
		
		params.add("la",10);    // scale of the bias dilution region
		params.add("lb",10); 
		params.add("lc",10); 
		params.add("R", 0);
		
		params.add("NJ", -6.0);
	    params.add("deadsites");

		params.add("percent", 0.10);
		params.add("biaspercent", 0.10);
		
		 		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));

	
		params.addm("T", 8.00);
		params.addm("H", 0.0);
		params.add("Emcs");    //MCS time for evolution
	
		    
		params.add("magnetization");

	}
	
	
	public void run(){
		
		
		L = (int)params.fget("L");

		la = (int)params.fget("la");
		lb = (int)params.fget("lb");
		lc = (int)params.fget("lc");
		R =(int)params.fget("R");
		M = L*L*L;
		NJ = params.fget("NJ");
		
		isingx=new int[L*L];
		isingy=new int[L*L];
		isingz=new int[L*L];
		
		clusterx=new int[L*L];
		clustery=new int[L*L];
		clusterz=new int[L*L];

		percent=params.fget("percent");
		biaspercent=params.fget("biaspercent");
		dynamics= params.sget("Dynamics");
		
		Dseed = 1;
		Bseed = 1;
		Sseed = 1;

		
	    IS=new IsingStructure3D(L,L,L,R,NJ,percent,biaspercent,"square");   
	    Istemp=new IsingStructure3D(L,L,L,R,NJ,percent,biaspercent,"square");

	    
	    Tools=new BasicTools();
	    T=params.fget("T");
	    H=params.fget("H");
	    
	    {//initialization
	    	
	    	IS.Dinitialization(Dseed, Bseed, la, lb, lc);
	    	params.set("deadsites",IS.deadsites);
	    	IS.Sinitialization(1, Sseed);
	    	Istemp=IS.clone();
	    	Istemp.largestcluster=new int[IS.M];
	    	
	    }
	    
	    Job.animate();
	   
	    singlerunTc(Istemp, T, 2000, true, 1);
	    
	    //multiplerunsTc(Istemp, T, 2000, true, 1, 10);
	    //singlerunHs(Istemp, T, H, 200, true, 1);
	    
	    //approachHs(Istemp, T, 1.26, 1.265, 0.001, 100, 1);
	    
	    Job.animate();

	}
	
}
