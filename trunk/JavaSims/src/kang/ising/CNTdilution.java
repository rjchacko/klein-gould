package kang.ising;

import java.awt.Color;

import java.text.DecimalFormat;



import kang.util.PrintUtil;
import chris.util.Random;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;



import kang.ising.BasicStructure.IsingStructure;
import kang.ising.BasicStructure.BasicTools;
import kang.ising.BasicStructure.Percolation;



public class CNTdilution extends Simulation{
	
	Grid grid1=new Grid("dilution configuration");     // the map to display the dilution configuration
	Grid grid2=new Grid("simulation");     // the map to display the simulation
	Grid grid3=new Grid("intervention copy");
	Grid grid4=new Grid("evolution display"); 
	Grid grid5=new Grid("total spin");
	Grid grid6=new Grid("dilution map");
	
	public IsingStructure IS;
	public IsingStructure Istemp;
	public IsingStructure Intervention;
	public IsingStructure Evolution;
	public Random Erand;
	
	public BasicTools Tools;
	
	private DecimalFormat fmt = new DecimalFormat("000");
	private DecimalFormat bmt = new DecimalFormat("0000");
	
	
	
	//initialization parameters
	public int L,la,lb;
	public double M;
	public double NJ;
	public int Dseed, Bseed, Sseed;
	public double percent;
	public double biaspercent;
	public int deadsite;
	public String dynamics;
	
	
	//dynamic parameters
	public double T, H;
	
	//intervention parameters
	public int totalruns;
	public int grownumber, decaynumber;
	public double threshold;
	public int breakpoint;
	public int steplimit;
	
	//percolation mapping parameters for determining the droplet size with pb=1
	public Percolation Droplet;
	public int dropletsize;
	
	//droplet distribution parameter
	public double spintotal[];
	public double dilutionmap[];
    public double meanNT;    //mean nucleation time
    public double SDNT;   //standard deviation of nucleation time
    public double nucleationevents[];
    public double lifetime[];
	
    public int surfacedilution[];
    public int surfacemetaspin[];
    public int surfacemap[];  //the map to record the configuration of surface of the droplet: 0-background 1-surface
	
	//public int ti, tf, tm;     //the MCS time of the intervention range[ti,tf] tm=(ti+tf)/2
	//public Random randi;  //the random number @ the beginning of the intervention range
	//public IsingStructure isingi;   //the ising configuration @ the beginning of the intervention range
	
	
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
		grid1.registerData(L, L, IS.spin);
		grid2.setColors(ising);
		grid2.registerData(L, L, Istemp.spin);
		grid3.setColors(ising);
		grid3.registerData(L, L, Intervention.spin);
		grid4.setColors(ising);
		grid4.registerData(L, L, Evolution.spin);
		
		
		ColorGradient heatmap = new ColorGradient();
		grid5.setColors(heatmap);
		grid5.registerData(L, L, spintotal);
		grid6.setColors(heatmap);
		grid6.registerData(L, L, dilutionmap);
	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
		grid3.clear();
		grid4.clear();
		grid5.clear();
		grid6.clear();
	}

	public static void main (String[] CNTdilution){
		new Control(new CNTdilution(), "Kang Liu's nearest neighbor site-diluted ising model's nucleation" );
	}
	
	
	
	public void load(Control CNTdilution)
	{

		CNTdilution.frameTogether("Display", grid1 ,grid2, grid3, grid4, grid5, grid6);

		params.add("L", 200);
		params.add("la",10);    // scale of the bias dilution region
		params.add("lb",10); 
		
		params.add("NJ",-4.0);
	    params.add("deadsites");

		params.add("percent", 0.111);
		params.add("biaspercent", 0.111);
		params.add("totalruns",20);     //the number of total intervention runs
		 

		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));

		//params.add("Dseed",1);    //seed for dilution configuration
		//params.add("Sseed",1);    //seed for spin flip
		
		params.addm("T", 0.826);
		params.addm("H", 0.18);
		params.add("Emcs");    //MCS time for evolution
		params.add("Imcs");     //MCS clock for each intervention run
		
		params.add("runs");    //intervention run number   
		params.add("grow");
		params.add("decay");
		//params.add("copies");     
		params.add("magnetization");
		params.add("Dropletsize");
		params.add("copies");    //ensemble copy for droplet distribution
		

	}
	
	public void testrun(IsingStructure ising)
	{
		Random trand= new Random(1);
		for(int tstep=0; tstep<9999999; tstep++)
		{
			T=params.fget("T");
			H=params.fget("H");
			ising.MCS(T, H, trand, 1, dynamics);
			Job.animate();
			params.set("Emcs", tstep);
			params.set("magnetization", ising.magnetization);
		}
	}
	
	
	public void Properh(IsingStructure ising, Random rand, double T, double minH, double maxH, double dH)  // output the time of nucleation for different h, threshold=90%
	{
		
		String run="<T="+fmt.format(T*1000)+", L= "+fmt.format(L)+", la= "+fmt.format(la)+", lb= "+fmt.format(lb)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String path = "/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/scanfield "+run+".txt";
		for(double field=maxH; field>minH; field-=dH)
		{
			String pic="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/pic/<H= "+fmt.format(field*1000)+">"+run;
			
			Evolution=ising.clone();
			Job.animate();
			Erand=rand.clone();
			params.set("H", field);
			//Evolve(Evolution, 90, Erand, T, field);
			for(int pres=0; pres<90; pres++)
			{
				
				Evolution.MCS(T, field, Erand, 1, dynamics);
				Job.animate();
				params.set("Emcs", pres);
				params.set("magnetization", ising.magnetization);
			}
			
			double totalM=0;
			for(int ps=0; ps<10; ps++)
			{
				totalM+=Evolution.magnetization;
				Evolution.MCS(T, field, Erand, 1, dynamics);
				Job.animate();
				params.set("Emcs", ps+90);
				params.set("magnetization", Evolution.magnetization);
			}
			double Ms=totalM/10;    //calculate the saturate magnetization
			params.set("H", -field);//flip the field;
			int ss=0;
			for(ss=0; Evolution.magnetization>(Ms*0.9);ss++)
			{
				Evolution.MCS(T, -field, Erand, 1, dynamics);
				Job.animate();
				params.set("Emcs", ss);
				params.set("magnetization", Evolution.magnetization);
			}
			PrintUtil.printlnToFile(path , field , ss);
			Tools.Picture(grid4, ss, ss, pic);
		}
		
	}
	
	public void Singlerun(IsingStructure ising, Random rand, double T, double H)
	{
		String singlerun="<T="+fmt.format(T*1000)+", H="+fmt.format(H*1000)+", la= "+fmt.format(la)+", lb= "+fmt.format(lb)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String singlepath = "/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/singlerun "+singlerun+".txt";
		String singlepic="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/singlerunpic/"+singlerun;
		
		Evolution=ising.clone();
		Job.animate();
		Erand=rand.clone();
		params.set("H", H);
		//Evolve(Evolution, 90, Erand, T, H);
		for(int pres=0; pres<90; pres++)
		{
			
			Evolution.MCS(T, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", pres);
			params.set("magnetization", Evolution.magnetization);
		}
		
		
		
		double totalM=0;
		for(int ps=0; ps<10; ps++)
		{
			totalM+=Evolution.magnetization;
			Evolution.MCS(T, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", ps+90);
			params.set("magnetization", Evolution.magnetization);
		}
		double Ms=totalM/10;    //calculate the saturate magnetization
		params.set("H", -H);//flip the field;
		int ss=0;
		for(ss=0; Evolution.magnetization>-(Ms*0.5);ss++)
		{
			Evolution.MCS(T, -H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", ss);
			params.set("magnetization", Evolution.magnetization);
			PrintUtil.printlnToFile(singlepath , ss , Evolution.magnetization/Ms);
			if(ss%2000==0)
			{
				Tools.Picture(grid4, ss, (int)(H*1000), singlepic);
			}
		}
		
		
	}
	
	public void Intervention(IsingStructure ising, Random rand, double T, double H, int breakpoint, int steplimit)
	{
		String Irun="<T="+fmt.format(T*1000)+", H="+fmt.format(H*1000)+", la= "+fmt.format(la)+", lb= "+fmt.format(lb)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String Ipath = "/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Intervention log.txt";
		String Ipic="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Interventionpic/<t= "+fmt.format(breakpoint)+">"+Irun;
		String Bpic="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Breakpoint/"+Irun;
		String Dpic="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Droplet/"+Irun;
		
		Evolution=ising.clone();
		Job.animate();
		Erand=rand.clone();
		params.set("H", H);
		//Evolve(Evolution, 90, Erand, T, H);
		for(int pres=0; pres<90; pres++)
		{
			
			Evolution.MCS(T, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", pres);
			params.set("magnetization", Evolution.magnetization);
		}
		
		
		double totalM=0;
		for(int ps=0; ps<10; ps++)
		{
			totalM+=Evolution.magnetization;
			Evolution.MCS(T, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", ps+90);
			params.set("magnetization", Evolution.magnetization);
		}
		double Ms=totalM/10;    //calculate the saturate magnetization
		params.set("H", -H);//flip the field;
		int ss=0;
		for(ss=0; ss<breakpoint;ss++)
		{
			Evolution.MCS(T, -H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", ss);
			params.set("magnetization", Evolution.magnetization);
		}
		Evolution.MCS(T, -H, Erand, 1, dynamics);
		params.set("magnetization", Evolution.magnetization);
		Job.animate();
		double Mc=Evolution.magnetization;
		
		
		
		
		Tools.Picture(grid4, breakpoint, 9999, Bpic);        //the snapshot at the intervention point
		{                       
			Droplet=new Percolation(Evolution,2);            //percolation mapping to determine the droplet
			Droplet.SetProbability(1);
			Droplet.fastNNMapping(47);
			dropletsize=Droplet.CS.maximumsize;
			for(int jj=0; jj<ising.M; jj++)
			{
				Istemp.spin[jj]=-1;
				if(Evolution.spin[jj]==0)
					Istemp.spin[jj]=0;
				if(Droplet.CS.set[Droplet.CS.maximumpin].lattice[jj]==2)
					Istemp.spin[jj]=2;
				
			}
			params.set("Dropletsize",dropletsize);
			Job.animate();
		}
		Tools.Picture(grid2, breakpoint, dropletsize, Dpic);        //the snapshot at the intervention point
		
		
		
		
		//now is the intervention time
		grownumber=0;
		decaynumber=0;
		
		for(int c=0; c<totalruns; c++)
		{
			params.set("runs", c+1);
			Intervention=Evolution.clone();
			Random irand=new Random(c+99);
			for(int is=0; is<steplimit; is++)
			{
				Intervention.MCS(T, -H, irand, 1, dynamics);
				Job.animate();
				params.set("Imcs", is);
				params.set("magnetization", Intervention.magnetization);
			}
			if(Intervention.magnetization>threshold*Mc)
			{
				decaynumber++;
				params.set("decay", decaynumber);
				Tools.Picture(grid3, c+1, 1111, Ipic);     //1111--decay
			}
			else
			{
				grownumber++;
				params.set("grow", grownumber);
				Tools.Picture(grid3, c+1, 8888, Ipic);     //8888--grow
			}
			
		}
		PrintUtil.printlnToFile(Ipath , Irun);
		PrintUtil.printlnToFile(Ipath , "breakpoint=  ",breakpoint);
		PrintUtil.printlnToFile(Ipath , "decay =  ", decaynumber);
		PrintUtil.printlnToFile(Ipath , "grow =  ", grownumber);
		PrintUtil.printlnToFile(Ipath , "droplet size =  ", dropletsize);
		PrintUtil.printlnToFile(Ipath , "magnetization =  ", Evolution.magnetization);
		PrintUtil.printlnToFile(Ipath , "ratio =  ", Evolution.magnetization/Ms);
		PrintUtil.printlnToFile(Ipath , "    ");
	
	
	
	}
	
	public void Dropletdistribution(IsingStructure ising, double T, double H, int copies, double thresholdM)
	{
		nucleationevents=new double[copies];
		
		
		String Drun="<T="+fmt.format(T*1000)+", H="+fmt.format(H*1000)+", la= "+fmt.format(la)+", lb= "+fmt.format(lb)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String Dpath = "/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/nucleation log.txt";
		String ddpic="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/map/"+Drun;
		String Npic="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/nucleationevents/"+Drun;
		
		for(int cc=0; cc<copies; cc++)
		{
			
			
			Evolution=ising.clone();
			Random crand=new Random(cc+99);
			params.set("H", H);
			params.set("copies", cc+1);
			
			for(int pres=0; pres<90; pres++)
			{
				Evolution.MCS(T, H, crand, 1, dynamics);
				Job.animate();
				params.set("Emcs", pres);
				params.set("magnetization", Evolution.magnetization);
			}
		
			double totalM=0;
			for(int ps=0; ps<10; ps++)
			{
				totalM+=Evolution.magnetization;
				Evolution.MCS(T, H, crand, 1, dynamics);
				Job.animate();
				params.set("Emcs", ps+90);
				params.set("magnetization", Evolution.magnetization);
			}
			
			double Ms=totalM/10;    //calculate the saturate magnetization
			params.set("H", -H);//flip the field;
			int ss=0;
			for(ss=0; Evolution.magnetization>(Ms*thresholdM);ss++)
			{
				Evolution.MCS(T, -H, crand, 1, dynamics);
				Job.animate();
				params.set("Emcs", ss);
				params.set("magnetization", Evolution.magnetization);
				
			}
			Tools.Picture(grid4, ss, (int)(H*1000), Npic);
			
			PrintUtil.printlnToFile(Dpath ,cc+1, ss , Evolution.magnetization);
			nucleationevents[cc]=ss;
			
			for(int jjj=0; jjj<ising.M; jjj++)
			{
				spintotal[jjj]+=Evolution.spin[jjj];
			}
			Job.animate();
			
			Tools.Picture(grid5, cc+1, ss, ddpic);
			
		}
		
		
		
		meanNT=Tools.Mean(nucleationevents, copies);
		SDNT=Tools.SD(nucleationevents, copies, meanNT);
		
		
		PrintUtil.printlnToFile(Dpath , Drun);
		PrintUtil.printlnToFile(Dpath , "meanNT=  ",meanNT);
		PrintUtil.printlnToFile(Dpath , "SDNT =  ", SDNT);

		PrintUtil.printlnToFile(Dpath , "    ");
		
		Tools.Picture(grid5, 9999, 9999, ddpic);   //the final totalspin distribution
		
		
	}
	
	
	
	public void Multihistogram(IsingStructure ising, double T, double H, int runs, int copies, double thresholdM)
	{
		lifetime=new double[runs];
		String Mrun="<T="+fmt.format(T*1000)+", H="+fmt.format(H*1000)+", la= "+fmt.format(la)+", lb= "+fmt.format(lb)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">"+"threshold= "+fmt.format(thresholdM*1000);
		String Mpath="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Multihistogram/lifetime"+Mrun+".txt";
		//String MSDpath="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Multihistogram/SDlifetime"+Mrun+".txt";
		String Mlog="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Multihistogram/multihistogramlog.txt";
		String Mpic="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Multihistogram/"+Mrun;
		
		
		
		for(int cc=0; cc<copies; cc++)
		{
			params.set("copies", cc+1);
			spintotal=new double[ising.M];
			Job.animate();
			
			for(int rr=0; rr<runs; rr++)
			{
				Evolution=ising.Dperturbation(cc+1);
				Evolution.Sinitialization(1,Sseed);
				
				Random rrand=new Random(rr+99);
				params.set("H", H);
				params.set("runs", rr+1);
				Job.animate();
				
				for(int pres=0; pres<90; pres++)
				{
					Evolution.MCS(T, H, rrand, 1, dynamics);
					Job.animate();
					params.set("Emcs", pres);
					params.set("magnetization", Evolution.magnetization);
				}
			
				double totalM=0;
				for(int ps=0; ps<10; ps++)
				{
					totalM+=Evolution.magnetization;
					Evolution.MCS(T, H, rrand, 1, dynamics);
					Job.animate();
					params.set("Emcs", ps+90);
					params.set("magnetization", Evolution.magnetization);
				}
				
				double Ms=totalM/10;    //calculate the saturate magnetization
				params.set("H", -H);//flip the field;
				int ss=0;
				for(ss=0; Evolution.magnetization>(Ms*thresholdM);ss++)
				{
					Evolution.MCS(T, -H, rrand, 1, dynamics);
					Job.animate();
					params.set("Emcs", ss);
					params.set("magnetization", Evolution.magnetization);
					
				}
				
				lifetime[rr]=ss;
				
				for(int jjj=0; jjj<ising.M; jjj++)
				{
					spintotal[jjj]+=Evolution.spin[jjj];
				}
				Job.animate();
							
			}
			meanNT=Tools.Mean(lifetime, runs);
			SDNT=Tools.SD(lifetime, runs, meanNT);
			
			
			PrintUtil.printlnToFile(Mlog , Mrun);
			
			PrintUtil.printlnToFile(Mpath , cc+1, meanNT, SDNT);
			//PrintUtil.printlnToFile(MSDpath, cc+1, SDNT);
			
			PrintUtil.printlnToFile(Mlog , "copies=  ", cc+1);
			PrintUtil.printlnToFile(Mlog , "deadsites=  ",Evolution.deadsites);

			PrintUtil.printlnToFile(Mlog , "    ");
			
			Tools.Picture(grid5, 9999, cc+1, Mpic);   //the final totalspin distribution
			
			
			
			
		}
		
	}
	
	public void Singlehistogram(IsingStructure ising, double T, double H, int runs, double thresholdM, int DPseed)
	{
		lifetime=new double[runs];
		
		String Srun="<T="+fmt.format(T*1000)+", H="+fmt.format(H*1000)+", la= "+fmt.format(la)+", lb= "+fmt.format(lb)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">"+"Seed= "+fmt.format(DPseed);
		String Spath="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Singlehistogram/"+Srun+".txt";
		String Slog="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Singlehistogram/singlehistogramlog.txt";
		String Spic="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Singlehistogram/"+Srun;

		Job.animate();
		
		for(int rr=0; rr<runs; rr++)
		{
			if(DPseed==0)
				Evolution=ising.clone();
			else
				Evolution=ising.Dperturbation(DPseed);
			
			Evolution.Sinitialization(1,Sseed);
			
			Random rrand=new Random(rr+99);
			params.set("H", H);
			params.set("runs", rr+1);
			Job.animate();
			
			for(int pres=0; pres<90; pres++)
			{
				Evolution.MCS(T, H, rrand, 1, dynamics);
				Job.animate();
				params.set("Emcs", pres);
				params.set("magnetization", Evolution.magnetization);
			}
		
			double totalM=0;
			for(int ps=0; ps<10; ps++)
			{
				totalM+=Evolution.magnetization;
				Evolution.MCS(T, H, rrand, 1, dynamics);
				Job.animate();
				params.set("Emcs", ps+90);
				params.set("magnetization", Evolution.magnetization);
			}
			
			double Ms=totalM/10;    //calculate the saturate magnetization
			params.set("H", -H);//flip the field;
			int ss=0;
			for(ss=0; Evolution.magnetization>(Ms*thresholdM);ss++)
			{
				Evolution.MCS(T, -H, rrand, 1, dynamics);
				Job.animate();
				params.set("Emcs", ss);
				params.set("magnetization", Evolution.magnetization);
				
			}
			
			PrintUtil.printlnToFile(Spath ,rr+1, ss);
			lifetime[rr]=ss;
			
			for(int jjj=0; jjj<ising.M; jjj++)
			{
				spintotal[jjj]+=Evolution.spin[jjj];
			}
			Job.animate();
						
		}
		meanNT=Tools.Mean(lifetime, runs);
		SDNT=Tools.SD(lifetime, runs, meanNT);
		
		
		PrintUtil.printlnToFile(Slog , Srun);
		
		PrintUtil.printlnToFile(Slog , "DPseed=  ", DPseed);
		PrintUtil.printlnToFile(Slog, "threshold= ", thresholdM);
		PrintUtil.printlnToFile(Slog , "meanNT=  ",meanNT);
		PrintUtil.printlnToFile(Slog , "SDNT =  ", SDNT);
		PrintUtil.printlnToFile(Slog , "deadsites=  ",Evolution.deadsites);

		PrintUtil.printlnToFile(Slog , "    ");
		
		Tools.Picture(grid5, 9999, DPseed, Spic);   //the final totalspin distribution
		
		
		
		
	}
	
	public void Singlegrowth(IsingStructure ising, double T, double H, int runs, int DPseed, int thNumber)  //single realization of dilution, multiple runs, thNumber (default-6) is the total number of thresholds
	{
		//int thNumber=20;
		double growtime[][]=new double[thNumber][runs];
		int inttime[][]=new int[thNumber][runs];  //integer form of growtime
		int printtemp[]=new int[thNumber];
		
		double threshold[]=new double[thNumber];
		double meanNT[]=new double[thNumber];
		double SDNT[]=new double[thNumber];
		
		int snapshottarget=1;
		for(int th=0; th<thNumber; th++)
		{
			threshold[th]=0.95-0.05*th;
		}
		
		
		
		//about the array of[thNumber]: 0---95%  i---(95%-i*5%)
		
		String Srun="<T="+fmt.format(T*1000)+", H="+fmt.format(H*1000)+", la= "+fmt.format(la)+", lb= "+fmt.format(lb)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">"+"Seed= "+fmt.format(DPseed);
		String Spath="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Singlegrowth/growth"+fmt.format(thNumber)+Srun+".txt";
		String Slog="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Singlegrowth/singlegrowthlog.txt";
		String Spic="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Singlegrowth/"+Srun;

		Job.animate();
		
		for(int rr=0; rr<runs; rr++)
		{
			if(DPseed==0)
				Evolution=ising.clone();
			else
				Evolution=ising.Dperturbation(DPseed);
			
			Evolution.Sinitialization(1,Sseed);
			
			Random rrand=new Random(rr+99);
			params.set("H", H);
			params.set("runs", rr+1);
			Job.animate();
			
			for(int pres=0; pres<90; pres++)
			{
				Evolution.MCS(T, H, rrand, 1, dynamics);
				Job.animate();
				params.set("Emcs", pres);
				params.set("magnetization", Evolution.magnetization);
			}
		
			double totalM=0;
			for(int ps=0; ps<10; ps++)
			{
				totalM+=Evolution.magnetization;
				Evolution.MCS(T, H, rrand, 1, dynamics);
				Job.animate();
				params.set("Emcs", ps+90);
				params.set("magnetization", Evolution.magnetization);
			}
			
			double Ms=totalM/10;    //calculate the saturate magnetization
			params.set("H", -H);//flip the field;
			int ss=0;
			int tempin=0;
			for(ss=0; tempin<thNumber; ss++)
			{
				Evolution.MCS(T, -H, rrand, 1, dynamics);
				Job.animate();
				params.set("Emcs", ss);
				params.set("magnetization", Evolution.magnetization);
				if(Evolution.magnetization<(Ms*threshold[tempin]))
				{
					if(tempin==snapshottarget)
					{
						for(int jjj=0; jjj<ising.M; jjj++)
						{
							spintotal[jjj]+=Evolution.spin[jjj];
						}
						Job.animate();
					}	
					growtime[tempin][rr]=ss;
					inttime[tempin][rr]=ss;
					printtemp[tempin]=ss;
					tempin++;
				}
			}
			
			PrintUtil.printScalarAndVectorToFile(Spath, rr+1, printtemp);
			//PrintUtil.printlnToFile(Spath ,rr+1, inttime[0][rr], inttime[1][rr], inttime[2][rr], inttime[3][rr], inttime[4][rr], inttime[5][rr]);
			
			

						
		}
		
		for(int ppp=0; ppp<thNumber; ppp++)
		{
			meanNT[ppp]=Tools.Mean(growtime[ppp], runs);
			SDNT[ppp]=Tools.SD(growtime[ppp], runs, meanNT[ppp]);
		}
		
		
		
		PrintUtil.printlnToFile(Slog , Srun);
		
		PrintUtil.printlnToFile(Slog , "DPseed=  ", DPseed);
		
		PrintUtil.printlnToFile(Slog , "threshold=  ");
		PrintUtil.printScalarAndVectorToFile(Slog, 0, threshold);
		//PrintUtil.printlnToFile(Slog , threshold[0], threshold[1], threshold[2], threshold[3], threshold[4], threshold[5]);
		
		PrintUtil.printlnToFile(Slog , "meanNT=  ");
		PrintUtil.printScalarAndVectorToFile(Slog, 0, meanNT);
		//PrintUtil.printlnToFile(Slog , meanNT[0], meanNT[1], meanNT[2], meanNT[3], meanNT[4], meanNT[5]);
		
		PrintUtil.printlnToFile(Slog , "SDNT =  ");
		PrintUtil.printScalarAndVectorToFile(Slog, 0, SDNT);
		//PrintUtil.printlnToFile(Slog , SDNT[0], SDNT[1], SDNT[2], SDNT[3], SDNT[4], SDNT[5]);
		
		PrintUtil.printlnToFile(Slog , "deadsites=  ",Evolution.deadsites);

		PrintUtil.printlnToFile(Slog , "    ");
		
		Tools.Picture(grid5, 9999, DPseed, Spic);   //the final totalspin distribution
		
		
		
	}
	
	public void Multigrowth(IsingStructure ising, double T, double H, int copies, int thNumber)  //single realization of dilution, multiple runs, thNumber (default-6) is the total number of thresholds
	{
		//int thNumber=20;
		double growtime[][]=new double[thNumber][copies];
		int inttime[][]=new int[thNumber][copies];  //integer form of growtime
		int printtemp[]=new int[thNumber];
		
		double threshold[]=new double[thNumber];
		double meanNT[]=new double[thNumber];
		double SDNT[]=new double[thNumber];
		
		int snapshottarget=1;
		for(int th=0; th<thNumber; th++)
		{
			threshold[th]=0.95-0.05*th;
		}
		
		
		
		//about the array of[thNumber]: 0---95%  i---(95%-i*5%)
		
		String Mrun="<T="+fmt.format(T*1000)+", H="+fmt.format(H*1000)+", la= "+fmt.format(la)+", lb= "+fmt.format(lb)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String Mpath="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Multigrowth/growth"+fmt.format(thNumber)+Mrun+".txt";
		String Mlog="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Multigrowth/Multigrowthlog.txt";
		String Mpic="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Multigrowth/"+Mrun;

		Job.animate();
		
		for(int cc=0; cc<copies; cc++)
		{
			
			Evolution=ising.Dperturbation(cc+1);
			
			Evolution.Sinitialization(1,Sseed);
			
			Random rrand=new Random(cc+99);
			params.set("H", H);
			params.set("copies", cc+1);
			params.set("deadsites", Evolution.deadsites);
			Job.animate();
			
			for(int pres=0; pres<90; pres++)
			{
				Evolution.MCS(T, H, rrand, 1, dynamics);
				Job.animate();
				params.set("Emcs", pres);
				params.set("magnetization", Evolution.magnetization);
			}
		
			double totalM=0;
			for(int ps=0; ps<10; ps++)
			{
				totalM+=Evolution.magnetization;
				Evolution.MCS(T, H, rrand, 1, dynamics);
				Job.animate();
				params.set("Emcs", ps+90);
				params.set("magnetization", Evolution.magnetization);
			}
			
			double Ms=totalM/10;    //calculate the saturate magnetization
			params.set("H", -H);//flip the field;
			int ss=0;
			int tempin=0;
			for(ss=0; tempin<thNumber; ss++)
			{
				Evolution.MCS(T, -H, rrand, 1, dynamics);
				Job.animate();
				params.set("Emcs", ss);
				params.set("magnetization", Evolution.magnetization);
				if(Evolution.magnetization<(Ms*threshold[tempin]))
				{
					if(tempin==snapshottarget)
					{
						for(int jjj=0; jjj<ising.M; jjj++)
						{
							spintotal[jjj]+=Evolution.spin[jjj];
						}
						Job.animate();
					}	
					growtime[tempin][cc]=ss;
					inttime[tempin][cc]=ss;
					printtemp[tempin]=ss;
					tempin++;
				}
			}
			
			PrintUtil.printScalarAndVectorToFile(Mpath, cc+1, printtemp);
			
			
			

						
		}
		
		for(int ppp=0; ppp<thNumber; ppp++)
		{
			meanNT[ppp]=Tools.Mean(growtime[ppp], copies);
			SDNT[ppp]=Tools.SD(growtime[ppp], copies, meanNT[ppp]);
		}
		
		
		
		PrintUtil.printlnToFile(Mlog , Mrun);
		
		
		PrintUtil.printlnToFile(Mlog , "threshold=  ");
		PrintUtil.printScalarAndVectorToFile(Mlog, 0, threshold);
		//PrintUtil.printlnToFile(Slog , threshold[0], threshold[1], threshold[2], threshold[3], threshold[4], threshold[5]);
		
		PrintUtil.printlnToFile(Mlog , "meanNT=  ");
		PrintUtil.printScalarAndVectorToFile(Mlog, 0, meanNT);
		//PrintUtil.printlnToFile(Slog , meanNT[0], meanNT[1], meanNT[2], meanNT[3], meanNT[4], meanNT[5]);
		
		PrintUtil.printlnToFile(Mlog , "SDNT =  ");
		PrintUtil.printScalarAndVectorToFile(Mlog, 0, SDNT);
		//PrintUtil.printlnToFile(Slog , SDNT[0], SDNT[1], SDNT[2], SDNT[3], SDNT[4], SDNT[5]);
		
		PrintUtil.printlnToFile(Mlog , "deadsites=  ",Evolution.deadsites);

		PrintUtil.printlnToFile(Mlog , "    ");
		
		Tools.Picture(grid5, 9999, snapshottarget, Mpic);   //the final totalspin distribution
		
		
		
	}
	
	public void Surface(IsingStructure ising, Random rand, double T, double H, int steplength)
	{
		String singlerun="<T="+fmt.format(T*1000)+", H="+fmt.format(H*1000)+", la= "+fmt.format(la)+", lb= "+fmt.format(lb)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String singlepath = "/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Surface/"+singlerun+".txt";
		String singlepic="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Surface/pic/"+singlerun;
	    surfacedilution=new int[ising.M];
		surfacemetaspin=new int[ising.M];
		surfacemap=new int [ising.M];  //the map to record the configuration of surface of the droplet: 0-background 1-surface
			
		
		
		
		
		Evolution=ising.clone();
		Job.animate();
		Erand=rand.clone();
		params.set("H", H);
		//Evolve(Evolution, 90, Erand, T, H);
		for(int pres=0; pres<90; pres++)
		{
			
			Evolution.MCS(T, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", pres);
			params.set("magnetization", Evolution.magnetization);
		}
		
		
		
		double totalM=0;
		for(int ps=0; ps<10; ps++)
		{
			totalM+=Evolution.magnetization;
			Evolution.MCS(T, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", ps+90);
			params.set("magnetization", Evolution.magnetization);
		}
		double Ms=totalM/10;    //calculate the saturate magnetization
		params.set("H", -H);//flip the field;
		for(int pp=0; Evolution.magnetization>(Ms*0.97); pp++)
		{
			Evolution.MCS(T, -H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", -pp);
			params.set("magnetization", Evolution.magnetization);
		}
		
		
		for(int ss=0; Evolution.magnetization>-(Ms*0.0);ss++)
		{
			Evolution.MCS(T, -H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", ss);
			params.set("magnetization", Evolution.magnetization);
			PrintUtil.printlnToFile(singlepath , ss , Evolution.magnetization/Ms);
			if(ss%steplength==0)
			{
				Tools.Picture(grid4, ss, (int)(H*1000), singlepic);
				
				
				{                       
					Droplet=new Percolation(Evolution,2);            //percolation mapping to determine the droplet
					Droplet.SetProbability(1);
					Droplet.fastNNMapping(47);
					dropletsize=Droplet.CS.maximumsize;
					surfacedilution=Droplet.CS.set[Droplet.CS.maximumpin].Surface(Evolution.spin,0);
					surfacemetaspin=Droplet.CS.set[Droplet.CS.maximumpin].Surface(Evolution.spin,1);
					
					
					for(int jj=0; jj<ising.M; jj++)
					{
						Istemp.spin[jj]=-1;
						if(Evolution.spin[jj]==0)
							Istemp.spin[jj]=0;
						if(Droplet.CS.set[Droplet.CS.maximumpin].lattice[jj]==2)
							Istemp.spin[jj]=2;
						
					}
					params.set("Dropletsize",dropletsize);
					Job.animate();
				}
				
				
				
			}
		}
	}
	
	public void run(){
		
		
		L = (int)params.fget("L");
		la = (int)params.fget("la");
		lb = (int)params.fget("lb");
		M = L * L;
		NJ = params.fget("NJ");

		percent=params.fget("percent");
		biaspercent=params.fget("biaspercent");
		dynamics= params.sget("Dynamics");
		totalruns=(int)params.fget("totalruns");
		Dseed = 1;
		Bseed = 1;
		Sseed = 1;

		
	    IS=new IsingStructure(L,L,0,NJ,percent,biaspercent,"square");   
	    Istemp=new IsingStructure(L,L,0,NJ,percent,biaspercent,"square");
	    Intervention=new IsingStructure(L,L,0,NJ,percent,biaspercent,"square");
	    Evolution=new IsingStructure(L,L,0,NJ,percent,biaspercent,"square");
	    Droplet=new Percolation();
	    spintotal= new double[IS.M];
	    dilutionmap= new double[IS.M];
	   
	    
	    
	    Tools=new BasicTools();
	    T=params.fget("T");
	    H=params.fget("H");
	    
	    {//initialization
	    	
	    	IS.Dinitialization(Dseed, Bseed, la, lb);
	    	params.set("deadsites",IS.deadsites);
	    	IS.Sinitialization(1, Sseed);
	    	Istemp=IS.clone();
	    	//Istemp.Sinitialization(1, Sseed);
	    	Intervention=IS.clone();
	    }
	    
	    Job.animate();
	    
	    //Random rand=new Random(Sseed);
	    
	    //testrun(Istemp);
	    
	    //Properh(IS, rand, T, 0.01, 0.6, 0.01);
	    
	    //Singlerun(IS, rand, T, H);
        
	    
	    breakpoint= 9335;
	    threshold= 0.99;
	    steplimit= 2000;
	    //Intervention(IS, rand, T, H, breakpoint, steplimit);
      
	    //Dropletdistribution(IS, T, H, 100, 0.97);
	    
	    //Singlehistogram(IS, T, H, 500, 0.9, 1);
	    
	    
	    
	    Singlegrowth(IS, T, H, 500, 4, 20);
	    
	    Singlegrowth(IS, T, H, 500, 5, 20);
	    //Multigrowth(IS, T, H, 500, 20);

	    //Multihistogram(IS, T, H, 10, 500, 0.9);
	    
	    //Multihistogram(IS, T, H, 10, 500, 0.8);
	    
	    //Multihistogram(IS, T, H, 10, 500, 0.7);
	    
	    //Multihistogram(IS, T, H, 10, 500, 0.5);
	    
	    
	    
	    
	    

	}
	
	
	
}