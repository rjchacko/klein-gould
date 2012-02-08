package kang.ising;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import chris.util.PrintUtil;
import chris.util.Random;

import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;


import kang.ising.BasicStructure.IsingStructure;
import kang.ising.BasicStructure.BasicTools;
import kang.ising.BasicStructure.Percolation;



public class CNTdilution extends Simulation{
	
	Grid grid1=new Grid("dilution configuration");     // the map to display the dilution configuration
	Grid grid2=new Grid("simulation");     // the map to display the simulation
	Grid grid3=new Grid("intervention copy");
	Grid grid4=new Grid("evolution display"); 
	
	public IsingStructure IS;
	public IsingStructure Istemp;
	public IsingStructure Intervention;
	public IsingStructure Evolution;
	public Random Erand;
	
	public BasicTools Tools;
	
	private DecimalFormat fmt = new DecimalFormat("000");
	private DecimalFormat bmt = new DecimalFormat("0000");
	private DecimalFormat qmt = new DecimalFormat("00000");
	
	
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
		
	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
		grid3.clear();
		grid4.clear();
	}

	public static void main (String[] CNTdilution){
		new Control(new CNTdilution(), "Kang Liu's nearest neighbor site-diluted ising model's nucleation" );
	}
	
	
	
	public void load(Control CNTdilution)
	{

		CNTdilution.frameTogether("Display", grid1 ,grid2, grid3, grid4);

		params.add("L", 200);
		params.add("la",10);    // scale of the bias dilution region
		params.add("lb",10); 
		
		params.add("NJ",-4.0);
	    params.add("deadsites");

		params.add("percent", 0.0);
		params.add("biaspercent", 1);
		params.add("totalruns",20);     //the number of total intervention runs
		 

		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));

		//params.add("Dseed",1);    //seed for dilution configuration
		//params.add("Sseed",1);    //seed for spin flip
		
		params.addm("T", 1.008);
		params.addm("H", 0.33);
		
		
		
		
		params.add("Emcs");    //MCS time for evolution
		params.add("Imcs");     //MCS clock for each intervention run
		
		params.add("runs");    //intervention copy number
		params.add("grow");
		params.add("decay");
		//params.add("copies");     
		params.add("magnetization");
		params.add("Dropletsize");
		

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
			params.set("MCS", tstep);
			params.set("magnetization", ising.magnetization);
		}
	}
	
	/*public void Evolve(IsingStructure ising, int step, Random rand, double T, double H)
	{
		Random flip= rand.clone();
		for(int s=0; s<step; s++)
		{
			ising.MCS(T, H, flip, 1, dynamics);
			Job.animate();
			params.set("Emcs", s);
			params.set("magnetization", ising.magnetization);
		}
	}  // the function to run ising N steps with T and H using specific random number rand 
	*/
	
	
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
			params.set("magnetization", ising.magnetization);
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
			Tools.Picture(grid4, ss, (int)(H*1000), singlepic);
		}
		
		
		
		
	}
	
	
	public void Intervention(IsingStructure ising, Random rand, double T, double H, int breakpoint, int steplimit)
	{
		String Irun="<T="+fmt.format(T*1000)+", H="+fmt.format(H*1000)+", la= "+fmt.format(la)+", lb= "+fmt.format(lb)+", p= "+fmt.format(percent*1000)+", pb= "+bmt.format(biaspercent*1000)+">";
		String Ipath = "/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Intervention log.txt";
		String Ipic="/Users/liukang2002507/Desktop/simulation/CNTdilution/"+dynamics+"/Interventionpic/<H= "+fmt.format(H*1000)+">"+Irun;
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
			params.set("magnetization", ising.magnetization);
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
		Job.animate();
		
		
		
		
		Tools.Picture(grid4, breakpoint, 9999, Bpic);        //the snapshot at the intervention point
		{                       
			Droplet=new Percolation(Evolution,2);            //percolation mapping to determine the droplet
			Droplet.SetProbability(1);
			Droplet.fastNNMapping(47);
			dropletsize=Droplet.CS.maximumsize;
			for(int jj=0; jj<ising.M; jj++)
			{
				Istemp.spin[jj]=Evolution.spin[jj];
				if(Droplet.CS.set[Droplet.CS.maximumpin].lattice[jj]==2)
					Istemp.spin[jj]=2;
				
			}
			params.set("Dropletsize",dropletsize);
			Job.animate();
		}
		Tools.Picture(grid2, breakpoint, 9999, Dpic);        //the snapshot at the intervention point
		
		
		
		
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
			if(Intervention.magnetization>threshold*Ms)
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
	    
	    
	    Tools=new BasicTools();
	    T=params.fget("T");
	    H=params.fget("H");
	    
	    {//initialization
	    	
	    	IS.Dinitialization(Dseed, Bseed, la, lb);
	    	params.set("deadsites",IS.deadsites);
	    	IS.Sinitialization(1, Sseed);
	    	Istemp=IS.clone();
	    	Intervention=IS.clone();
	    }
	    
	    Job.animate();
	    
	    Random rand=new Random(Sseed);
	    
	    //testrun(Istemp);
	    
	    //Properh(IS, rand, T, 0.1, 0.6, 0.02);
	    
	    //Singlerun(IS, rand, T, H);

	    
	    breakpoint= 2365;
	    threshold= 0.98;
	    steplimit= 200;
	    Intervention(IS, rand, T, H, breakpoint, steplimit);
      
        
	    
	    
	    
	    
	    
	    

	}
	
	
	
}