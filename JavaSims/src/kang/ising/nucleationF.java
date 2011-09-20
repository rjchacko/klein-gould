package kang.ising;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import chris.util.Random;
import chris.util.PrintUtil;

import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.params.DoubleValue;
import kang.ising.BasicStructure.IsingBasic;
import kang.ising.BasicStructure.IsingStructure;
import kang.ising.BasicStructure.Percolation;

import scikit.graphics.ColorGradient;

import scikit.jobs.Job;
import scikit.jobs.Simulation;


public class nucleationF extends Simulation{
	
	private DecimalFormat fmt = new DecimalFormat("000000");
	
	Grid grid1 = new Grid("Grid 1");        // for the evolution of IBf.IS's evolution
	Grid grid2 = new Grid("Grid 2");        // for the intervention
	Grid grid3 = new Grid("Grid 3");        // for the movement of the centers of clusters 
	Grid grid4 = new Grid("Grid 4");        // to check the 50th step of cluster evolution
	Grid grid5 = new Grid("Grid 5");        // for display of the dilution map
	
	public int L1,L2,M,R,deadsites,Dseed,Bseed,Sseed;
	public double percent,biaspercent,NJ;
	public double T,H;
	public IsingStructure ISf;
	public IsingStructure IVcopy;
	public IsingStructure ISCE;  //for cluster evolution
	public Random CEflip;
	
	public IsingBasic IBf;
	public Percolation ISfP;
	
	
	
	public int clustermap[];    //the array for the display of the clusters
	
	
	public int step;
	public double IVstart;
	public double Ms;  //saturated magnetization using for intervention criteria
	public int decay;
	public int grow;
	
    public void Intervention(int copies, int steplimit,double threshold)
    {
    	int IVstep=0;
 
    	grow=0;
    	decay=0;
    	
    	for(int c=0;c<copies;c++)
    	{
            IVcopy= IBf.IS.clone();
    		Random Iflip= new Random(c+2);
    		
    		for(IVstep=0; IVstep<steplimit; IVstep++)
    			{
    			    IVcopy.MCS(IBf.QuenchT, IBf.QuenchH, Iflip, 1);
    			    Job.animate();
    			}
    		if((IVcopy.Magnetization())<(Ms*threshold))
    			grow++;
    		if((IVcopy.Magnetization())>(Ms*threshold))
    			decay++;
    		Job.animate();
    		movie(grid2, (int)(IVstart*100), c);
    	}
    }
    
    public boolean findnucleation()
    {
    	boolean event=false;
    	if((IBf.H<0)&(IBf.IS.totalspin<0))
    		{
    		event=true;
    		String eventpath="/Users/liukang2002507/Desktop/simulation/Fnucleation/event.txt";
    		PrintUtil.printlnToFile(eventpath, "L=", IBf.IS.L1);
    		PrintUtil.printlnToFile(eventpath, "q=", IBf.IS.percent);
    		PrintUtil.printlnToFile(eventpath, "biasp=", IBf.IS.biaspercent);
    		PrintUtil.printlnToFile(eventpath, "R=", IBf.IS.R);
    		PrintUtil.printlnToFile(eventpath, "T=", T);
    		PrintUtil.printlnToFile(eventpath, "H=", H);
    		PrintUtil.printlnToFile(eventpath, "Nucleation happens at around  ",step);
    		PrintUtil.printlnToFile(eventpath, "           ");
    		PrintUtil.printlnToFile(eventpath, "           ");
    		}
    	return event;
    }
    
	public void movie(Grid grid, int copynumber, int number)   //function to capture the grid
	{
		
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/Fnucleation/pic_"+fmt.format(copynumber)+"_"+fmt.format(number)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
	}
	
	public void clustermovie(Grid grid, int clusternumber, int step)
	{
		String SaveAs = "/Users/liukang2002507/Desktop/simulation/Fnucleation/clusters/pic_"+fmt.format(clusternumber)+"_"+fmt.format(step)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
	}
	
	public void evolution(IsingStructure EIS, Random Eflip, double T, double H)
	{

		
		for(int step=0; step<150; step++)      // step=51 is the critical droplet
		{

		    movie(grid4,6666,step); 
			Percolation P1= new Percolation(EIS,1);
			Percolation P3= new Percolation(EIS,3);
			Percolation P5= new Percolation(EIS,5);
			Percolation P10= new Percolation(EIS,10);
			
			
			P1.probability(T); 
			P1.Mapping();
			clustermap= new int [M];
			for(int o=0; o<M; o++)
			{
				clustermap[o]=-1;
			}
	        P1.CS.ClustersDisplay(clustermap);
	        P1.CS.CentersDisplay(clustermap);
	        Job.animate();
	        clustermovie(grid3, 1, step);
	        
	        
			P3.probability(T); 
			P3.Mapping();
			clustermap= new int [M];
			for(int o=0; o<M; o++)
			{
				clustermap[o]=-1;
			}
	        P3.CS.ClustersDisplay(clustermap);
	        P3.CS.CentersDisplay(clustermap);
	        Job.animate();
	        clustermovie(grid3, 3, step);
	        
	        
			P5.probability(T); 
			P5.Mapping();
			clustermap= new int [M];
			for(int o=0; o<M; o++)
			{
				clustermap[o]=-1;
			}
	        P5.CS.ClustersDisplay(clustermap);
	        P5.CS.CentersDisplay(clustermap);
	        Job.animate();
	        clustermovie(grid3, 5, step);
	        
	        
			P10.probability(T); 
			P10.Mapping();
			clustermap= new int [M];
			for(int o=0; o<M; o++)
			{
				clustermap[o]=-1;
			}
	        P10.CS.ClustersDisplay(clustermap);
	        P10.CS.CentersDisplay(clustermap);
	        Job.animate();
	        clustermovie(grid3, 10, step);
	        
	        EIS.MCS(T, H, Eflip, 0.1);
	        
		}
	  
		
		
		
	}
	
	public void animate()
	{
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);      //up spin
		ising.setColor(-1, Color.WHITE);     //down spin
		ising.setColor(0, Color.RED);        //normal dilution
		ising.setColor(2, Color.BLUE);       //clusters
		ising.setColor(-2, Color.GREEN);     //
		ising.setColor(3, Color.darkGray);    // the centers of the clusters
		
		ColorGradient heatmap = new ColorGradient();
		
		grid1.setColors(ising);
		grid1.registerData(IBf.IS.L1, IBf.IS.L2, IBf.IS.spin);
		grid2.setColors(ising);
		grid2.registerData(IVcopy.L1, IVcopy.L2, IVcopy.spin);
		grid3.setColors(ising);
		grid3.registerData(L1, L2, clustermap);
		grid4.setColors(ising);
		grid4.registerData(L1, L2, ISCE.spin);
		grid5.setColors(heatmap);
		grid5.registerData(L1, L2, ISf.dilutionmap);
		
		
		params.set("grow", grow);
		params.set("decay", decay);

	}
	public void clear()
	{
		grid1.clear();
		grid2.clear();
		grid3.clear();
		grid4.clear();
		grid5.clear();
	}
	
	public static void main (String[] nucleationF){
		new Control(new nucleationF(), "Kang Liu's Ferro nucleation" );
	}
	
	public void load(Control nucleationF){
		nucleationF.frame (grid1);
		nucleationF.frame (grid2);
		nucleationF.frame (grid3);
		nucleationF.frame (grid4);
		nucleationF.frame (grid5);


		params.add("L1", 200);
		params.add("L2", 200);
		params.add("R", 10);
		params.add("NJ",-4.0);	
		params.add("percent", new DoubleValue(0.111,0,1).withSlider());
		params.add("biaspercent", new DoubleValue(0.111,0,1).withSlider());
		params.add("deadsites");	
		//params.add("Dseed",1);
		//params.add("Bseed",1);
		//params.add("Sseed",1);
		
		params.add("Ti", 1.591);                 //set Ti=4/9Tc
		params.add("Hi", 1.021);              //set Hi>0
		params.add("Tf", 1.591);         //set Tf=Ti=4/9Tc  
		params.add("Hf", -1.021);              //set Hf<0 and Hf=-Hi
		
		params.addm("T", new DoubleValue(2, 0, 10).withSlider());
		params.addm("H", new DoubleValue(0, -2, 2).withSlider());
		
		params.add("MCS");
		params.add("magnetization");

		params.add("Intervention Start",2511.90);
		params.add("grow");
		params.add("decay");
		params.add("# of clusters");


	}
	
	   public void run(){
			
			
			L1 = (int)params.fget("L1");
			L2 = (int)params.fget("L2");
			M = L1 * L2 ;
			
			clustermap= new int[M];
			for(int o=0; o<M; o++)
			{
				clustermap[o]=-1;
			}
			
			R = (int)params.fget("R");
			NJ = params.fget("NJ");

			percent=params.fget("percent");
			biaspercent=params.fget("biaspercent");
			
			Dseed = 1;//(int)params.fget("Dseed");
			Bseed = 1;//(int)params.fget("Bseed");
			Sseed = 1;//(int)params.fget("Sseed");
			
		    ISf=new IsingStructure(L1,L2,R,NJ,percent,biaspercent,"square");
		    IVcopy=new IsingStructure(L1,L2,R,NJ,percent,biaspercent,"square");
		    ISf.Dinitialization(Dseed, Bseed, 10, 10);
		    params.set("deadsites",ISf.deadsites);
		    ISf.Sinitialization(0, Sseed);
		    
		    ISCE= ISf.clone();
		  
		    
		    T=params.fget("Ti");
		    H=params.fget("Hi");
		    
		    IBf=new IsingBasic(ISf,T,H,1);
		    IBf.QuenchT=params.fget("Tf");
		    IBf.QuenchH=params.fget("Hf");
		    IVstart=params.fget("Intervention Start");
		    Job.animate();
		    
		    
		    //movie(grid5, 8888, (int)(ISf.percent*100));      //draw the heatmap
		    
		    // now run the system at high temperature for a while
		    for(int heatstep=0; heatstep<20; heatstep++)
		    {
		    	IBf.T=9;
		    	IBf.H=0;
		    	IBf.IS.MCS(IBf.T, IBf.H,IBf.flip, 1);
		    	Job.animate();
		    	params.set("magnetization", IBf.IS.magnetization);

		    }
		    
		    // now run the system at 4/9 Tc with h until the field flipping
		    params.set("T",IBf.InitialT);         
		    params.set("H",IBf.InitialH);
		    for(int prestep=0; prestep<100; prestep++)
		    {
			    IBf.T=params.fget("T");
			    IBf.H=params.fget("H");
				IBf.IS.MCS(IBf.T, IBf.H,IBf.flip, 1);
				Job.animate();
				params.set("MCS", prestep-100);
				params.set("magnetization", IBf.IS.magnetization);

		    }
		    Ms=IBf.IS.Magnetization();      //now record the saturated magnetization
		    
		    
		    // now keep the temperature at constant and flip the field to begin the nucleation
		    params.set("T",IBf.QuenchT);
		    params.set("H",IBf.QuenchH);
		    
		    for(step=0; (step<=(int)IVstart)&(!findnucleation()); step++)
		    {
		    	IBf.T=params.fget("T");
			    IBf.H=params.fget("H");
			    if(step==(int)IVstart-4)   //prepare for the cluster evolution
			    {
			    	CEflip=IBf.flip.clone();    //copy the Random number
			    	ISCE=IBf.IS.clone();
			    }
				IBf.IS.MCS(IBf.T, IBf.H,IBf.flip, 1);
				Job.animate();
				params.set("MCS", step);
				params.set("magnetization", IBf.IS.magnetization);

		    }
		    
		    if(IVstart!=(int)IVstart)      //intervention starts at rational MCS, now run down the residue step
		    {
		    	IBf.T=params.fget("T");
		    	IBf.H=params.fget("H");
			    IBf.IS.MCS(IBf.T, IBf.H, IBf.flip, (IVstart-(int)IVstart));
			    
			    ISCE.MCS(IBf.T, IBf.H, CEflip, (IVstart-(int)IVstart));    //evolve the clusters to exactly IVstart-5 MCS
		    }
		    Job.animate();
		    movie(grid1,9999,(int)(IVstart*100));    //record the critical droplet
		    
		    
		    // now make the copy of IBf.IS and run the intervention
		    //Intervention(20, 30, 0.7);
		    
		        

		    
		    evolution(ISCE, CEflip, IBf.T, IBf.H);
		   
		    /*ISfP=new Percolation(IBf.IS,1);
		    ISfP.probability(IBf.T);
		    ISfP.Mapping();
		    ISfP.CS.ClustersDisplay(clustermap);
		    ISfP.CS.CentersDisplay(clustermap);
		    Job.animate(); 
		    */
		    
	
	
	
	   }
	
}