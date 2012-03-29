package kang.ising;

import java.awt.Color;

import java.text.DecimalFormat;



import kang.util.PrintUtil;
import chris.util.Random;


import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;




import kang.ising.BasicStructure.J1J2Structure;
import kang.ising.BasicStructure.BasicTools;


public class J1J2Nucleation extends Simulation{
	

	Grid grid1=new Grid("Spin simulation");     // the map to display the simulation
	Grid grid2=new Grid("Display");
	Grid grid3=new Grid("Vertical vs Horizontal");
	Grid grid4=new Grid("Order vs Disorder");
	
	public J1J2Structure JJS;
	public J1J2Structure JJstemp;

	public Random Erand;
	
	public BasicTools Tools;
	
	
	//initialization parameters
	public int L,la,lb;
	public double M;
	public double NJ1, NJ2;
	public int Dseed, Bseed, Sseed;
	public double percent;
	public double biaspercent;
	public int deadsite;
	public String dynamics;
	
	
	//dynamic parameters
	public double T;
	public double h;   //the amplitude for the field
	public double[] H;
	
	public void animate()
	{
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);      //up spin
		ising.setColor(-1, Color.WHITE);     //down spin
		ising.setColor(0, Color.RED);        //normal dilution
		ising.setColor(2, Color.BLUE);       //clusters
		ising.setColor(-2, Color.GREEN);     //
		ising.setColor(3, Color.darkGray);    // the centers of the clusters
			
		
		//J1-J2 color code
		ColorPalette JJising = new ColorPalette ();
		JJising.setColor(3, Color.BLACK);      
		JJising.setColor(-1, Color.WHITE);    
		JJising.setColor(1, Color.RED);        
		JJising.setColor(2, Color.BLUE);       
		JJising.setColor(4, Color.GREEN);    
		JJising.setColor(-2, Color.darkGray);  
		JJising.setColor(0, Color.CYAN);

		//J1-J2 color code vertical vs Horizontal
		ColorPalette VHJising = new ColorPalette ();
		VHJising.setColor(3, Color.BLUE);      
		VHJising.setColor(-1, Color.WHITE);    
		VHJising.setColor(1, Color.RED);        
		VHJising.setColor(2, Color.BLUE);       
		VHJising.setColor(4, Color.RED);    
		VHJising.setColor(-2, Color.darkGray);  
		VHJising.setColor(0, Color.CYAN);
		
		//J1-J2 color code Order vs Disorder
		ColorPalette ODJising = new ColorPalette ();
		ODJising.setColor(3, Color.BLUE);      
		ODJising.setColor(-1, Color.WHITE);    
		ODJising.setColor(1, Color.BLUE);        
		ODJising.setColor(2, Color.BLUE);       
		ODJising.setColor(4, Color.BLUE);    
		ODJising.setColor(-2, Color.darkGray);  
		ODJising.setColor(0, Color.CYAN);
		
		
		grid1.setColors(ising);
		grid1.registerData(L, L, JJstemp.spin);
		grid2.setColors(JJising);
		grid2.registerData(L, L, JJstemp.display);
		grid3.setColors(VHJising);
		grid3.registerData(L, L, JJstemp.display);
		grid4.setColors(ODJising);
		grid4.registerData(L, L, JJstemp.display);
	

	}

	public void clear()
	{
		grid1.clear();
		grid2.clear();
	
	}
	
	public static void main (String[] J1J2Nucleation){
		new Control(new J1J2Nucleation(), "Kang Liu's J1-J2 ising model's nucleation" );
	}
	
	
	
	public void load(Control J1J2Nucleation)
	{

		J1J2Nucleation.frameTogether("Display", grid1 ,grid2, grid3, grid4);

		params.add("L", 300);
		params.add("la",10);    // scale of the bias dilution region
		params.add("lb",10); 
		
		params.add("NJ1",-4.0);     //ferromagnetic NJ1
		params.add("NJ2", 2.2);      //antiferromagnetic NJ2  
		params.add("g", 0.55);
	    params.add("deadsites");

		params.add("percent", 0.0);
		params.add("biaspercent", 0.0);
		
		//params.add("totalruns",20);     //the number of total intervention runs
		 

		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));

		//params.add("Dseed",1);    //seed for dilution configuration
		//params.add("Sseed",1);    //seed for spin flip
		
		params.addm("T", 0.826);
		params.addm("h", 0.0);
		params.add("Emcs");    //MCS time for evolution
		//params.add("Imcs");     //MCS clock for each intervention run
		
		params.add("runs");    //intervention run number   
		params.add("mx");
		params.add("my");
		    
		params.add("magnetization");
		params.add("mm2");
		//params.add("Dropletsize");
		//params.add("copies");    //ensemble copy for droplet distribution
		

	}
	
	public void testrun(J1J2Structure jjising)
	{
		Random trand= new Random(1);
		for(int tstep=0; tstep<9999999; tstep++)
		{
			T=params.fget("T");
			h=params.fget("h");
			
			for(int hj=0; hj<L*L; hj++)
			{
				H[hj]=h;
			}
			
			jjising.MCS(T, H, trand, 1, dynamics);
			Job.animate();
			params.set("Emcs", tstep);
			params.set("magnetization", jjising.magnetization);
			params.set("mx", jjising.mx);
			params.set("my", jjising.my);
			params.set("mm2", jjising.mm2);
			
		}
	}
	
	public void run(){
		
		
		L = (int)params.fget("L");
		la = (int)params.fget("la");
		lb = (int)params.fget("lb");
		M = L * L;
		NJ1 = params.fget("NJ1");
		NJ2 = params.fget("NJ2");
		H= new double[L*L];

		

		percent=params.fget("percent");
		biaspercent=params.fget("biaspercent");
		dynamics= params.sget("Dynamics");
		
		Dseed = 1;
		Bseed = 1;
		Sseed = 1;

		
	    JJS=new J1J2Structure(L,L,NJ1,NJ2,percent,biaspercent);   
	    JJstemp=new J1J2Structure(L,L,NJ1,NJ2,percent,biaspercent);
	    
	    Tools=new BasicTools();
	    T=params.fget("T");
	    h=params.fget("h");
	    
	    {//initialization
	    	
	    	JJS.Dinitialization(Dseed, Bseed, la, lb);
	    	params.set("deadsites",JJS.deadsites);
	    	JJS.Sinitialization(0, Sseed);
	        JJstemp=JJS.clone();
	    
	    }
	    
	    Job.animate();
	    
	    //Random rand=new Random(Sseed);
	    
	    testrun(JJstemp);
	}
	
}
