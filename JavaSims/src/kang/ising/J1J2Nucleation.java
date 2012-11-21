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




import kang.ising.BasicStructure.IsingStructure;
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
	public int rseed;
	private DecimalFormat fmt = new DecimalFormat("000");
	private DecimalFormat pmt = new DecimalFormat("0000");
	public BasicTools Tools;
	
	
	//initialization parameters
	public int L,la,lb;
	public double M;
	public double NJ1, NJ2;
	public double g;
	public int Dseed, Bseed, Sseed;
	public double percent;
	public double biaspercent;
	public int deadsite;
	public String dynamics;
	
	
	//dynamic parameters
	public double T;
	public double h;   //the amplitude for the field
	public double hx, hy;
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
		grid3.clear();
		grid4.clear();
	
	}
	
	public static void main (String[] J1J2Nucleation){
		new Control(new J1J2Nucleation(), "Kang Liu's J1-J2 ising model's nucleation" );
	}
	
	
	
	public void load(Control J1J2Nucleation)
	{

		J1J2Nucleation.frameTogether("Display", grid1 ,grid2, grid3, grid4);

		params.add("L", 200);
		params.add("la",10);    // scale of the bias dilution region
		params.add("lb",10); 
		
		params.add("NJ1",-4.0);     //ferromagnetic NJ1
		params.add("NJ2", 2.2);      //antiferromagnetic NJ2  
		params.add("g", 0.0);
	    params.add("deadsites");

		params.add("percent", 0.0);
		params.add("biaspercent", 0.0);
		
		//params.add("totalruns",20);     //the number of total intervention runs
		 

		
		params.addm("Dynamics", new ChoiceValue("Metropolis","Glauber"));

		//params.add("Dseed",1);    //seed for dilution configuration
		//params.add("Sseed",1);    //seed for spin flip
		
		params.addm("T", 0.343);
		params.addm("h", 0.0);
		params.addm("hx", 0.0);
		params.addm("hy", 0.0);
		params.add("Emcs");    //MCS time for evolution
		//params.add("Imcs");     //MCS clock for each intervention run
		
		params.add("runs");    //intervention run number   
		params.add("mx");
		params.add("my");
		    
		params.add("magnetization");
		params.add("mm2");
		params.add("InteractionE");
		//params.add("Dropletsize");
		//params.add("copies");    //ensemble copy for droplet distribution
		

	}
	
	public void Getfield(J1J2Structure jjising, double h, double hx, double hy)
	{
		for(int hj=0; hj<L*L; hj++)
		{
			if(hx!=0)
				H[hj]=hx*jjising.Xsign(hj);
			else if(hy!=0)
				H[hj]=hy*jjising.Ysign(hj);
			else
				H[hj]=h;
			
		}
	}
		
	public void testrun(J1J2Structure jjising)
	{
		Random trand= new Random(1);
		for(int tstep=0; tstep<9999999; tstep++)
		{
			T=params.fget("T");
			h=params.fget("h");
			hx=params.fget("hx");
			hy=params.fget("hy");
			
			Getfield(jjising, h, hx, hy);
			
			jjising.MCS(T, H, trand, 1, dynamics);
			Job.animate();
			params.set("Emcs", tstep);
			params.set("magnetization", jjising.magnetization);
			params.set("mx", jjising.mx);
			params.set("my", jjising.my);
			params.set("mm2", jjising.mm2);
			params.set("InteractionE",jjising.totalintenergy);
			
		}
	}
	
	public void Multipleruns(J1J2Structure jjising, Random rand, double Ti, double Tf, double dT)
	{
		for(double t=Ti; t>Tf; t-=dT)
		{
			Singlerun(jjising, rand, 9, t);
		}
	}
	

	public void Singlerun(J1J2Structure jjising, Random rand, double Ti, double Tf)
	{
		String singlerun="g="+fmt.format(g*1000)+" L= "+fmt.format(L) +"<Ti="+pmt.format(Ti*10000)+", Tf="+pmt.format(Tf*10000)+">"+"seed"+fmt.format(rseed);
		String singlepath = "/Users/liukang2002507/Desktop/simulation/J1J2/"+dynamics+"/singlerun "+singlerun+".txt";
		String singlepic="/Users/liukang2002507/Desktop/simulation/J1J2/"+dynamics+"/singlerunpic/"+singlerun;
		
		Job.animate();
		Erand=rand.clone();
		params.set("T", Ti);
		
		for(int pres=0; pres<90; pres++)
		{
			
			jjising.MCS(Ti, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", pres);
			params.set("magnetization", jjising.magnetization);
			params.set("mx", jjising.mx);
			params.set("my", jjising.my);
			params.set("mm2", jjising.mm2);
			params.set("InteractionE",jjising.totalintenergy);
			
		}
		
		double totalM=0;
		
		for(int ps=0; ps<10; ps++)
		{
			totalM+=jjising.magnetization;
			jjising.MCS(Ti, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", ps+90);
			params.set("magnetization", jjising.magnetization);
			params.set("mx", jjising.mx);
			params.set("my", jjising.my);
			params.set("mm2", jjising.mm2);
			params.set("InteractionE",jjising.totalintenergy);
			
		}
		//double Ms=totalM/10;    //calculate the saturate magnetization
		
		params.set("T", Tf);//flip the field;
		int ss=0;
		for(ss=0; (jjising.mm2<0.7)&(ss<200000);ss++)
		{
			jjising.MCS(Tf, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", ss);
			params.set("magnetization", jjising.magnetization);
			params.set("mx", jjising.mx);
			params.set("my", jjising.my);
			params.set("mm2", jjising.mm2);
			params.set("InteractionE",jjising.totalintenergy);
			
			PrintUtil.printlnToFile(singlepath , ss , jjising.magnetization, jjising.mx,jjising.my,jjising.mm2, jjising.totalintenergy);
			if(ss%20==0)
			{
				Tools.Picture(grid2, ss, (int)(Tf*1000), singlepic);
			}
		}
		
		
	}
	
	public void XtoX(J1J2Structure jjising, Random rand, double T, double Hx)
	{
		String singlerun="g="+fmt.format(g*1000)+" L= "+fmt.format(L) +"<T="+pmt.format(T*10000)+", Hx="+pmt.format(Hx*10000)+">"+"seed"+fmt.format(rseed);
		String singlepath = "/Users/liukang2002507/Desktop/simulation/J1J2/"+dynamics+"/XtoX "+singlerun+".txt";
		String singlepic="/Users/liukang2002507/Desktop/simulation/J1J2/"+dynamics+"/XtoXpic/"+singlerun;
		
		Job.animate();
		Erand=rand.clone();
		params.set("T", 9);
		params.set("h", 0);
		params.set("hy", 0);
		params.set("hx", Hx);
		Getfield(jjising, 0, Hx, 0);
		
		for(int heat=0; heat<20; heat++)
		{
			jjising.MCS(9, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", -heat);
			params.set("magnetization", jjising.magnetization);
			params.set("mx", jjising.mx);
			params.set("my", jjising.my);
			params.set("mm2", jjising.mm2);
			params.set("InteractionE",jjising.totalintenergy);
		}
		
		params.set("T", T);
		for(int pres=0; pres<90; pres++)
		{
			
			jjising.MCS(T, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", pres);
			params.set("magnetization", jjising.magnetization);
			params.set("mx", jjising.mx);
			params.set("my", jjising.my);
			params.set("mm2", jjising.mm2);
			params.set("InteractionE",jjising.totalintenergy);
			
		}
		
		double totalMx=0;
		
		for(int ps=0; ps<10; ps++)
		{
			totalMx+=jjising.mx;
			jjising.MCS(T, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", ps+90);
			params.set("magnetization", jjising.magnetization);
			params.set("mx", jjising.mx);
			params.set("my", jjising.my);
			params.set("mm2", jjising.mm2);
			params.set("InteractionE",jjising.totalintenergy);
			
		}
		
		//double Mxs=totalMx/10;    //calculate the saturate magnetization
		params.set("hx", -Hx);//flip the field;
		Getfield(jjising, 0, -Hx, 0);
		
		int ss=0;
		for(ss=0; (jjising.mx>0.0)&(ss<200000);ss++)
		{
			jjising.MCS(T, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", ss);
			params.set("magnetization", jjising.magnetization);
			params.set("mx", jjising.mx);
			params.set("my", jjising.my);
			params.set("mm2", jjising.mm2);
			params.set("InteractionE",jjising.totalintenergy);
			
			PrintUtil.printlnToFile(singlepath , ss , jjising.magnetization, jjising.mx,jjising.my,jjising.mm2, jjising.totalintenergy);
			if(ss%50==0)
			{
				Tools.Picture(grid2, ss, (int)(Hx*1000), singlepic);
			}
		}
		
		
		
		
	}
	
	public void XtoY(J1J2Structure jjising, Random rand, double T, double Hxy)
	{
		String singlerun="g="+fmt.format(g*1000)+" L= "+fmt.format(L) +"<T="+pmt.format(T*10000)+", Hxy="+pmt.format(Hxy*10000)+">"+"seed"+fmt.format(rseed);
		String singlepath = "/Users/liukang2002507/Desktop/simulation/J1J2/"+dynamics+"/XtoY "+singlerun+".txt";
		String singlepic="/Users/liukang2002507/Desktop/simulation/J1J2/"+dynamics+"/XtoYpic/"+singlerun;
		
		Job.animate();
		Erand=rand.clone();
		params.set("T", 9);
		params.set("h", 0);
		params.set("hy", 0);
		params.set("hx", Hxy);
		Getfield(jjising, 0, Hxy, 0);
		
		for(int heat=0; heat<20; heat++)
		{
			jjising.MCS(9, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", -heat);
			params.set("magnetization", jjising.magnetization);
			params.set("mx", jjising.mx);
			params.set("my", jjising.my);
			params.set("mm2", jjising.mm2);
			params.set("InteractionE",jjising.totalintenergy);
		}
		
		params.set("T", T);
		for(int pres=0; pres<90; pres++)
		{
			
			jjising.MCS(T, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", pres);
			params.set("magnetization", jjising.magnetization);
			params.set("mx", jjising.mx);
			params.set("my", jjising.my);
			params.set("mm2", jjising.mm2);
			params.set("InteractionE",jjising.totalintenergy);
			
		}
		
		double totalMx=0;
		
		for(int ps=0; ps<10; ps++)
		{
			totalMx+=jjising.mx;
			jjising.MCS(T, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", ps+90);
			params.set("magnetization", jjising.magnetization);
			params.set("mx", jjising.mx);
			params.set("my", jjising.my);
			params.set("mm2", jjising.mm2);
			params.set("InteractionE",jjising.totalintenergy);
			
		}
		
		//double Mxs=totalMx/10;    //calculate the saturate magnetization
		params.set("hx", 0);//flip the field;
		params.set("hy", Hxy);//flip the field;
		Getfield(jjising, 0, 0, Hxy);
		
		int ss=0;
		for(ss=0; (jjising.my<0.9)&(ss<200000);ss++)
		{
			jjising.MCS(T, H, Erand, 1, dynamics);
			Job.animate();
			params.set("Emcs", ss);
			params.set("magnetization", jjising.magnetization);
			params.set("mx", jjising.mx);
			params.set("my", jjising.my);
			params.set("mm2", jjising.mm2);
			params.set("InteractionE",jjising.totalintenergy);
			
			PrintUtil.printlnToFile(singlepath , ss , jjising.magnetization, jjising.mx,jjising.my,jjising.mm2, jjising.totalintenergy);
			if(ss%5==0)
			{
				Tools.Picture(grid2, ss, (int)(Hxy*1000), singlepic);
			}
		}
		
		
		
		
	}
	
	
	
	public void run(){
		
		
		L = (int)params.fget("L");
		la = (int)params.fget("la");
		lb = (int)params.fget("lb");
		M = L * L;
		NJ1 = params.fget("NJ1");
		NJ2 = params.fget("NJ2");
		if(NJ1!=0)
			{
			g=-NJ2/NJ1;
			params.set("g",g);
			}
		
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
	    hx=params.fget("hx");
		hy=params.fget("hy");
	    
	    {//initialization
	    	
	    	JJS.Dinitialization(Dseed, Bseed, la, lb);
	    	params.set("deadsites",JJS.deadsites);
	    	JJS.Sinitialization(0, Sseed);
	        JJstemp=JJS.clone();
	    
	    }
	    
	    Job.animate();
	    
	    rseed=1;
	    Random rand=new Random(rseed);
	    
	    //testrun(JJstemp);
	    
	    //Singlerun(JJstemp, rand, 9, 0.7710);
	    //Multipleruns(JJstemp, rand, 1.2799, 1.2781, 0.0002);
	    
	    //XtoX(JJstemp, rand, T, 0.20); 
	    XtoY(JJstemp, rand, T, 0.40); 
	    
	}
	
	
}
