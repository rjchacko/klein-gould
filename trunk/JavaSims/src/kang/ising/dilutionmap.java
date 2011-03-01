package kang.ising;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.Simulation;

import kang.ising.BasicStructure.IsingStructure;
import kang.ising.BasicStructure.TemperatureField;
import chris.util.Random;


public class dilutionmap extends Simulation{
	
	Grid grid1=new Grid("evolution");
	Grid grid2=new Grid("dilutionmap");
	Grid grid3=new Grid("spintotal");
	Grid grid4=new Grid("spinmap");
	Grid grid5=new Grid("temperaturefieldmodel");
	Grid grid6=new Grid("fluctuation");
	
	private DecimalFormat fmt = new DecimalFormat("000");
	public int isingspin[];
	public int L1,L2,M;
    public int R;
    public int Dseed, Bseed, Sseed;
    public double percent,biaspercent;
    public double NJ;
    public int spintemp[];
    public double spintotal[];
    public double spinmap[];
    public double newspinmap[];
    
	public IsingStructure Initialising;
	public IsingStructure IS;
	public IsingStructure FL;  //isingstructure for the fluctuation
	public double fluctuation[];
	public TemperatureField TF;
	//public IsingStructure IScq;
	
	//public double InitialT;
	//public double InitialH;
	//public double QuenchT;
	//public double QuenchH;
	public double threshold;
	
	public void movie(Grid grid, int number, int copynumber)   //function to capture the grid
	{
		
			String SaveAs = "/Users/liukang2002507/Desktop/simulation/dilutionmap/pic_"+fmt.format(copynumber)+"_"+fmt.format(number)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
	}
	
	public void FluctuationMap(IsingStructure Ising, double T, int steplimit, int copies)
	{
		params.set("H",0);
		double fltotal[]= new double[M];
		for(int k=0;k<M; k++)
		{
			fluctuation[k]=0;
			fltotal[k]=0;
		}

		for(int c=0; c<copies; c++)
		{
			double spintotal[]= new double[M];
			double spinsqrtotal[]= new double[M];
			
			for(int i=0;i<M;i++)
			{
				spintotal[i]=0;
				spinsqrtotal[i]=0;
			}
			
			
			params.set("copies", c+1);
			IS=Ising.clone();
			Random cflip= new Random(c);
			for(int heat=0; heat<10; heat++)
			{
				params.set("T",9);     
				IS.MCS(9, 0, cflip, 1);
			    params.set("MCS", heat-10);
	     	    params.set("magnetization", IS.Magnetization());
		        Job.animate();
			}
			for(int prestep=0; prestep<2000; prestep++)
			{
				IS.MCS(T, 0, cflip, 1);
				params.set("MCS",prestep-2000);
				params.set("magnetization",IS.Magnetization());
				Job.animate();
			}
			
			for(int step=0; step<steplimit; step++)
			{
				params.set("T",T);
				IS.MCS(T, 0, cflip, 1);
				Job.animate();
				params.set("MCS",step);
				params.set("magnetization",IS.Magnetization());
				for(int j=0;j<M; j++)
				{
					if(IS.spin[j]!=0)
					{
						spintotal[j]+=IS.spin[j];
						spinsqrtotal[j]+=((IS.spin[j])*(IS.spin[j]));
						fluctuation[j]=(spinsqrtotal[j]/(step+1))-(spintotal[j]/(step+1))*(spintotal[j]/(step+1));
	                    Job.animate();				
					}
					
				}
			}
			
			
			for(int r=0; r<M; r++)
			{
				fltotal[r]+=fluctuation[r];
			}
			
			
		}
		for(int w=0; w<M; w++)
		{
			fluctuation[w]=fltotal[w]/copies;
		}
		double[] newmap= new double[M];
		newmap=takeoutzeroes(Initialising, fluctuation);
		for(int i=0; i<M; i++)
			fluctuation[i]=newmap[i];
		Job.animate();
		Random flip= new Random(1);
		for(int af=0; af<1000; af++)
		{
			IS.MCS(T, 0, flip, 1);
			Job.animate();
		}
		
		
		
	}
	
	public void CriticalTquench(IsingStructure Ising, double Ti, double Tf, double steplimit, int copies)
	{
		int limit=(int)steplimit;
		
		for(int a=0; a<M; a++)
		{
			spintemp[a]=0;
			spintotal[a]=0;
		}
		params.set("H",0);
		
		for(int c=0; c<copies; c++)
		{
			params.set("copies", c+1);
			IS=Ising.clone();
			double mag=0;
			Random cflip= new Random(c);
			for(int heat=0; heat<10; heat++)
			{
				params.set("T",Ti);     
				IS.MCS(Ti, 0, cflip, 1);
			    params.set("MCS", heat-10);
	     	    params.set("magnetization", IS.Magnetization());
		        Job.animate();
			}
			
			for(int step=0; step<limit; step++)
			{
				params.set("T",Tf);
				IS.MCS(Tf, 0, cflip, 1);
				Job.animate();
				mag=IS.Magnetization();
				params.set("MCS",step);
				params.set("magnetization",mag);
			}
			IS.MCS(Tf, 0, cflip, steplimit-limit);
			for(int j=0; j<M; j++)
			{
				spintemp[j]=IS.spin[j];
			}
			
			for(int astep=0; (astep<100)&(Math.abs(mag)<threshold); astep++)
			{
				params.set("T",Tf);
				IS.MCS(Tf, 0, cflip, 1);
				Job.animate();
				mag=IS.Magnetization();
				params.set("MCS",astep+limit);
				params.set("magnetization",mag);
			}
			
			if(mag>0)// +1 is the stable direction
			{
				for(int l=0; l<M; l++)
					spintotal[l]+=spintemp[l];
			}
			if(mag<0)  //-1 is the stable direction
			{
				for(int n=0; n<M; n++)
					spintotal[n]+=(-spintemp[n]);
			}
			
			if(c==copies-10)
			{
				newspinmap=takeoutzeroes(Initialising, spintotal);
				for(int s=0; s<M; s++)
					spinmap[s]=newspinmap[s];
				Job.animate();
			}
				
		}
		
	}
	
	public void PCriticalTquench(IsingStructure Ising, double Ti, double Tf, double steplimit, int copies)
	{
		int limit=(int)steplimit;
		
		for(int a=0; a<M; a++)
		{
			spintemp[a]=0;
			spintotal[a]=0;
		}
		params.set("H",0);
		
		for(int c=0; c<copies; c++)
		{
			params.set("copies", c+1);
			TF=new TemperatureField(Ising);
			double mag=0;
			Random cflip= new Random(c);
			for(int heat=0; heat<10; heat++)
			{
				params.set("T",Ti);     
				TF.PMCS(Ti, 0, cflip, 1);
			    params.set("MCS", heat-10);
	     	    params.set("magnetization", TF.PIS.Magnetization());
		        Job.animate();
			}
			
			for(int step=0; step<limit; step++)
			{
				params.set("T",Tf);
				TF.PMCS(Tf, 0, cflip, 1);
				Job.animate();
				mag=TF.PIS.Magnetization();
				params.set("MCS",step);
				params.set("magnetization",mag);
			}
			TF.PMCS(Tf, 0, cflip, steplimit-limit);
			for(int j=0; j<M; j++)
			{
				spintemp[j]=TF.PIS.spin[j];
			}
			
			for(int astep=0; (astep<100)&(Math.abs(mag)<threshold); astep++)
			{
				params.set("T",Tf);
				TF.PMCS(Tf, 0, cflip, 1);
				Job.animate();
				mag=TF.PIS.Magnetization();
				params.set("MCS",astep+limit);
				params.set("magnetization",mag);
			}
			
			if(mag>0)// +1 is the stable direction
			{
				for(int l=0; l<M; l++)
					spintotal[l]+=spintemp[l];
			}
			if(mag<0)  //-1 is the stable direction
			{
				for(int n=0; n<M; n++)
					spintotal[n]+=(-spintemp[n]);
			}
			
			if(c==copies-10)
			{
				newspinmap=takeoutzeroes(Initialising, spintotal);
				for(int s=0; s<M; s++)
					spinmap[s]=newspinmap[s];
				Job.animate();
			}	
		}
	}
	
	public void TCriticalTquench(IsingStructure Ising, double Ti, double Tf, double steplimit, int copies)
	{
		int limit=(int)steplimit;
		Ising.dilutionmap(Ising.R);
		for(int a=0; a<M; a++)
		{
			spintemp[a]=0;
			spintotal[a]=0;
		}
		params.set("H",0);
		
		for(int c=0; c<copies; c++)
		{
			params.set("copies", c+1);
			TF=new TemperatureField(Ising);
			double mag=0;
			Random cflip= new Random(c);
			
			params.set("T",Ti); 
			TF.GeneratingT(Ti,Ising.dilutionmap);
			
			for(int heat=0; heat<10; heat++)
			{
				TF.TMCS(Ti, 0, cflip, 1);
			    params.set("MCS", heat-10);
	     	    params.set("magnetization", TF.PIS.Magnetization());
		        Job.animate();
			}
			
			params.set("T",Tf); 
			TF.GeneratingT(Tf, Ising.dilutionmap);
			
			for(int step=0; step<limit; step++)
			{
				TF.TMCS(Tf, 0, cflip, 1);
				Job.animate();
				mag=TF.PIS.Magnetization();
				params.set("MCS",step);
				params.set("magnetization",mag);
			}
			TF.TMCS(Tf, 0, cflip, steplimit-limit);
			for(int j=0; j<M; j++)
			{
				spintemp[j]=TF.PIS.spin[j];
			}
			
			
			
			for(int astep=0; (astep<100)&(Math.abs(mag)<threshold); astep++)
			{

				TF.TMCS(Tf, 0, cflip, 1);
				Job.animate();
				mag=TF.PIS.Magnetization();
				params.set("MCS",astep+limit);
				params.set("magnetization",mag);
			}
			
			if(mag>0)// +1 is the stable direction
			{
				for(int l=0; l<M; l++)
					spintotal[l]+=spintemp[l];
			}
			if(mag<0)  //-1 is the stable direction
			{
				for(int n=0; n<M; n++)
					spintotal[n]+=(-spintemp[n]);
			}
			
			if(c==copies-10)
			{
				newspinmap=takeoutzeroes(Initialising, spintotal);
				for(int s=0; s<M; s++)
					spinmap[s]=newspinmap[s];
				Job.animate();
			}	
		}
	}

	
	public void OffCriticalTquench(IsingStructure Ising, double Ti, double Tf, double h,double steplimit, int copies)
	{
		int limit=(int)steplimit;
		
		for(int a=0; a<M; a++)
		{
			spintemp[a]=0;
			spintotal[a]=0;
		}
		params.set("H",h);
		
		for(int c=0; c<copies; c++)
		{
			params.set("copies", c+1);
			IS=Ising.clone();
			double mag=0;
			Random cflip= new Random(c);
			for(int heat=0; heat<10; heat++)
			{
				params.set("T",Ti);     
				IS.MCS(Ti, h, cflip, 1);
			    params.set("MCS", heat-10);
	     	    params.set("magnetization", IS.Magnetization());
		        Job.animate();
			}
			
			for(int step=0; step<limit; step++)
			{
				params.set("T",Tf);
				IS.MCS(Tf, h, cflip, 1);
				Job.animate();
				mag=IS.Magnetization();
				params.set("MCS",step);
				params.set("magnetization",mag);
			}
			IS.MCS(Tf, h, cflip, steplimit-limit);
			for(int j=0; j<M; j++)
			{
				spintemp[j]=IS.spin[j];
			}
			

				for(int l=0; l<M; l++)
					spintotal[l]+=spintemp[l];

			if(c==copies-10)
			{
				newspinmap=takeoutzeroes(Initialising, spintotal);
				for(int s=0; s<M; s++)
					spinmap[s]=newspinmap[s];
				Job.animate();
			}
			
			
			
		}
		
		

		
	}
	
	public double[] takeoutzeroes(IsingStructure Ising, double map[])
	{
		double total=0;
		double newmap[]= new double[Ising.M];
		for(int j=0; j<Ising.M; j++)
		{
			if(Ising.spin[j]!=0)
			total+=map[j];
		}
		double average=0;
		average= total/(Ising.M-Ising.deadsites);
		for(int x=0; x<Ising.M; x++)
		{
			newmap[x]=map[x];
			if(Ising.spin[x]==0)
				newmap[x]=average;
		}
		
		return newmap;
		
	}

	
	public void animate()
	{
		ColorPalette ising = new ColorPalette ();
		ising.setColor(1, Color.BLACK);      //up spin
		ising.setColor(-1, Color.WHITE);     //down spin
		ising.setColor(0, Color.RED);        //normal dilution
		ising.setColor(2, Color.BLUE);       //clusters
		ising.setColor(-2, Color.GREEN);     //
		
		ColorGradient heatmap = new ColorGradient();
		
		
		grid1.setColors(ising);
		grid1.registerData(IS.L1, IS.L2, IS.spin);
		grid2.setColors(heatmap);
		grid2.registerData(L1, L2, Initialising.dilutionmap);

		grid3.setColors(heatmap);
		grid3.registerData(L1, L2, spintotal);
		grid4.setColors(heatmap);
		grid4.registerData(L1, L2, spinmap);
		grid5.setColors(heatmap);
		grid5.registerData(L1, L2, TF.PIS.spin);
		grid6.setColors(heatmap);
		grid6.registerData(L1, L2, fluctuation);
	
	}

	public void clear()
	{
		grid1.clear();
		grid2.clear();
		grid3.clear();
		grid4.clear();
		grid5.clear();
	}
	
	public static void main (String[] dilutionmap){
		new Control(new dilutionmap(), "Kang Liu's dilutionmap" );
	}

	public void load(Control dilutionmap){
		/*dilutionmap.frame (grid1);
		dilutionmap.frame (grid2);
		dilutionmap.frame (grid3);
		dilutionmap.frame (grid4);
		dilutionmap.frame (grid5);*/
		dilutionmap.frameTogether("Display", grid1,grid2,grid3,grid4,grid5,grid6);
		params.add("L1", 100);
		params.add("L2", 100);
		params.add("R", 5);
		params.add("NJ",-4.0);	
		params.add("percent", new DoubleValue(0.05,0,1).withSlider());
		params.add("biaspercent", new DoubleValue(0.05,0,1).withSlider());
		params.add("deadsites");	

		
		params.addm("T", new DoubleValue(1.547, 0, 10).withSlider());
		params.addm("H", new DoubleValue(0, -2, 2).withSlider());
		params.add("threshold", 0.5);
		
		params.add("MCS");
		params.add("copies");
		params.add("magnetization");
	}
	
	public void run(){
		
		
		L1 = (int)params.fget("L1");
		L2 = (int)params.fget("L2");
		M = L1 * L2;
		R = (int)params.fget("R");
		NJ = params.fget("NJ");
		percent=params.fget("percent");
		biaspercent=params.fget("biaspercent");
		threshold= params.fget("threshold");
		
		isingspin= new int[M];
		spintotal= new double[M];
		spintemp=new int[M];
		spinmap= new double[M];
		newspinmap=new double[M];
		fluctuation=new double[M];
		
		
		Dseed = 1;
		Bseed = 1;
		Sseed = 1;
		
	    Initialising=new IsingStructure(L1,L2,R,NJ,percent,biaspercent);
	    Initialising.Dinitialization(Dseed, Bseed, R, R);
	    params.set("deadsites",Initialising.deadsites);
	    Initialising.Sinitialization(0, Sseed);
	    Initialising.dilutionmap(R);
	    IS= Initialising.clone();
	    TF= new TemperatureField(IS);
	    Job.animate();
	    FluctuationMap(Initialising, 4, 10000, 5);

	    
	    
	    double estep=0;
	    estep=1;
	    //TCriticalTquench(Initialising, 9, 0.5, estep, 30010);
	    
	    //OffCriticalTquench(Initialising, 9, 0.5, 0.1, estep, 100010);
	    
	    //movie(grid2, 0000,0000);
	    //movie(grid3, 3333,(int)(estep*1000));
	    //movie(grid4, 4444,(int)(estep*1000));
	    movie(grid6, 6666,6666);

	    

      
        
	    
	    
	    
	    
	    
	    

	}
	
	
	
	
	
}