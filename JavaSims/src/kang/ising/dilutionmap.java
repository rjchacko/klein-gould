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
import chris.util.PrintUtil;


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
    public double correlation;
    public double errorC;
    public double Cor[];    // the array for the correlation value
    
    public int spintemp[];
    public double spintotal[];
    public double spinmap[];
    public double newspinmap[];
    
    public double RSDmap[];
    public double RSSmap[];
    
    
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
	
	public void FluctuationMap(IsingStructure Ising, double T, double H, int steplimit, int copies)
	{
		params.set("H",H);
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
			for(int heat=0; heat<11; heat++)
			{
				params.set("T",9);     
				IS.MCS(9, H, cflip, 1);
			    params.set("MCS", heat-10);
	     	    params.set("magnetization", IS.Magnetization());
		        Job.animate();
			}
			for(int prestep=0; prestep<2000; prestep++)
			{
				IS.MCS(T, H, cflip, 1);
				params.set("MCS",prestep-2000);
				params.set("magnetization",IS.Magnetization());
				Job.animate();
			}
			
			for(int step=0; step<steplimit; step++)
			{
				params.set("T",T);
				IS.MCS(T, H, cflip, 1);
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
			IS.MCS(T, H, flip, 1);
			Job.animate();
		}
		
		
		
	}
	
	public void CriticalTquench(IsingStructure Ising, double Ti, double Tf, double steplimit, int copies, int seed)
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
			Random cflip= new Random(c+100*seed);
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
			
			for(int astep=0; (astep<50)&(Math.abs(mag)<threshold); astep++)
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
				//correlation=Correlation(Initialising,Initialising.dilutionmap,spinmap);
				//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/dilutionmap/correlation.txt", Tf , 0 , correlation, steplimit, copies);
				
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

	
	public void OffCriticalTquench(IsingStructure Ising, double Ti, double Tf, double h,double steplimit, int copies, int seed)
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
			Random cflip= new Random(c+100*seed);
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
				//correlation=Correlation(Initialising,Initialising.dilutionmap,spinmap);
				//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/dilutionmap/correlation.txt", Tf , h , correlation, steplimit, copies);
			}
			
			
			
		}
		
		

		
	}
	
	public void SmallfieldTquench(IsingStructure Ising, double Ti, double Tf, double h,double steplimit, int copies, int seed)
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
			Random cflip= new Random(c+100*seed);
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
			
			for(int astep=0; (astep<50)&(Math.abs(mag)<threshold); astep++)
			{
				params.set("T",Tf);
				IS.MCS(Tf, h, cflip, 1);
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
				//correlation=Correlation(Initialising,Initialising.dilutionmap,spinmap);
				//PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/dilutionmap/correlation.txt", Tf , h , correlation, steplimit, copies);
				
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

	public double Correlation(IsingStructure Ising, double map1[], double map2[])
	{
		double c=0;
		double total1=0;
		double total2=0;
		double totalc=0;
		double fl1=0;
		double fl2=0;
		for(int j=0; j<Ising.M; j++)
		{
			if(Ising.spin[j]!=0)
			{
				total1+=map1[j];
				total2+=map2[j];
			}
		}
		double avg1=total1/(Ising.M-Ising.deadsites);
		double avg2=total2/(Ising.M-Ising.deadsites);
		for(int i=0; i<Ising.M; i++)
		{
			if(Ising.spin[i]!=0)
			{
				totalc+=(map1[i]-avg1)*(map2[i]-avg2);
				fl1+=(map1[i]-avg1)*(map1[i]-avg1);
				fl2+=(map2[i]-avg2)*(map2[i]-avg2);
			}
		}
		c=totalc/(Math.sqrt(fl1*fl2));
		
		return c;
	}
	
	public double average(double data[], int size)
	{
		double total=0;
		for(int j=0;j<size;j++)
			total+=data[j];
		return total/size;
	}
	
	public double errorbar(double data[], int size)
	{
		double error=0;
		double total=0;
		double total2=0;
		for(int j=0; j<size; j++)
		{
			total+=data[j];
			total2+=((data[j])*(data[j]));
		}
		error=Math.sqrt((total2/size)-(total/size)*(total/size));
		return error;
		
	}
	
	public double MIN(double data[], int size)
	{
		double min=data[0];
		for(int j=0; j<size; j++)
		{
			if(data[j]<=min)
				min=data[j];
		}
		
		return min;
	}
	
	public double MAX(double data[], int size)
	{
		double max=data[0];
		for(int j=0; j<size; j++)
		{
			if(data[j]>=max)
				max=data[j];
		}
		
		return max;
	}
	
	public double[] Rescale(double data[], int size)  // the function to rescale a double map data[] to be within range 0~1
	{
		double newdata[]= new double [size];
		double min=MIN(data, size);
		double max=MAX(data, size);
		for(int j=0; j<size; j++)
		{
			newdata[j]=(data[j]-min)/(max-min);
		}
		return newdata;
	}
	
	public void CalculateCor(IsingStructure Ising, double MinH, double MaxH, double dH, double smallH, int runs, double Ti, double Tf, double steplimit, int copies)
	{
		Cor=new double [runs];

		for(double h=MinH; h<=MaxH; h+=dH)
		{
			for(int z=0; z<runs; z++)
				Cor[z]=0;
			
			if(h==0)
			{
				for(int rr=0; rr<runs; rr++)
				{
					params.set("runs", rr+1);
					CriticalTquench(Ising, Ti, Tf, steplimit, copies, rr+1);
					RSSmap=Rescale(spinmap,Ising.M);
					Cor[rr]=Correlation(Initialising,RSDmap,RSSmap);	
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/dilutionmap/check.txt", Tf , h , Cor[rr], copies, rr);
				}
			}
			else if(h<=smallH)
			{
				for(int rr=0; rr<runs; rr++)
				{
					params.set("runs", rr+1);
					SmallfieldTquench(Ising, Ti, Tf,h, steplimit, copies, rr+1);
					RSSmap=Rescale(spinmap,Ising.M);
					Cor[rr]=Correlation(Initialising,RSDmap,RSSmap);
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/dilutionmap/check.txt", Tf , h , Cor[rr], copies, rr);
				}
			}
			if(h>smallH)
			{
				for(int rr=0; rr<runs; rr++)
				{
					params.set("runs", rr+1);
					OffCriticalTquench(Ising, Ti, Tf, h, steplimit, copies, rr+1);
					RSSmap=Rescale(spinmap,Ising.M);
					Cor[rr]=Correlation(Initialising,RSDmap,RSSmap);
					PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/dilutionmap/check.txt", Tf , h , Cor[rr], copies, rr);
				}
			}
		correlation=average(Cor,runs);	
		errorC=errorbar(Cor,runs);
		PrintUtil.printlnToFile("/Users/liukang2002507/Desktop/simulation/dilutionmap/correlation.txt", Tf , h , correlation, errorC, steplimit, copies, runs);
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
		grid6.clear();
		
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
		params.add("runs");
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
		RSDmap= new double [M];   //rescaled dilution map
		RSSmap= new double [M];   //rescaled spin map
		
		
		Dseed = 1;
		Bseed = 1;
		Sseed = 1;
		
	    Initialising=new IsingStructure(L1,L2,R,NJ,percent,biaspercent);
	    Initialising.Dinitialization(Dseed, Bseed, R, R);
	    params.set("deadsites",Initialising.deadsites);
	    Initialising.Sinitialization(0, Sseed);
	    Initialising.dilutionmap(R);
	    RSDmap=Rescale(Initialising.dilutionmap, M);
	    IS= Initialising.clone();
	    TF= new TemperatureField(IS);
	    Job.animate();
	    
	    
	    //FluctuationMap(Initialising, 2, 0, 10000, 1);

	    
	    
	    double estep=0;
	    estep=1;

	    CalculateCor(Initialising, 0, 0.3, 0.002, 0.15, 10, 9, 0.5, estep, 10010);
	    
	    //movie(grid2, 0000,0000);
	    //movie(grid3, 3333,(int)(estep*1000));
	    //movie(grid4, 4444,(int)(estep*1000));
	    //movie(grid6, 6666,6666);

	    

      
        
	    
	    
	    
	    
	    
	    

	}
	
	
	
	
	
}