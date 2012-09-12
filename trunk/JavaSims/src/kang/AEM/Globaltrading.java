package kang.AEM;

import java.awt.Color;
import java.text.DecimalFormat;


import kang.AEM.BasicStructure.AEMStructure;
import kang.util.PrintUtil;
import kang.ising.BasicStructure.BasicTools;
import kang.ising.BasicStructure.IsingStructure;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import chris.util.Random;


public class Globaltrading extends Simulation{
	
	Grid grid1=new Grid("wealth distribution");     // the map to display the wealth distribution	
	
	
	public AEMStructure AH;
	public AEMStructure AHtemp;


	
	public BasicTools Tools;
	
	private DecimalFormat fmt = new DecimalFormat("000");
	private DecimalFormat bmt = new DecimalFormat("0000");
	

	//initialization parameters
	public int L,R;
	public double M;


	public double percent;
	public double tax;
	public double growth;
	public double Ngrowth;
	public double mu;      //the exponential growth coeffcient

	public int time;
	public double order;
	

	
	public void animate()
	{
	
	    ColorGradient heatmap = new ColorGradient();
		
		
		grid1.setColors(heatmap);
		grid1.registerData(L, L, AHtemp.wealth);
		
		params.set("totalwealth", AHtemp.totalwealth);
		params.set("meanwealth", AHtemp.meanwealth);
		params.set("order", AHtemp.order);
		
	}
	
	public void clear()
	{
		grid1.clear();
	
	}
	
	public static void main (String[] Globaltrading){
		new Control(new Globaltrading(), "Kang Liu's global trading AEM" );
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
	
	public void singlerunfast(AEMStructure aem, int steplimit, int seed)
	{
		Random rand= new Random(seed);
		for(int t=0; t<steplimit; t++)
		{
			aem.TSfast(rand, aem.percent, aem.tax, aem.alpha, aem.Ngrowth);
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"fast wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"fast wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"fast wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"fast wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"fast wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}	
			if(t==100000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"fast wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}	
		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void singlerunexp(AEMStructure aem, int steplimit, int seed, double mu)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.TSfast(rand, aem.percent, aem.tax, aem.alpha, ngrowth);
			ngrowth=ngrowth*(1+mu);
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			
			
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"exp wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"exp wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"exp wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"exp wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"exp wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}	

		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void singlerun(AEMStructure aem, int steplimit, int seed)
	{
		Random rand= new Random(seed);
		for(int t=0; t<steplimit; t++)
		{
			aem.TS(rand, aem.percent, aem.tax, aem.alpha, aem.growth);
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}	
			if(t==100000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}	
		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void singletrajectory(AEMStructure aem, int steplimit, int seed, int[] index)
	{
		Random rand= new Random(seed);
		
		for(int t=0; t<steplimit; t++)
		{
			aem.TS(rand, aem.percent, aem.tax, aem.alpha, aem.growth);
			params.set("time", t+1);
			Job.animate();
			if(t%200==0)
			{
				for(int j=0; j<index.length; j++)
				{
					String saveas="/Users/liukang2002507/Desktop/simulation/AEM/trajectory/"+"singletrajectory <L="+fmt.format(aem.L1)+", j="+fmt.format(j)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
					PrintUtil.printlnToFile(saveas, t+1, aem.wealth[index[j]], aem.meanwealth);
				}
			}
			
		}
	
	}
	
	
	public double[] findrich(double data[], int number)   // the function to find the richest agents in the array
	{
		double[] richest=new double[number];
		double threshold=data[0];
		
		int doorman=0;     //the index for the poorest in the rich ones
		for(int i=0; i<number; i++)
		{
			richest[i]=data[i];
			if(data[i]<threshold)    
				{
				threshold=data[i];
				doorman=i;                  //find the doorman, the poorest in the first n agents
				}
		}
		
		for(int j=number; j<data.length; j++)
		{
			
			if(data[j]>threshold)
			{
				richest[doorman]=data[j];    //replace the doorman with the new member
				//now find the new doorman
				threshold=richest[0];
				doorman=0;
				for(int k=0; k<number; k++)
				{
					if(richest[k]<threshold)
						{
						threshold=richest[k];
						doorman=k;
						}
				}
			}
			
		}
		
		return richest;
	}
	
	public double[] findpoor(double data[], int number)   // the function to find the richest agents in the array
	{
		double[] poorest=new double[number];
		double threshold=data[0];
		
		int doorman=0;     //the index for the poorest in the rich ones
		for(int i=0; i<number; i++)
		{
			poorest[i]=data[i];
			if(data[i]>threshold)
				{
				threshold=data[i];
				doorman=i;                  //find the doorman
				}
		}
		for(int j=number; j<data.length; j++)
		{
			
			if(data[j]<threshold)
			{
				poorest[doorman]=data[j];    //replace the doorman with the new member
				
				
				//now find the new doorman
				threshold=poorest[0];
				doorman=0;
				for(int k=0; k<number; k++)
				{
					if(poorest[k]>threshold)
						{
						threshold=poorest[k];
						doorman=k;
						}
				}
			}
			
		}
		
		return poorest;
	}
	
	
	public void extremerun(AEMStructure aem, int steplimit, int number, int seed)    // keep track the richest and the poorest several agents total wealth over time
	{
		double poor[]=new double[number];
		double rich[]=new double[number];
		
		double poorP, richP;
		double totalpoor=0;
		double totalrich=0;
		double totalwealth=0;
		Random rand = new Random(seed);
		for(int t=0; t<steplimit; t++)
		{
			aem.TS(rand, aem.percent, aem.tax, aem.alpha, aem.Ngrowth);
			params.set("time", t+1);
			Job.animate();
			if(t%20==0)
			{
				totalpoor=0;
				totalrich=0;
				
				poor=findpoor(aem.wealth, number);
				rich=findrich(aem.wealth, number);
				for(int j=0; j<number; j++)
				{
					totalpoor+=poor[j];
					totalrich+=rich[j];
					
				}
				totalwealth=0;
				for(int k=0; k<aem.M; k++)
				{
				    totalwealth+=aem.wealth[k];
				}
				//totalwealth=aem.totalwealth;
				poorP=totalpoor/totalwealth;
				richP=totalrich/totalwealth;
				
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/trajectory/"+"extrem <L="+fmt.format(aem.L1)+", n="+fmt.format(number)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
				PrintUtil.printlnToFile(saveas, t+1, poorP, richP, totalpoor, totalrich, totalwealth);
			}
			
		}
	}
	
	
	
	public void extremerunexp(AEMStructure aem, int steplimit, int number, int seed, double mu)    // keep track the richest and the poorest several agents total wealth over time
	{
		double poor[]=new double[number];
		double rich[]=new double[number];
		
		double poorP, richP;
		double totalpoor=0;
		double totalrich=0;
		double totalwealth=0;
		Random rand = new Random(seed);
		double ngrowth=aem.Ngrowth;
		for(int t=0; t<steplimit; t++)
		{
			aem.TS(rand, aem.percent, aem.tax, aem.alpha, ngrowth);
			params.set("time", t+1);
			Job.animate();
			ngrowth=ngrowth*(1+mu);
			
			if(t%20==0)
			{
				totalpoor=0;
				totalrich=0;
				
				poor=findpoor(aem.wealth, number);
				rich=findrich(aem.wealth, number);
				for(int j=0; j<number; j++)
				{
					totalpoor+=poor[j];
					totalrich+=rich[j];
					
				}
				totalwealth=0;
				for(int k=0; k<aem.M; k++)
				{
				    totalwealth+=aem.wealth[k];
				}
				//totalwealth=aem.totalwealth;
				poorP=totalpoor/totalwealth;
				richP=totalrich/totalwealth;
				
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/trajectory/"+"extremexp <L="+fmt.format(aem.L1)+", n="+fmt.format(number)+", p= "+fmt.format(aem.percent*1000)+", mu= "+bmt.format(mu*10000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
				PrintUtil.printlnToFile(saveas, t+1, poorP, richP, totalpoor, totalrich, totalwealth);
			}
			
		}
	}
	
	public void output(AEMStructure aem, String path)
	{
		for(int i=0; i<aem.M; i++)
		{
			PrintUtil.printlnToFile(path, i+1, aem.wealth[i]);
		}
	}
	
	
	public void load(Control Globaltrading)
	{

		Globaltrading.frameTogether("Display", grid1);

		params.add("L", 50);
		params.add("R", 50);
		
	    params.add("tax",0.00);
		params.add("percent", 0.10);
		params.add("Ngrowth", 1.0);
		params.add("mu",0.0005);
		
		params.addm("time", 0);
		params.add("totalwealth");
		params.add("meanwealth");
		params.add("order");

	}
	
	public void run(){
		
		
		L = (int)params.fget("L");
		R =(int)params.fget("R");
		M = L * L;


		percent=params.fget("percent");
		tax=params.fget("tax");
	
		Ngrowth=params.fget("Ngrowth");
		mu=params.fget("mu");
		
	    AH=new AEMStructure(L,L,R,percent,tax, 0, Ngrowth);   
	    AHtemp=new AEMStructure(L,L,R,percent,tax, 0,Ngrowth);

	    
	    Tools=new BasicTools();

	    
	    {//initialization
	    	
	    	AH.Uinitialization(100);

	    	AHtemp=AH.clone();
	    	
	    }
	    
	    Job.animate();
	   
	    //singlerun(AHtemp, 101000,1);
	    //singlerunfast(AHtemp, 101000,1);
	    //singlerunexp(AHtemp, 51000, 1, mu);
	    
	    //extremerun(AHtemp, 20000, 25, 1);
	    extremerunexp(AHtemp, 20000, 25, 1, mu);
	    
	    /*{
	    	int[] index=new int[5];
	    	for(int j=0; j<index.length; j++)
	    	{
	    		index[j]=j;
	    	}
	    	
	    	singletrajectory(AHtemp, 100000, 1, index);
	    }*/
	    
	    Job.animate();

	}
	
}