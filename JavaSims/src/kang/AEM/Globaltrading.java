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
		
	public void entropyrunEXP(AEMStructure aem, int steplimit, int seed, double mu)
	{
		Random rand= new Random(seed);
		double ngrowth=aem.Ngrowth;
		for(int t=0; t<steplimit; t++)
		{
			aem.TSfast(rand, aem.percent, aem.tax, aem.alpha, ngrowth);
			ngrowth=ngrowth*(1+mu);
			String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"entropy <L="+fmt.format(aem.L1)+", p= "+fmt.format(aem.percent*1000)+",  mu="+fmt.format(mu*10000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
			
			params.set("time", t+1);
			Job.animate();
			if(t%100==0)
			{
				PrintUtil.printlnToFile(saveas, t+1, aem.order);
			}

		}
	}
	
	public void BiasentropyrunEXP(double biasp, AEMStructure aem, int steplimit, int seed, double mu)
	{
		Random rand= new Random(seed);
		double ngrowth=aem.Ngrowth;
		for(int t=0; t<steplimit; t++)
		{
			aem.BiasTSfast(biasp, rand, aem.percent, aem.tax, aem.alpha, ngrowth);
			ngrowth=ngrowth*(1+mu);
			String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Biasentropy <L="+fmt.format(aem.L1)+", biasp= "+fmt.format(biasp*1000)+", p= "+fmt.format(aem.percent*1000)+",  mu="+fmt.format(mu*10000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
			
			params.set("time", t+1);
			Job.animate();
			if(t%100==0)
			{
				PrintUtil.printlnToFile(saveas, t+1, aem.order);
			}

		}
	}
	
	public void FIXentropyrunEXP(AEMStructure aem, int steplimit, int seed, double mu, double fixamount)    //entropy run for the model with fixed trading amount
	{
		Random rand= new Random(seed);
		double ngrowth=aem.Ngrowth;
		for(int t=0; t<steplimit; t++)
		{
			aem.TSfastFix(rand, fixamount, aem.tax, aem.alpha, ngrowth);
			ngrowth=ngrowth*(1+mu);
			String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"FIXentropy <L="+fmt.format(aem.L1)+", W= "+fmt.format(fixamount)+",  mu="+fmt.format(mu*10000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
			
			params.set("time", t+1);
			Job.animate();
			if(t%100==0)
			{
				PrintUtil.printlnToFile(saveas, t+1, aem.order);
			}

		}
	}
	
	public void PFIXentropyrunEXP(AEMStructure aem, int steplimit, int seed, double mu, double fixamountpercent)    //entropy run for the model with fixed trading amount
	{
		Random rand= new Random(seed);
		double ngrowth=aem.Ngrowth;
		for(int t=0; t<steplimit; t++)
		{
			aem.TSfastFixExp(rand, fixamountpercent, aem.tax, aem.alpha, ngrowth);
			ngrowth=ngrowth*(1+mu);
			String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"PFIXentropy <L="+fmt.format(aem.L1)+", fixp= "+fmt.format(fixamountpercent*10000)+",  mu="+fmt.format(mu*10000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
			
			params.set("time", t+1);
			Job.animate();
			if(t%100==0)
			{
				PrintUtil.printlnToFile(saveas, t+1, aem.order);
			}

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

	public void negsubexp(AEMStructure aem, int steplimit, int seed, double mu, double gamma)  //negative sub-linear growth favoring the poor 
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"neg check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		String entropy="/Users/liukang2002507/Desktop/simulation/AEM/"+"neg order <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		
		double entropyC=Math.log(aem.M);
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, -gamma);
			ngrowth=ngrowth*(1+mu);
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth, aem.sumtrading, aem.sumflow);
			if(t%50==0)
			{
				PrintUtil.printlnToFile(entropy, t+1, aem.order, aem.order/M+entropyC);
			}
			
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"neg wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"neg wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"neg wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"neg wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"neg wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}	
			if(t==100000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"neg wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}

		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}	
	
	public void CompareSeed(AEMStructure aem, int steplimit, int totalseed, double mu, double gamma)  //the function to test the effect of different random number seeds on the wealth distribution
	{
		for(int seed=1; seed<=totalseed; seed++)
		{
			AHtemp=aem.clone();
			double ngrowth=aem.Ngrowth;
			Random rand= new Random(seed);
			//String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"sub check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
			//String entropy="/Users/liukang2002507/Desktop/simulation/AEM/"+"sub order <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
            
			
			double entropyC=Math.log(aem.M);
			
			
			for(int t=0; t<steplimit; t++)
			{
				
				AHtemp.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
				ngrowth=ngrowth*(1+mu);
				//PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
				/*if(t%50==0)
				{
					PrintUtil.printlnToFile(entropy, t+1, aem.order, aem.order/M+entropyC);
				}*/
				
				params.set("time", t+1);
				Job.animate();
	
				if(t==10000)
				{
					String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Compare Seed <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", seed= "+fmt.format(seed)+">.txt";
			        output(AHtemp, saveas);
				}

			}
		}
	}
	
	public void sublinearexp(AEMStructure aem, int steplimit, int seed, double mu, double gamma)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"sub check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		String entropy="/Users/liukang2002507/Desktop/simulation/AEM/"+"sub order <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";

		double entropyC=Math.log(aem.M);
		
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			if(t%50==0)
			{
				PrintUtil.printlnToFile(entropy, t+1, aem.order, aem.order/M+entropyC);
			}
			
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"sub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"sub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"sub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"sub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"sub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}	
			if(t==100000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"sub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}

		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void Biassublinearexp(AEMStructure aem, int steplimit, int seed, double mu, double gamma, double biasp)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"biassub check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", biasp= "+fmt.format(biasp*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		String entropy="/Users/liukang2002507/Desktop/simulation/AEM/"+"biassub order <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", biasp= "+fmt.format(biasp*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";

		double entropyC=Math.log(aem.M);
		
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.BiasSubTS(rand, biasp, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			if(t%50==0)
			{
				PrintUtil.printlnToFile(entropy, t+1, aem.order, aem.order/M+entropyC);
			}
			
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"biassub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", biasp= "+fmt.format(biasp*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"biassub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", biasp= "+fmt.format(biasp*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"biassub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", biasp= "+fmt.format(biasp*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"biassub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", biasp= "+fmt.format(biasp*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"biassub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", biasp= "+fmt.format(biasp*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}	
			if(t==100000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"biassub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", biasp= "+fmt.format(biasp*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}

		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void IndexSkill(AEMStructure aem, int steplimit, int seed, double mu, double gamma)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"IndexSkill check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		String entropy="/Users/liukang2002507/Desktop/simulation/AEM/"+"IndexSkill order <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";

		double entropyC=Math.log(aem.M);
		
		//generate the skill matrix
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.IndexSkillTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			if(t%50==0)
			{
				PrintUtil.printlnToFile(entropy, t+1, aem.order, aem.order/M+entropyC);
			}
			
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"IndexSkill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"IndexSkill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"IndexSkill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"IndexSkill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"IndexSkill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}	
			if(t==100000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"IndexSkill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}

		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void Volsubexp(AEMStructure aem, int steplimit, int seed, double mu, double gamma, double sigma)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"Volsub check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", sigma= "+fmt.format(sigma*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		String entropy="/Users/liukang2002507/Desktop/simulation/AEM/"+"Volsub order <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", sigma= "+fmt.format(sigma*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";

		double entropyC=Math.log(aem.M);
		
		//generate the skill matrix

		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.VolSubTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma,sigma);
			ngrowth=ngrowth*(1+mu);
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			if(t%50==0)
			{
				PrintUtil.printlnToFile(entropy, t+1, aem.order, aem.order/M+entropyC);
			}
			
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Volsub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", sigma= "+fmt.format(sigma*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Volsub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", sigma= "+fmt.format(sigma*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Volsub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", sigma= "+fmt.format(sigma*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Volsub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", sigma= "+fmt.format(sigma*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Volsub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", sigma= "+fmt.format(sigma*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}	
			if(t==100000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Volsub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", sigma= "+fmt.format(sigma*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}

		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void Skillsublinearexp(AEMStructure aem, int steplimit, int seed, double mu, double gamma, double kappa)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"skillsub check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		String entropy="/Users/liukang2002507/Desktop/simulation/AEM/"+"skillsub order <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";

		double entropyC=Math.log(aem.M);
		
		//generate the skill matrix
		for(int j=0; j<aem.M; j++)
		{
			aem.skill[j]=Math.pow(j+1, kappa);
		}
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.SkillSubTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			if(t%50==0)
			{
				PrintUtil.printlnToFile(entropy, t+1, aem.order, aem.order/M+entropyC);
			}
			
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"skillsub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"skillsub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"skillsub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"skillsub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"skillsub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}	
			if(t==100000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"skillsub wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}

		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void EXPSkillsublinearexp(AEMStructure aem, int steplimit, int seed, double mu, double gamma, double kappa)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"EXPskill check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		String entropy="/Users/liukang2002507/Desktop/simulation/AEM/"+"EXPskill order <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";

		double entropyC=Math.log(aem.M);
		
		//generate the skill matrix
		
		aem.skill[0]=1;
		for(int j=1; j<aem.M; j++)
		{
			aem.skill[j]=aem.skill[j-1]*kappa;
		}
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.SkillSubTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			if(t%50==0)
			{
				PrintUtil.printlnToFile(entropy, t+1, aem.order, aem.order/M+entropyC);
			}
			
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"EXPskill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"EXPskill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"EXPskill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"EXPskill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"EXPskill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}	
			if(t==100000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"EXPskill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}

		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void GaussianSkill(AEMStructure aem, int steplimit, int seed, double mu, double gamma, double kappa)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		Random grand= new Random(seed+1);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"Gskill check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		String entropy="/Users/liukang2002507/Desktop/simulation/AEM/"+"Gskill order <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";

		double entropyC=Math.log(aem.M);
		
		//generate the skill matrix
		
		
		for(int j=0; j<aem.M; j++)
		{
			aem.skill[j]=1+kappa*grand.nextGaussian();
			if(aem.skill[j]<0)
				aem.skill[j]=0;
		}
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.SkillSubTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			if(t%50==0)
			{
				PrintUtil.printlnToFile(entropy, t+1, aem.order, aem.order/M+entropyC);
			}
			
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Gskill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        skilloutput(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Gskill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        skilloutput(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Gskill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        skilloutput(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Gskill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        skilloutput(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Gskill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        skilloutput(aem, saveas);
			}	
			if(t==100000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Gskill wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        skilloutput(aem, saveas);
			}

		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void CorrelationRun(AEMStructure aem, int steplimit, int seed, double mu, double gamma, int stepsize)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String path="/Users/liukang2002507/Desktop/simulation/AEM/"+"correlation <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", step= "+fmt.format(stepsize)+">.txt";
		
		double[] wi= new double[aem.M];
		double[] wf= new double[aem.M];
		double correlation=0;
		
		
		for(int j=0; j<aem.M; j++)
		{
			wi[j]=aem.wealth[j];
		}
		
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			
			if(t%stepsize==0)
			{
				
				for(int jj=0; jj<aem.M; jj++)
				{
					wf[jj]=aem.wealth[jj];
					
				}
				correlation=Tools.Correlation(wi, wf, aem.M);
				PrintUtil.printlnToFile(path, t+1, correlation, aem.totalwealth, aem.sumtrading, aem.sumflow);
				for(int jf=0; jf<aem.M; jf++)
				{
					wi[jf]=aem.wealth[jf];           //reset the initial wealth distribution
				}
				
			}
			
			params.set("time", t+1);
			Job.animate();

            

		}
		
	}
	
	public void SkillCorrelationRun(AEMStructure aem, int steplimit, int seed, double mu, double gamma, int stepsize, double kappa)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String path="/Users/liukang2002507/Desktop/simulation/AEM/"+"skillcorrelation <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*100)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", step= "+fmt.format(stepsize)+">.txt";
		
		double[] wi= new double[aem.M];
		double[] wf= new double[aem.M];
		double correlation=0;
		double skillC=0;  // the correlation with the skill distribution
		
		//generate the skill matrix
		for(int j=0; j<aem.M; j++)
		{
			aem.skill[j]=Math.pow(j+1, kappa);
		}
		
		
		for(int j=0; j<aem.M; j++)
		{
			wi[j]=aem.wealth[j];
		}
		
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.SkillSubTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			
			if(t%stepsize==0)
			{
				
				for(int jj=0; jj<aem.M; jj++)
				{
					wf[jj]=aem.wealth[jj];
					
				}
				correlation=Tools.Correlation(wi, wf, aem.M);
				skillC=Tools.Correlation(aem.skill, wf, aem.M);
				PrintUtil.printlnToFile(path, t+1, correlation, aem.totalwealth, aem.sumtrading, aem.sumflow, skillC);
				for(int jf=0; jf<aem.M; jf++)
				{
					wi[jf]=aem.wealth[jf];           //reset the initial wealth distribution
				}
				
			}
			
			params.set("time", t+1);
			Job.animate();

            

		}
		
	}
	
	public void EXPSkillCorrelationRun(AEMStructure aem, int steplimit, int seed, double mu, double gamma, int stepsize, double kappa)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String path="/Users/liukang2002507/Desktop/simulation/AEM/"+"EXPskillcorrelation <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", kappa= "+fmt.format(kappa*10000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", step= "+fmt.format(stepsize)+">.txt";
		
		double[] wi= new double[aem.M];
		double[] wf= new double[aem.M];
		double correlation=0;
		double skillC=0;  // the correlation with the skill distribution
		
		//generate the skill matrix
		aem.skill[0]=1;
		for(int j=1; j<aem.M; j++)
		{
			aem.skill[j]=aem.skill[j-1]*kappa;
		}
		
		
		for(int j=0; j<aem.M; j++)
		{
			wi[j]=aem.wealth[j];
		}
		
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.SkillSubTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			
			if(t%stepsize==0)
			{
				
				for(int jj=0; jj<aem.M; jj++)
				{
					wf[jj]=aem.wealth[jj];
					
				}
				correlation=Tools.Correlation(wi, wf, aem.M);
				skillC=Tools.Correlation(aem.skill, wf, aem.M);
				PrintUtil.printlnToFile(path, t+1, correlation, aem.totalwealth, aem.sumtrading, aem.sumflow, skillC);
				for(int jf=0; jf<aem.M; jf++)
				{
					wi[jf]=aem.wealth[jf];           //reset the initial wealth distribution
				}
				
			}
			
			params.set("time", t+1);
			Job.animate();

            

		}
		
	}
	
	public void BiasCorrelationRun(AEMStructure aem, int steplimit, int seed, double mu, double gamma, int stepsize, double biasp)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String path="/Users/liukang2002507/Desktop/simulation/AEM/"+"biascorrelation <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", biasp= "+fmt.format(biasp*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", step= "+fmt.format(stepsize)+">.txt";
		
		double[] wi= new double[aem.M];
		double[] wf= new double[aem.M];
		double correlation=0;
		
		
		for(int j=0; j<aem.M; j++)
		{
			wi[j]=aem.wealth[j];
		}
		
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.BiasSubTS(rand, biasp, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			
			if(t%stepsize==0)
			{
				
				for(int jj=0; jj<aem.M; jj++)
				{
					wf[jj]=aem.wealth[jj];
					
				}
				correlation=Tools.Correlation(wi, wf, aem.M);
				PrintUtil.printlnToFile(path, t+1, correlation, aem.totalwealth, aem.sumtrading, aem.sumflow);
				for(int jf=0; jf<aem.M; jf++)
				{
					wi[jf]=aem.wealth[jf];           //reset the initial wealth distribution
				}
				
			}
			
			params.set("time", t+1);
			Job.animate();

            

		}
		
	}
	
	public void CAPexp(AEMStructure aem, int steplimit, int seed, int richpeople, double mu, double gamma)
	{
		
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"CAP check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", rich= "+bmt.format(richpeople)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		String entropy="/Users/liukang2002507/Desktop/simulation/AEM/"+"CAP order <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", rich= "+bmt.format(richpeople)+", gamma= "+fmt.format(gamma*1000)+">.txt";

		double entropyC=Math.log(aem.M);
		
		
		for(int t=0; t<steplimit; t++)
		{
			
			
			aem.CapTS(rand, aem.percent, aem.tax,  aem.alpha, richpeople, mu, gamma);
			
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth, aem.sumtrading, aem.sumflow);
			if(t%50==0)
			{
				PrintUtil.printlnToFile(entropy, t+1, aem.order, aem.order/M+entropyC);
			}
			
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"CAP wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", rich= "+bmt.format(richpeople)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==1000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"CAP wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", rich= "+bmt.format(richpeople)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"CAP wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", rich= "+bmt.format(richpeople)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"CAP wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", rich= "+bmt.format(richpeople)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"CAP wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", rich= "+bmt.format(richpeople)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"CAP wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", rich= "+bmt.format(richpeople)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}	
			if(t==100000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"CAP wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", rich= "+bmt.format(richpeople)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
			}

		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void singlerunexp(AEMStructure aem, int steplimit, int seed, double mu)   // geometric growth with uniform distribution
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
			if(t==100000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"exp wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}

			
		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void Unfairsinglerunexp(AEMStructure aem, int steplimit, int seed, double pmu)    //growth linearly depending on the wealth
	{
		
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"Unfair check <L="+fmt.format(aem.L1)+", pmu="+fmt.format(pmu*1000000)+", p= "+fmt.format(aem.percent*1000)+">.txt";
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.UnfairTSfast(rand, aem.percent, aem.tax, aem.alpha, pmu);
		
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			
			
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Unfair wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(pmu*1000000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Unfair wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(pmu*1000000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Unfair wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(pmu*1000000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Unfair wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(pmu*1000000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Unfair wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(pmu*1000000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}	

		}
		
		
		
	}
	
	public void FIXsinglerunexp(AEMStructure aem, int steplimit, int seed, double mu, double fixamount)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"FIX check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", W= "+fmt.format(fixamount)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.TSfastFix(rand, fixamount, aem.tax, aem.alpha, ngrowth);
			ngrowth=ngrowth*(1+mu);
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			
			
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"FIX wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", W= "+fmt.format(fixamount)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"FIX wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", W= "+fmt.format(fixamount)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"FIX wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", W= "+fmt.format(fixamount)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"FIX wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", W= "+fmt.format(fixamount)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"FIX wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", W= "+fmt.format(fixamount)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}	

		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void BiasFIXsinglerunexp(double biasp, AEMStructure aem, int steplimit, int seed, double mu, double fixamount)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"BiasFIX check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", biasp= "+fmt.format(biasp*1000)+", W= "+fmt.format(fixamount)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.BiasTSfastFix(biasp, rand, fixamount, aem.tax, aem.alpha, ngrowth);
			ngrowth=ngrowth*(1+mu);
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			
			
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"BiasFIX wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", biasp= "+fmt.format(biasp*1000)+", W= "+fmt.format(fixamount)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"BiasFIX wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", biasp= "+fmt.format(biasp*1000)+", W= "+fmt.format(fixamount)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"BiasFIX wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", biasp= "+fmt.format(biasp*1000)+", W= "+fmt.format(fixamount)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"BiasFIX wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", biasp= "+fmt.format(biasp*1000)+", W= "+fmt.format(fixamount)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"BiasFIX wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", biasp= "+fmt.format(biasp*1000)+", W= "+fmt.format(fixamount)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}	

		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void Biassinglerunexp(double biasp, AEMStructure aem, int steplimit, int seed, double mu)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"biascheck <L="+fmt.format(aem.L1)+", biasp= "+fmt.format(biasp*1000)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.BiasTSfast(biasp, rand, aem.percent, aem.tax, aem.alpha, ngrowth);
			ngrowth=ngrowth*(1+mu);
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			
			
			params.set("time", t+1);
			Job.animate();
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"biasexp wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"biasexp wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"biasexp wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"biasexp wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
		        output(aem, saveas);
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"biasexp wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+">.txt";
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
	
	public void ScanP(AEMStructure aem, int steplimit, double minP, double maxP, double dP, double gamma, Boolean oneseed, int number)
	{
		
		//String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"ScanP <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(0)+", gamma= "+fmt.format(gamma*1000)+", n= "+bmt.format(number)+">.txt";
		
		String saveas;
		if(oneseed)
		{
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanP data/1seed/"+"ScanP <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(0)+", gamma= "+fmt.format(gamma*1000)+", n= "+bmt.format(number)+">.txt";
		}
		else
		{
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanP data/nseed/"+"ScanP <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(0)+", gamma= "+fmt.format(gamma*1000)+", n= "+bmt.format(number)+">.txt";
		}
		
		
		Random rand;
		int seed=1;
		
		for(double percent=minP; percent<maxP; percent+=dP)
		{
			params.set("percent", percent);
			AHtemp=new AEMStructure(aem.L1,aem.L2,R,percent, tax, 0, aem.Ngrowth);
			

				AHtemp.Uinitialization(1);
				params.set("mu", 0);
				double ngrowth=0;
				rand=new Random(seed);
				if(oneseed)
				{
					
				}
				else
				{
					seed++;
				}
				
				
				Job.animate();
				
				for(int t=0; t<steplimit; t++)
				{
					
					AHtemp.SublinearTS(rand, percent, AHtemp.tax, AHtemp.alpha, ngrowth, gamma);
					ngrowth=0;
					params.set("time", t+1);
					Job.animate();
				
				}
				
				double poordata[]=new double[number];
				poordata=findpoor(AHtemp.wealth, number);
				double poorest=Tools.Mean(poordata, number);
				//double poorest=poorest(AHtemp.wealth);
				PrintUtil.printlnToFile(saveas, percent, poorest);
				
			
		}
	}
	
	public void ScanMuP(AEMStructure aem, int steplimit, double minMu, double maxMu, double dMu, double gamma, Boolean oneseed, int number)
	{
		double minP=0.01;
		double maxP=0.40;
		double dP=0.01;
		
		
		String saveas;
		if(oneseed)
		{
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanMuP data/1seed/"+"Scan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", n= "+bmt.format(number)+">.txt";
		}
		else
		{
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanMuP data/nseed/"+"Scan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", n= "+bmt.format(number)+">.txt";
		}
		
		int randseed=1;
		
		for(double percent=minP; percent<maxP; percent+=dP)
		{
			params.set("percent", percent);
			AHtemp=new AEMStructure(aem.L1,aem.L2,R,percent, tax, 0, aem.Ngrowth);
			
			Random rand;
			
			for(double mutemp=minMu; mutemp<maxMu; mutemp+=dMu)
			{
				
				AHtemp.Uinitialization(1);
				params.set("mu", mutemp);
				double ngrowth=AHtemp.Ngrowth;
				
				if(oneseed)
				{
					rand=new Random(1);
				}
				else
				{
					rand=new Random(randseed);
					randseed++;
				}
		

				Job.animate();
				
				for(int t=0; t<steplimit; t++)
				{
					
					AHtemp.SublinearTS(rand, percent, AHtemp.tax, AHtemp.alpha, ngrowth, gamma);
					ngrowth=ngrowth*(1+mutemp);
					params.set("time", t+1);
					Job.animate();
				
				}
				
				double poordata[]=new double[number];
				poordata=findpoor(AHtemp.wealth, number);
				double poorest=Tools.Mean(poordata, number);
				//double poorest=poorest(AHtemp.wealth);
				PrintUtil.printlnToFile(saveas, percent, mutemp, poorest);
				
			}
		}

	}
	
	public double poorest(double data[])
	{
		double poorest=data[0];
		for(int i=0; i<data.length; i++)
		{
			if(data[i]<=poorest)
				poorest=data[i];
		}
		return poorest;
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
	
	public double[] findpoor(double data[], int number)   // the function to find the poorest agents in the array
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
	
	public void skilloutput(AEMStructure aem, String path)
	{
		for(int i=0; i<aem.M; i++)
		{
			PrintUtil.printlnToFile(path, i+1, aem.wealth[i], aem.skill[i]);
		}
	}
	
	
	public void load(Control Globaltrading)
	{

		Globaltrading.frameTogether("Display", grid1);

		params.add("L", 50);
		params.add("R", 50);
		
	    params.add("tax",0.00);
		params.add("percent", 0.01);
		params.add("Ngrowth", 1.0);
		params.add("mu",0.0010);     //for mu=0.0001    for pmu=0.000001
		
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
	    	
	    	AH.Uinitialization(1);

	    	AHtemp=AH.clone();
	    	
	    }
	    
	    Job.animate();
	   
	    //singlerun(AHtemp, 101000,1);
	    //singlerunfast(AHtemp, 101000,1);
	    
	    //singlerunexp(AHtemp, 101000, 1, mu);
	    
	    //entropyrunEXP(AHtemp, 101000, 1, mu);
	    //FIXentropyrunEXP(AHtemp, 101000, 1, mu, 5);
	    //PFIXentropyrunEXP(AHtemp, 101000, 1, mu, 0.0005);
	    //BiasentropyrunEXP(0.6, AHtemp, 101000, 1, mu);
	    //Biassinglerunexp(0.6, AHtemp, 51000, 1, mu);
	    
	    
	    //FIXsinglerunexp(AHtemp, 51000, 1, mu, 10);
	    //BiasFIXsinglerunexp(0.6, AHtemp, 51000, 1, mu, 10);
	    //Unfairsinglerunexp(AHtemp, 101000, 1, mu);
	    
	    
	    //sublinearexp(AHtemp, 101000, 9, mu, 1.100);
	    
	    //CompareSeed(AH, 11000, 5, mu, 0.500);
	    
	    
	    //Biassublinearexp(AHtemp, 101000, 1, mu, 1.000, 0.50);
	    //Skillsublinearexp(AHtemp, 101000, 1, mu, 0.000, 10.00);
	    
	    //EXPSkillsublinearexp(AHtemp, 101000, 1, mu, 0.000, 1.1000);
	    //IndexSkill(AHtemp, 101000, 1, mu, 0.000);
	    //Volsubexp(AHtemp, 101000, 1, mu, 0.000, 0.00);
	    
	    //GaussianSkill(AHtemp, 101000, 1, mu, 0.000, 2.0000);
	    
	    //CorrelationRun(AHtemp, 101000, 1, mu, 0.000, 100);
	    //BiasCorrelationRun(AHtemp, 101000, 1, mu, 0.000, 100, 0.55);
	    //SkillCorrelationRun(AHtemp, 101000, 1, mu, 0.000, 100, 2.00);
	    
	    
	    //CAPexp(AHtemp, 101000, 1, 2500, mu, 0.000);
	    
	    //negsubexp(AHtemp, 101000, 1, mu, 9.000);
	    
	    
	    ScanMuP(AH, 2000, 0.0000, 0.4000, 0.0100, 1.000, true, 25);
	    
	    //ScanP(AH, 20000, 0.02, 0.40, 0.02, 1, false, 1);
	    
	    //extremerun(AHtemp, 10000, 1, 1);
	    //extremerunexp(AHtemp, 50000, 25, 1, mu);
	    
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