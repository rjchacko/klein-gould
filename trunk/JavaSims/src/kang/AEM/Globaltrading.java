package kang.AEM;

import java.awt.Color;
import java.text.DecimalFormat;


import kang.AEM.BasicStructure.AEMStructure;
import kang.util.PrintUtil;
import kang.util.BasicTools;
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
	Grid grid2=new Grid("log wealth");     // the map to display the wealth distribution	
	
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
		double logwealth[]=new double[AHtemp.M];
		
		for(int i=0; i<AHtemp.M; i++)
		{
			logwealth[i]=Math.log(AHtemp.wealth[i]);
		}
		
		grid1.setColors(heatmap);
		grid1.registerData(AHtemp.L1, AHtemp.L2, AHtemp.wealth);
		
		grid2.setColors(heatmap);
		grid2.registerData(AHtemp.L1, AHtemp.L2, logwealth);
		
		params.set("totalwealth", AHtemp.totalwealth);
		params.set("meanwealth", AHtemp.meanwealth);
		params.set("order", AHtemp.order);
		
	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
	
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
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth, aem.psi);
			if(t%50==0)
			{
				PrintUtil.printlnToFile(entropy, t+1, aem.order, aem.order/M+entropyC, aem.psi);
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
	
	public void findsteadystate(AEMStructure aem, int steplimit, int seed, double mu, double gamma)   //the function to determine which is the best parameter to represent the steady state
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"steady state check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		String entropy="/Users/liukang2002507/Desktop/simulation/AEM/"+"steady state order <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";

		double entropyC=Math.log(aem.M);
		
		double[] ranki= new double[aem.M];
		int[] rankf= new int[aem.M];
		
		
		int index[]=new int[aem.M];
		//deltastep=100;
			
		for(int j=0; j<aem.M; j++)
		{
			index[j]=j;
			ranki[j]=j;
		}
		
		{
			index=Tools.BubbleSortIndex(aem.wealth, index, false);
			rankf=Tools.IndexToRank(index);
			for(int j=0; j<aem.M; j++)
			{
				ranki[j]=rankf[j];
			}
			
		}
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			
			PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			if(t%50==0)
			{
				
				index=Tools.BubbleSortIndex(aem.wealth, index, false);
				rankf=Tools.IndexToRank(index);
				double sumDrank=0;
				double sumLnw=0;
				for(int j=0; j<aem.M; j++)
				{
					sumDrank+=((ranki[j]-rankf[j])*(ranki[j]-rankf[j]));
					ranki[j]=rankf[j];
					sumLnw+=Math.log(aem.wealth[j]/aem.totalwealth);
				}
				PrintUtil.printlnToFile(entropy, t+1, aem.order, aem.order/M+entropyC, sumDrank, sumLnw, aem.wealth[index[0]]/aem.totalwealth);
			}
			
			params.set("time", t+1);
			Job.animate();

		}
		
		
		
	}
	
	public int steadystatetime(AEMStructure aem, int steplimit, int dt, int seed, double mu, double gamma)   //the function to determine which is the best parameter to represent the steady state
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		
		double lnw[]=new double[steplimit/dt];
		int turningpoint=0;
				
		for(int t=0; t<steplimit; t++)
		{
			
			aem.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			
			if(t%dt==0)
			{
				
				double sumLnw=0;
				for(int j=0; j<aem.M; j++)
				{
					sumLnw+=Math.log(aem.wealth[j]/aem.totalwealth);
				}
				lnw[t/dt]=sumLnw;
			}
			
			params.set("time", t+1);
			Job.animate();
		}
		turningpoint=Tools.turningpoint(lnw, steplimit-50*dt, steplimit, dt, true);
		
		params.set("Ngrowth", turningpoint);
		
		return turningpoint;
	}
	
	public void Rankchange(AEMStructure aem, int steplimit, int dt, int seed, double mu, double gamma)   //the function to determine which is the best parameter to represent the steady state
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		

		//double sumDrank=0;
		int number=0;
		
		int change[]=new int[aem.M]; 
		
		int displacement[]=new int[aem.M];
		int up[]=new int[aem.M];
		int down[]=new int[aem.M];
		
		double[] ranki= new double[aem.M];
		int[] rankf= new int[aem.M];
		
		
		int index[]=new int[aem.M];
		int oldindex[]=new int[aem.M];
		//deltastep=100;
			
		for(int j=0; j<aem.M; j++)
		{
			index[j]=j;
			ranki[j]=j;
			oldindex[j]=index[j];
			change[j]=0;
			up[j]=0;
			down[j]=0;
			displacement[j]=0;
		}
		
		{
			index=Tools.BubbleSortIndex(aem.wealth, index, false);
			rankf=Tools.IndexToRank(index);
			for(int j=0; j<aem.M; j++)
			{
				ranki[j]=rankf[j];
			}
			
		}	
		for(int t=0; t<steplimit; t++)
		{
			aem.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			
			if(t>=(steplimit-dt))
			{
				for(int j=0; j<aem.M; j++)
				{
					oldindex[j]=index[j];
				}
				
				index=Tools.BubbleSortIndex(aem.wealth, index, false);
				rankf=Tools.IndexToRank(index);

				for(int j=0; j<aem.M; j++)
				{
					
					if(index[j]!=oldindex[j])
						{
						change[j]++;
						int temp=0;
						temp=rankf[oldindex[j]]-j-1;
						displacement[j]+=temp;
						
						if(temp>0)
							down[j]+=temp;
						else
							up[j]+=temp;
						}
					ranki[j]=rankf[j];
				}
							
				number++;
			}
			
			params.set("time", t+1);
			Job.animate();
		}

		String path="/Users/liukang2002507/Desktop/simulation/AEM/"+"rank change <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", mu="+fmt.format(mu*10000)+", gamma= "+fmt.format(gamma*1000)+", dt= "+bmt.format(dt)+">.txt";
		
		for(int i=0; i<aem.M; i++)
		{
			PrintUtil.printlnToFile(path, i+1, change[i], number, displacement[i], up[i], down[i]);
		}
	
	}
	
	public int[] RankChange(AEMStructure aem, int steplimit, int dt, int seed, double mu, double gamma)   //the function to determine which is the best parameter to represent the steady state
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		int data[]=new int[6];    //data[0]=richest, data[1]=rank 5%, data[2]=rank 10%, data[3]=rank 50%, data[4]=rank 90%, data[5]=poorest
		
		

		//double sumDrank=0;
		int number=0;
		
		int change[]=new int[aem.M]; 
		
		int displacement[]=new int[aem.M];
		int up[]=new int[aem.M];
		int down[]=new int[aem.M];
		
		double[] ranki= new double[aem.M];
		int[] rankf= new int[aem.M];
		
		
		int index[]=new int[aem.M];
		int oldindex[]=new int[aem.M];
		//deltastep=100;
			
		for(int j=0; j<aem.M; j++)
		{
			index[j]=j;
			ranki[j]=j;
			oldindex[j]=index[j];
			change[j]=0;
			up[j]=0;
			down[j]=0;
			displacement[j]=0;
		}
		
		{
			index=Tools.BubbleSortIndex(aem.wealth, index, false);
			rankf=Tools.IndexToRank(index);
			for(int j=0; j<aem.M; j++)
			{
				ranki[j]=rankf[j];
			}
			
		}	
		for(int t=0; t<steplimit; t++)
		{
			aem.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			
			if(t>=(steplimit-dt))
			{
				for(int j=0; j<aem.M; j++)
				{
					oldindex[j]=index[j];
				}
				
				index=Tools.BubbleSortIndex(aem.wealth, index, false);
				rankf=Tools.IndexToRank(index);

				for(int j=0; j<aem.M; j++)
				{
					
					if(index[j]!=oldindex[j])
						{
						change[j]++;
						int temp=0;
						temp=rankf[oldindex[j]]-j-1;
						displacement[j]+=temp;
						
						if(temp>0)
							down[j]+=temp;
						else
							up[j]+=temp;
						}
					ranki[j]=rankf[j];
				}
							
				number++;
			}
			
			params.set("time", t+1);
			Job.animate();
		}
		data[0]=change[0];
		data[1]=change[aem.M/20-1];
		data[2]=change[aem.M/10-1];
		data[3]=change[aem.M/2-1];
		data[4]=change[aem.M*9/10-1];
		data[5]=change[aem.M-1];

		
        return data;
	
	}
	
	public double[] DistributionShape(AEMStructure aem, int steplimit, int dt, int seed, double mu, double gamma)   //the function to determine which is the best parameter to represent the steady state
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		
		double totallnw=0;
		double totalratio=0;
		double total10ratio=0;
		double totalmedian=0;
		double totalrich=0;
		double totalpoor=0;
		//double sumDrank=0;
		int number=0;
		
		double returndata[]=new double[6];  //[0]--ln of ratio between the richest and poorest, [1]--sumlnw, [2]--medium, [3]--ln of ratio of top 10% to 90%, [4] richest [5]poorest
		
		double[] ranki= new double[aem.M];
		int[] rankf= new int[aem.M];
		
		
		int index[]=new int[aem.M];
		//deltastep=100;
			
		for(int j=0; j<aem.M; j++)
		{
			index[j]=j;
			ranki[j]=j;
		}
		
		{
			index=Tools.BubbleSortIndex(aem.wealth, index, false);
			rankf=Tools.IndexToRank(index);
			for(int j=0; j<aem.M; j++)
			{
				ranki[j]=rankf[j];
			}
			
		}	
		for(int t=0; t<steplimit; t++)
		{
			aem.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			
			if(t>=(steplimit-dt))
			{
				index=Tools.BubbleSortIndex(aem.wealth, index, false);
				rankf=Tools.IndexToRank(index);
				double sumLnw=0;
				for(int j=0; j<aem.M; j++)
				{
					sumLnw+=Math.log(aem.wealth[j]/aem.totalwealth);
					ranki[j]=rankf[j];
					//sumDrank+=((ranki[j]-rankf[j])*(ranki[j]-rankf[j]));
				}
							
				
				totallnw+=sumLnw;
				totalratio+=(Math.log(aem.wealth[index[0]])-Math.log(aem.wealth[index[aem.M-1]]));
				total10ratio+=(Math.log(aem.wealth[index[aem.M/10-1]])-Math.log(aem.wealth[index[aem.M*9/10-1]]));
				totalmedian+=Math.log(aem.wealth[index[aem.M/2-1]]/aem.totalwealth);
				totalrich+=Math.log(aem.wealth[index[0]]/aem.totalwealth);
				totalpoor+=Math.log(aem.wealth[index[aem.M-1]]/aem.totalwealth);
				number++;
			}
			
			
			
			params.set("time", t+1);
			Job.animate();
		}

		returndata[0]=totalratio/number;
		returndata[1]=totallnw/number;
		returndata[2]=totalmedian/number;
		returndata[3]=total10ratio/number;
		returndata[4]=totalrich/number;
		returndata[5]=totalpoor/number;

		return returndata;
	
	}
	
	public void TauScanMu(AEMStructure aem, int steplimit, int dt, double minMu, double maxMu, double dMu, double gamma, Boolean oneseed)
	{
		
		String saveas;
		if(oneseed)
		{
			
				saveas="/Users/liukang2002507/Desktop/simulation/AEM/SteadyState/ScanMu data/1seed/"+"TauScanMu <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", gamma= "+fmt.format(gamma*1000)+", dt= "+bmt.format(dt)+">.txt";
		}
		else
		{
			
				saveas="/Users/liukang2002507/Desktop/simulation/AEM/SteadyState/ScanMu data/nseed/"+"TauScanMu <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", gamma= "+fmt.format(gamma*1000)+", dt= "+bmt.format(dt)+">.txt";
			
		}
		
		int randseed=1;
		
		
		{
			params.set("percent", percent);
			AHtemp=new AEMStructure(aem.L1,aem.L2,R,percent, tax, 0, aem.Ngrowth);
			
			for(double mutemp=minMu; mutemp<maxMu; mutemp+=dMu)
			{
				
				AHtemp.Uinitialization(1);
				params.set("mu", mutemp);
				
				if(oneseed)
				{
					randseed=1;
				}
				else
				{
					
					randseed++;
				}
		
				Job.animate();
				
                int tau=steadystatetime(AHtemp, steplimit, dt, randseed, mutemp, gamma);	
				PrintUtil.printlnToFile(saveas,  mutemp, tau, dt);
				
			}
		}

	}
	
	public void ShapeScanMu(AEMStructure aem, int steplimit, int dt, double minMu, double maxMu, double dMu, double gamma, Boolean oneseed)
	{
		
		String saveas;
		if(oneseed)
		{
	
				saveas="/Users/liukang2002507/Desktop/simulation/AEM/SteadyState/ScanMu data/1seed/"+"ShapeScanMu <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", gamma= "+fmt.format(gamma*1000)+", dt= "+bmt.format(dt)+">.txt";
			
		}
		else
		{
			
				saveas="/Users/liukang2002507/Desktop/simulation/AEM/SteadyState/ScanMu data/nseed/"+"ShapeScanMu <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", gamma= "+fmt.format(gamma*1000)+", dt= "+bmt.format(dt)+">.txt";
			
		}
		
		int randseed=1;
		
		
		{
			params.set("percent", percent);
			AHtemp=new AEMStructure(aem.L1,aem.L2,R,percent, tax, 0, aem.Ngrowth);
			
			for(double mutemp=minMu; mutemp<maxMu; mutemp+=dMu)
			{
				
				AHtemp.Uinitialization(1);
				params.set("mu", mutemp);
				
				if(oneseed)
				{
					randseed=1;
				}
				else
				{
					
					randseed++;
				}
		
				Job.animate();
				
				double shape[]=new double[4];
                shape=DistributionShape(AHtemp, steplimit, dt, randseed, mutemp, gamma);
                
				PrintUtil.printlnToFile(saveas,  mutemp, shape[0], shape[1], shape[2], shape[3]);
				
			}
		}

	}
	
	public void TauScanGamma(AEMStructure aem, int steplimit, int dt, double ming, double maxg, double dg, double mu, Boolean oneseed)
	{
		
		String saveas;
		if(oneseed)
		{
			
	    	saveas="/Users/liukang2002507/Desktop/simulation/AEM/SteadyState/ScanGamma data/1seed/"+"TauScanGamma <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", mu="+fmt.format(mu*10000)+", dt= "+bmt.format(dt)+">.txt";
		}
		else
		{
		
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/SteadyState/ScanGamma data/nseed/"+"TauScanGamma <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", mu="+fmt.format(mu*10000)+", dt= "+bmt.format(dt)+">.txt";

		}
		
		int randseed=1;
		
		
		{
			params.set("percent", percent);
			AHtemp=new AEMStructure(aem.L1,aem.L2,R,percent, tax, 0, aem.Ngrowth);
			
			for(double gtemp=ming; gtemp<maxg; gtemp+=dg)
			{
				
				AHtemp.Uinitialization(1);
				params.set("mu", gtemp);
				
				if(oneseed)
				{
					randseed=1;
				}
				else
				{
					
					randseed++;
				}
		
				Job.animate();
				
                int tau=steadystatetime(AHtemp, steplimit, dt, randseed, mu, gtemp);
				
				PrintUtil.printlnToFile(saveas,  gtemp, tau, dt);
				
			}
		}

	}
	
	public void RankChangeScanSize(AEMStructure aem, int steplimit, int dt, int minL, int maxL, int multiplier, double mu, double gamma, Boolean oneseed)
	{
		String saveas;
		//String path;
		if(oneseed)
		{
			
	    	saveas="/Users/liukang2002507/Desktop/simulation/AEM/SteadyState/ScanSize data/1seed/"+"RankChangeScanSize <t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", mu="+fmt.format(mu*10000)+", gamma= "+fmt.format(gamma*1000)+", dt= "+bmt.format(dt)+">.txt";
		}
		else
		{
		
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/SteadyState/ScanSize data/nseed/"+"RankChangeScanSize <t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", mu="+fmt.format(mu*10000)+", gamma= "+fmt.format(gamma*1000)+", dt= "+bmt.format(dt)+">.txt";

		}
		
		int randseed=1;
		
		for(int l=minL; l<maxL; l=(l*multiplier))
		{
			AHtemp=new AEMStructure(l,l,l,aem.percent,aem.tax, 0,aem.Ngrowth);
			AHtemp.Uinitialization(1);
			
			params.set("L", l);
			params.set("R", l);
			
			if(oneseed)
			{
				randseed=1;
			}
			else
			{
				
				randseed++;
			}	
			Job.animate();
			
			int rankchange[]=new int[6];
            
            rankchange=RankChange(AHtemp, steplimit, dt, randseed, mu, gamma);         
			
			PrintUtil.printlnToFile(saveas, l, rankchange[0], rankchange[1],rankchange[2],rankchange[3],rankchange[4],rankchange[5], dt);
			

			
			
			
		}
		
		
		
		
	}
	
	public void ShapeScanSize(AEMStructure aem, int steplimit, int dt, int minL, int maxL, int multiplier, double mu, double gamma, Boolean oneseed)
	{
		String saveas;
		String path;
		if(oneseed)
		{
			
	    	saveas="/Users/liukang2002507/Desktop/simulation/AEM/SteadyState/ScanSize data/1seed/"+"ShapeScanSize <t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", mu="+fmt.format(mu*10000)+", gamma= "+fmt.format(gamma*1000)+", dt= "+bmt.format(dt)+">.txt";
		}
		else
		{
		
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/SteadyState/ScanSize data/nseed/"+"ShapeScanSize <t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", mu="+fmt.format(mu*10000)+", gamma= "+fmt.format(gamma*1000)+", dt= "+bmt.format(dt)+">.txt";

		}
		
		int randseed=1;
		
		for(int l=minL; l<maxL; l=(l*multiplier))
		{
			AHtemp=new AEMStructure(l,l,l,aem.percent,aem.tax, 0,aem.Ngrowth);
			AHtemp.Uinitialization(1);
			
			params.set("L", l);
			params.set("R", l);
			
			if(oneseed)
			{
				randseed=1;
			}
			else
			{
				
				randseed++;
			}	
			Job.animate();
			
			double shape[]=new double[6];
            shape=DistributionShape(AHtemp, steplimit, dt, randseed, mu, gamma);
                       
			PrintUtil.printlnToFile(saveas, l, shape[0], shape[1], shape[2], shape[3], shape[4], shape[5]);
			
			if(oneseed)
			{
				
		    	path="/Users/liukang2002507/Desktop/simulation/AEM/SteadyState/ScanSize data/1seed/"+"size wealth <t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", mu="+fmt.format(mu*10000)+", gamma= "+fmt.format(gamma*1000)+", L= "+bmt.format(l)+">.txt";
			}
			else
			{
			
				path="/Users/liukang2002507/Desktop/simulation/AEM/SteadyState/ScanSize data/nseed/"+"size wealth <t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", mu="+fmt.format(mu*10000)+", gamma= "+fmt.format(gamma*1000)+", L= "+bmt.format(l)+">.txt";

			}
			
			output(AHtemp, path);
			
			
			
		}
		
		
		
		
	}
	
	public void ShapeScanGamma(AEMStructure aem, int steplimit, int dt, double ming, double maxg, double dg, double mu, Boolean oneseed)
	{
		
		String saveas;
		if(oneseed)
		{
			
	    	saveas="/Users/liukang2002507/Desktop/simulation/AEM/SteadyState/ScanGamma data/1seed/"+"ShapeScanGamma <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", mu="+fmt.format(mu*10000)+", dt= "+bmt.format(dt)+">.txt";
		}
		else
		{
		
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/SteadyState/ScanGamma data/nseed/"+"ShapeScanGamma <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", p= "+bmt.format(aem.percent*1000)+", mu="+fmt.format(mu*10000)+", dt= "+bmt.format(dt)+">.txt";

		}
		
		int randseed=1;
		
		
		{
			params.set("percent", percent);
			AHtemp=new AEMStructure(aem.L1,aem.L2,R,percent, tax, 0, aem.Ngrowth);
			
			for(double gtemp=ming; gtemp<maxg; gtemp+=dg)
			{
				
				AHtemp.Uinitialization(1);
				params.set("mu", gtemp);
				
				if(oneseed)
				{
					randseed=1;
				}
				else
				{
					
					randseed++;
				}
		
				Job.animate();
				

				double shape[]=new double[6];
                shape=DistributionShape(AHtemp, steplimit, dt, randseed, mu, gtemp);
                
				PrintUtil.printlnToFile(saveas,  gtemp, shape[0], shape[1], shape[2], shape[3], shape[4], shape[5]);
				
			}
		}

	}

	public void GINIrun(AEMStructure aem, int steplimit, int seed, double mu, double gamma)   //the function to determine which is the best parameter to represent the steady state
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String Gpath="/Users/liukang2002507/Desktop/simulation/AEM/"+"Gini <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit-100)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";

		double G=0;
		double[] Ldata=new double[aem.M];
		
		double[] ranki= new double[aem.M];
		int[] rankf= new int[aem.M];
		
		
		int index[]=new int[aem.M];
		//deltastep=100;
			
		for(int j=0; j<aem.M; j++)
		{
			index[j]=j;
			ranki[j]=j;
		}
		
		{
			index=Tools.BubbleSortIndex(aem.wealth, index, true);   // here the sorting needs to be ascending
			rankf=Tools.IndexToRank(index);
			
			
			for(int j=0; j<aem.M; j++)
			{
				ranki[j]=rankf[j];
			}
			
		}
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			
			
			if(t%1==0)
			{
				
				index=Tools.BubbleSortIndex(aem.wealth, index, true);
				Ldata=Tools.LorentzData(aem.wealth, index, true);
				
				G=Tools.Gini(Ldata);
				rankf=Tools.IndexToRank(index);
				
				PrintUtil.printlnToFile(Gpath, t+1, aem.totalwealth, G);
				
			}
			
			if(t==(steplimit-100))
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"Gini wealth <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
				
			}
			
			params.set("time", t+1);
			Job.animate();

		}
		
		
		
	}
	
	public void populationrun(AEMStructure aem, int steplimit, int seed, double mu, double gamma)   //the function to determine which is the best parameter to represent the steady state
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"population check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		String entropy="/Users/liukang2002507/Desktop/simulation/AEM/"+"population order <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";

		double entropyC=Math.log(aem.M);
		
		double[] ranki= new double[aem.M];
		int[] rankf= new int[aem.M];
		
		
		int index[]=new int[aem.M];
		//deltastep=100;
			
		for(int j=0; j<aem.M; j++)
		{
			index[j]=j;
			ranki[j]=j;
		}
		
		{
			index=Tools.BubbleSortIndex(aem.wealth, index, false);
			rankf=Tools.IndexToRank(index);
			for(int j=0; j<aem.M; j++)
			{
				ranki[j]=rankf[j];
			}
			
		}
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			
			
			if(t%50==0)
			{
				
				index=Tools.BubbleSortIndex(aem.wealth, index, false);
				rankf=Tools.IndexToRank(index);
				double sumDrank=0;
				double sumLnw=0;
				for(int j=0; j<aem.M; j++)
				{
					sumDrank+=((ranki[j]-rankf[j])*(ranki[j]-rankf[j]));
					ranki[j]=rankf[j];
					sumLnw+=Math.log(aem.wealth[j]/aem.totalwealth);
				}
				PrintUtil.printlnToFile(total, t+1, aem.totalwealth, aem.wealth[index[0]],  aem.wealth[index[aem.M/20-1]], aem.wealth[index[aem.M/10-1]],  aem.wealth[index[aem.M/2-1]], aem.wealth[index[aem.M*9/10-1]], aem.wealth[index[aem.M-1]]);
				
				PrintUtil.printlnToFile(entropy, t+1, aem.order, aem.order/M+entropyC, sumDrank, sumLnw, aem.wealth[index[0]]/aem.totalwealth);
			}
			
			params.set("time", t+1);
			Job.animate();

		}
		
		
		
	}
	
	public void IncomeRun(AEMStructure aem, int steplimit, int time, int seed, double mu, double gamma)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/AEM/"+"Income check <L="+fmt.format(aem.L1)+", dt="+bmt.format(time)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		String entropy="/Users/liukang2002507/Desktop/simulation/AEM/"+"Income order <L="+fmt.format(aem.L1)+", dt="+bmt.format(time)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
    
		double entropyC=Math.log(aem.M);
		
		double wealthrecord[]=new double[aem.M];
		//double income[]=new double[aem.M];
		
		
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
			if(t==(500-time))
			{
				for(int jj=0; jj<aem.M; jj++)
				{
					wealthrecord[jj]=aem.wealth[jj];
				}
			}
			if(t==500)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"income <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", dt="+bmt.format(time)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        incomeoutput(aem,  wealthrecord, saveas);
			}
			
			
			if(t==(2000-time))
			{
				for(int jj=0; jj<aem.M; jj++)
				{
					wealthrecord[jj]=aem.wealth[jj];
				}
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"income <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", dt="+bmt.format(time)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
				incomeoutput(aem,  wealthrecord, saveas);
			}
			
			if(t==(5000-time))
			{
				for(int jj=0; jj<aem.M; jj++)
				{
					wealthrecord[jj]=aem.wealth[jj];
				}
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"income <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", dt="+bmt.format(time)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
				incomeoutput(aem,  wealthrecord, saveas);
			}
			
			if(t==(10000-time))
			{
				for(int jj=0; jj<aem.M; jj++)
				{
					wealthrecord[jj]=aem.wealth[jj];
				}
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"income <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", dt="+bmt.format(time)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
				incomeoutput(aem,  wealthrecord, saveas);
			}
			
			if(t==(50000-time))
			{
				for(int jj=0; jj<aem.M; jj++)
				{
					wealthrecord[jj]=aem.wealth[jj];
				}
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"income <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", dt="+bmt.format(time)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
				incomeoutput(aem,  wealthrecord, saveas);
			}	
			
			if(t==(100000-time))
			{
				for(int jj=0; jj<aem.M; jj++)
				{
					wealthrecord[jj]=aem.wealth[jj];
				}
			}
			if(t==100000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"income <L="+fmt.format(aem.L1)+", t="+bmt.format(t)+", dt="+bmt.format(time)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
				incomeoutput(aem,  wealthrecord, saveas);
			}

		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
	}
	
	public void RankFluctuation(AEMStructure aem, int steplimit, int deltastep, int seed, double mu, double gamma)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String rankpath="/Users/liukang2002507/Desktop/simulation/AEM/"+"Rank Fluctuation <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", dt="+fmt.format(deltastep)+">.txt";
        double[] ranki= new double[aem.M];
		int[] rankf= new int[aem.M];
		
		double totaldrank=0;
		double totaldRMI=0; // rank mobility index=(Rf-Ri)\(Rf+Ri)
		double SDrank=0;
		double SDRMI=0;
		
		double correlation=0;
		
		
		int index[]=new int[aem.M];
		//deltastep=100;
			
		for(int j=0; j<aem.M; j++)
		{
			index[j]=j;
			ranki[j]=j;
		}
		
		for(int t=0; t<steplimit; t++)
		{
			aem.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			//PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			if(t%deltastep==0)
			{
				index=Tools.BubbleSortIndex(aem.wealth, index, false);
				rankf=Tools.IndexToRank(index);

				totaldrank=0;
				totaldRMI=0;
				for(int jj=0; jj<aem.M; jj++)
				{
					totaldrank+=((rankf[jj]-ranki[jj])*(rankf[jj]-ranki[jj]));
					totaldRMI+=((rankf[jj]-ranki[jj])*(rankf[jj]-ranki[jj]))/((rankf[jj]+ranki[jj])*(rankf[jj]+ranki[jj]));
				}
				
				SDrank=Math.sqrt(totaldrank/aem.M);
				SDRMI=Math.sqrt(totaldRMI/aem.M);
				
				correlation=Tools.Correlation(ranki, rankf, aem.M);
				PrintUtil.printlnToFile(rankpath, t, correlation, SDrank, SDRMI, aem.totalwealth, aem.sumtrading, aem.sumflow);
				for(int jf=0; jf<aem.M; jf++)
				{
					ranki[jf]=rankf[jf];           //reset the initial rank
				}

			}
			params.set("time", t+1);
			Job.animate();
		}
				
	}
	
	public void RankCorrelation(AEMStructure aem, int start, int steplimit, int deltastep, int seed, double mu, double gamma)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String rankpath="/Users/liukang2002507/Desktop/simulation/AEM/"+"Rank Correlation <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", start="+fmt.format(start)+", dt="+fmt.format(deltastep)+">.txt";
		
		//double[] ranks= new double[aem.M];
		double[] ranki= new double[aem.M];
		int[] rankf= new int[aem.M];
		
		double totaldrank=0;
		double totaldRMI=0; // rank mobility index=(Rf-Ri)\(Rf+Ri)
		double SDrank=0;
		double SDRMI=0;
		
		double correlation=0;
		
		
		int index[]=new int[aem.M];
		//deltastep=100;
			
		for(int j=0; j<aem.M; j++)
		{
			index[j]=j;
			ranki[j]=j;
		}
		
		for(int t=0; t<steplimit; t++)
		{
			aem.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			//PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			if(t==start)
			{
				index=Tools.BubbleSortIndex(aem.wealth, index, false);
				rankf=Tools.IndexToRank(index);
				for(int j=0; j<aem.M; j++)
				{
					ranki[j]=rankf[j];
				}
				
			}
			
			
			if((t>start)&&((t-start)%deltastep==0))
			{
				index=Tools.BubbleSortIndex(aem.wealth, index, false);
				rankf=Tools.IndexToRank(index);

				totaldrank=0;
				totaldRMI=0;
				for(int jj=0; jj<aem.M; jj++)
				{
					totaldrank+=((rankf[jj]-ranki[jj])*(rankf[jj]-ranki[jj]));
					totaldRMI+=((rankf[jj]-ranki[jj])*(rankf[jj]-ranki[jj]))/((rankf[jj]+ranki[jj])*(rankf[jj]+ranki[jj]));
				}
				
				SDrank=Math.sqrt(totaldrank/aem.M);
				SDRMI=Math.sqrt(totaldRMI/aem.M);
				
				correlation=Tools.Correlation(ranki, rankf, aem.M);
				PrintUtil.printlnToFile(rankpath, t, correlation, SDrank, SDRMI, aem.totalwealth, aem.sumtrading, aem.sumflow);

			}
			params.set("time", t+1);
			Job.animate();
		}
				
	}
	
	public void Metric(AEMStructure aem, int start, int steplimit, int deltastep, int seed, double mu, double gamma)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String path="/Users/liukang2002507/Desktop/simulation/AEM/"+"Metric <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", start="+fmt.format(start)+", dt="+fmt.format(deltastep)+">.txt";
		
		//double[] ranks= new double[aem.M];
		double[] ranki= new double[aem.M];
		int[] rankf= new int[aem.M];
		
		double totaldrank=0;
		double totaldRMI=0; // rank mobility index=(Rf-Ri)\(Rf+Ri)
		//double SDrank=0;
		//double SDRMI=0;
		//double correlation=0;
		
	
		//metric variables
		
		double[] totalR=new double[aem.M];
		double SumR=0;
		double SavgR=0;
		double metricR=0;
		
		double[] totalW=new double[aem.M];      //rescaled wealth
		double SumW=0;
		double SavgW=0;
		double metricW=0;
		
		double[] totalAW=new double[aem.M];      //absolute value of wealth
		double SumAW=0;
		double SavgAW=0;
		double metricAW=0;
		
		
		double[] totallnW=new double[aem.M];
		double SumlnW=0;
		double SavglnW=0;
		double metriclnW=0;
		
		double time=0;
		
		
		int index[]=new int[aem.M];
		//deltastep=100;
			
		for(int j=0; j<aem.M; j++)
		{
			index[j]=j;
			ranki[j]=j;
			totalR[j]=0;
			totalW[j]=0;
			totalAW[j]=0;
			totallnW[j]=0;
		}
		
		for(int t=0; t<steplimit; t++)
		{
			aem.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			//PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			if(t==start)
			{
				index=Tools.BubbleSortIndex(aem.wealth, index, false);
				rankf=Tools.IndexToRank(index);
				for(int j=0; j<aem.M; j++)
				{
					ranki[j]=rankf[j];
				}
				
			}
			
			if((t>start)&&((t-start)%deltastep==0))
			{
				index=Tools.BubbleSortIndex(aem.wealth, index, false);
				rankf=Tools.IndexToRank(index);

				totaldrank=0;
				totaldRMI=0;
				for(int jj=0; jj<aem.M; jj++)
				{
				
					totalR[jj]+=rankf[jj];
					totalW[jj]+=(aem.wealth[jj]/aem.totalwealth);
					totalAW[jj]+=(aem.wealth[jj]);
					totallnW[jj]+=Math.log(aem.wealth[jj]/aem.totalwealth);
				    
					time++;
				}
				SumR=0;
				SumW=0;
				SumAW=0;
				SumlnW=0;
				
				
				for(int ii=0; ii<aem.M; ii++)
				{
					SumR+=(totalR[ii]/time);
					SumW+=(totalW[ii]/time);
					SumAW+=(totalAW[ii]/time);
					SumlnW+=(totallnW[ii]/time);
				}
				
				SavgR=SumR/aem.M;
				SavgW=SumW/aem.M;
				SavgAW=SumAW/aem.M;
				SavglnW=SumlnW/aem.M;
				
				double SMR=0;
				double SMW=0;
				double SMAW=0;
				double SMlnW=0;
				
			    for(int kk=0; kk<aem.M; kk++)
			    {
			    	SMR+=((totalR[kk]/time-SavgR)*(totalR[kk]/time-SavgR));
			    	SMW+=((totalW[kk]/time-SavgW)*(totalW[kk]/time-SavgW));
			    	SMAW+=((totalAW[kk]/time-SavgAW)*(totalAW[kk]/time-SavgAW));
			    	SMlnW+=((totallnW[kk]/time-SavglnW)*(totallnW[kk]/time-SavglnW));
			    }
				
			    metricR=SMR/aem.M;
			    metricW=SMW/aem.M;
			    metricAW=SMAW/aem.M;
			    metriclnW=SMlnW/aem.M;
				
				
				//SDrank=Math.sqrt(totaldrank/aem.M);
				//SDRMI=Math.sqrt(totaldRMI/aem.M);
				//correlation=Tools.Correlation(ranki, rankf, aem.M);
			    
				PrintUtil.printlnToFile(path, t, metricR, metricW, metriclnW, metricAW,aem.totalwealth, aem.sumtrading, aem.sumflow);

			}
			params.set("time", t+1);
			Job.animate();
		}
				
	}
	
	public void Rankrun(AEMStructure aem, int steplimit, int seed, double mu, double gamma)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String rankpath="/Users/liukang2002507/Desktop/simulation/AEM/"+"Rank check <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		String wealth="/Users/liukang2002507/Desktop/simulation/AEM/"+"Rank order <L="+fmt.format(aem.L1)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";

		
		int index[]=new int[aem.M];
		int deltastep=100;
		
		int rank[]=new int[aem.M];
		//int rank[][]=new int[steplimit/deltastep][aem.M];
		
		for(int j=0; j<aem.M; j++)
		{
			index[j]=j;
		}
		
		
		
		for(int t=0; t<steplimit; t++)
		{
			
			aem.SublinearTS(rand, aem.percent, aem.tax, aem.alpha, ngrowth, gamma);
			ngrowth=ngrowth*(1+mu);
			//PrintUtil.printlnToFile(total, t+1, aem.totalwealth);
			if(t%deltastep==0)
			{
				index=Tools.BubbleSortIndex(aem.wealth, index, false);
			}
			
			params.set("time", t+1);
			Job.animate();
			
		}
		
		index=Tools.BubbleSortIndex(aem.wealth, index, false);
		
		rank=Tools.IndexToRank(index);
		for(int i=0; i<aem.M; i++)
		{
			PrintUtil.printlnToFile(rankpath, i, rank[i], aem.wealth[i]);
			PrintUtil.printlnToFile(wealth, i, aem.wealth[index[i]]);
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
			if(maxMu<0.05)
				saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanMuP data/1seed/"+"Smallrange Scan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", n= "+bmt.format(number)+">.txt";
			else
				saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanMuP data/1seed/"+"Scan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", n= "+bmt.format(number)+">.txt";
		}
		else
		{
			if(maxMu<0.05)
				saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanMuP data/nseed/"+"Smallrange Scan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", n= "+bmt.format(number)+">.txt";
			else
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
	
	
    public void RankScanMuP(AEMStructure aem, int steplimit, double minMu, double maxMu, double dMu, double gamma, Boolean oneseed, int Index)
	{
		double minP=0.01;
		double maxP=0.40;
		double dP=0.01;
		
		
		String saveas;
		if(oneseed)
		{
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanMuP data/1seed/"+"RankScan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", I= "+bmt.format(Index)+">.txt";
		}
		else
		{
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanMuP data/nseed/"+"RankScan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", I= "+bmt.format(Index)+">.txt";
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
				
				
				int index[]=new int[aem.M];
								
				for(int j=0; j<aem.M; j++)
				{
					index[j]=j;
				}
				
				index=Tools.BubbleSortIndex(AHtemp.wealth, index, false);
				
				
				double Iwealth=AHtemp.wealth[index[Index-1]];
				PrintUtil.printlnToFile(saveas, percent, mutemp, Iwealth);
				
			}
		}

	}
    
    public void RankScanMuP(AEMStructure aem, int steplimit, double minMu, double maxMu, double dMu, double gamma, Boolean oneseed)
	{
		double minP=0.01;
		double maxP=0.40;
		double dP=0.01;
		
		
		String saveas;
		if(oneseed)
		{
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanMuP data/1seed/"+"RankScan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		}
		else
		{
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanMuP data/nseed/"+"RankScan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
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
				
				
				int index[]=new int[aem.M];
								
				for(int j=0; j<aem.M; j++)
				{
					index[j]=j;
				}
				
				index=Tools.BubbleSortIndex(AHtemp.wealth, index, false);
				
				
				
				PrintUtil.printlnToFile(saveas, percent, mutemp, AHtemp.wealth[index[0]],  AHtemp.wealth[index[aem.M/20-1]], AHtemp.wealth[index[aem.M/10-1]],  AHtemp.wealth[index[aem.M/2-1]], AHtemp.wealth[index[aem.M*9/10-1]], AHtemp.wealth[index[aem.M-1]]);
				
			}
		}

	}
    
	
    public void GINIScanMuP(AEMStructure aem, int steplimit, double minMu, double maxMu, double dMu, double gamma, Boolean oneseed)
	{
		double minP=0.01;
		double maxP=0.40;
		double dP=0.01;
		
		String saveas;
		if(oneseed)
		{
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanMuP data/1seed/"+"GiniScan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		}
		else
		{
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanMuP data/nseed/"+"GiniScan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
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
				
				double G=0;
				double Psi=0;
				double totalW=0;
				double[] Ldata= new double[aem.M];
				
				int index[]=new int[aem.M];
								
				for(int j=0; j<aem.M; j++)
				{
					index[j]=j;
					totalW+=AHtemp.wealth[j];    //calculate the total wealth
				}
				
				for(int k=0; k<aem.M; k++)
				{
					Psi+=Math.log(AHtemp.wealth[k]/totalW);
				}
				
				
				index=Tools.BubbleSortIndex(AHtemp.wealth, index, true);		
				Ldata=Tools.LorentzData(AHtemp.wealth, index, true);
								
				G=Tools.Gini(Ldata);
				
				
				PrintUtil.printlnToFile(saveas, percent, mutemp, G, Psi, AHtemp.wealth[index[0]], AHtemp.wealth[index[aem.M-1]]);
				
			}
		}

	}
    
    public void richScanMuP(AEMStructure aem, int steplimit, double minMu, double maxMu, double dMu, double gamma, Boolean oneseed, int number)
	{
		double minP=0.01;
		double maxP=0.40;
		double dP=0.01;
		
		
		String saveas;
		if(oneseed)
		{
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanMuP data/1seed/"+"richScan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", n= "+bmt.format(number)+">.txt";
		}
		else
		{
			saveas="/Users/liukang2002507/Desktop/simulation/AEM/ScanMuP data/nseed/"+"richScan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", n= "+bmt.format(number)+">.txt";
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
				poordata=findrich(AHtemp.wealth, number);
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
	
	public void incomeoutput(AEMStructure aem, double[] record, String path)
	{
		for(int i=0; i<aem.M; i++)
		{
			PrintUtil.printlnToFile(path, i+1, aem.wealth[i]-record[i]);
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

		Globaltrading.frameTogether("Display", grid1, grid2);

		params.add("L", 50);
		params.add("R", 50);
		
	    params.add("tax",0.00);
		params.add("percent", 0.10);
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
	    
	    
	    //sublinearexp(AHtemp, 51000, 9, mu, 0.000);
	    //sublinearexp(AHtemp, 101000, 9, mu, 1.100);
	    
	    //findsteadystate(AHtemp, 101000, 9, mu, 0.500);
	    //populationrun(AHtemp, 100000, 9, mu, 0.000);
	    
	    GINIrun(AHtemp, 3000, 9, mu, 1.000);
	    //GINIScanMuP(AH, 1000, 0.0000, 0.4000, 0.0100, 0.000, true);
	    
	    
	    //Rankchange(AHtemp, 150000, 20000, 9, mu, 1.100);
	    
	    
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
	    
	    //RankScanMuP(AH, 100, 0.0000, 0.4000, 0.0100, 1.000, true, 0001);
	    //RankScanMuP(AH, 1000, 0.0000, 0.4000, 0.0100, 1.100, true);
	    
	    
	    //ScanMuP(AH, 100, 0.0000, 0.04000, 0.0010, 1.100, true, 25);
	    //richScanMuP(AH, 1000, 0.0000, 0.4000, 0.0100, 1.100, true, 25);
	    //ScanP(AH, 20000, 0.02, 0.40, 0.02, 1, false, 1);
	    //TauScanMu(AH, 100000, 20, 0.0001, 0.0020, 0.0001, 0.50, true);
	    
	    //ShapeScanMu(AH, 100000, 20, 0.0001, 0.0020, 0.0001, 0.50, true);
	    //ShapeScanSize(AH, 50000, 20, 5, 400, 2, mu, 1.100, true);
	    //RankChangeScanSize(AH, 50000, 20000, 9, 100, 2, mu, 0.000, true);
	    
	    
	    //TauScanGamma(AH, 100000, 10, 0.000, 0.800, 0.100, mu, true);
	    //TauScanGamma(AH, 100000, 10, 0.820, 0.940, 0.020, mu, true);
	    //TauScanGamma(AH, 100000, 10, 0.945, 0.985, 0.005, mu, true);
	    
	    //ShapeScanGamma(AH, 100000, 10, 0.000, 0.800, 0.100, mu, true);
	    //ShapeScanGamma(AH, 100000, 10, 0.810, 0.940, 0.010, mu, true);
	    //ShapeScanGamma(AH, 100000, 10, 0.945, 0.985, 0.005, mu, true);
	    
	    
	    //Rankrun(AHtemp, 10000, 1, mu, 0.000);
	    //RankCorrelation(AHtemp, 100000, 200100, 1000, 1, mu, 0.900);
	    
	   
	   // Metric(AHtemp, 50000, 500000, 1000, 1, mu, 0.500);
	    
	    //IncomeRun(AHtemp, 50100, 200, 1, mu, 0.100);
	    
	    //RankFluctuation(AHtemp, 100100, 100, 1, mu, 0.300);
	    
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