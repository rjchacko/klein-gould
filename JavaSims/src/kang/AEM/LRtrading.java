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


public class LRtrading extends Simulation{
	
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
		
		
		double logwealth[]=new double[L*L];
		
		for(int i=0; i<L*L; i++)
		{
			logwealth[i]=Math.log(AHtemp.wealth[i]);
		}
		
		grid1.setColors(heatmap);
		grid1.registerData(L, L, AHtemp.wealth);
		
		grid2.setColors(heatmap);
		grid2.registerData(L, L, logwealth);
		
		params.set("totalwealth", AHtemp.totalwealth);
		params.set("meanwealth", AHtemp.meanwealth);
		params.set("order", AHtemp.order);
		
	}
	
	public void clear()
	{
		grid1.clear();
		grid2.clear();
	}
	
	public static void main (String[] LRtrading){
		new Control(new LRtrading(), "Kang Liu's long range trading AEM" );
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
	
	public void output(AEMStructure aem, String path)
	{
		for(int i=0; i<aem.M; i++)
		{
			PrintUtil.printlnToFile(path, i+1, aem.wealth[i]);
		}
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
	
	public void sublinearexp(AEMStructure aem, int steplimit, int seed, double mu, double gamma, boolean picture)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String total="/Users/liukang2002507/Desktop/simulation/LRAEM/"+"sub check <L="+fmt.format(aem.L1)+", R="+fmt.format(aem.R)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		String entropy="/Users/liukang2002507/Desktop/simulation/LRAEM/"+"sub order <L="+fmt.format(aem.L1)+", R="+fmt.format(aem.R)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
        String picpath="/Users/liukang2002507/Desktop/simulation/LRAEM/"+"<L="+fmt.format(aem.L1)+", R="+fmt.format(aem.R)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+">/";
		
		double entropyC=Math.log(aem.M);
		
		int intgamma=(int)(gamma*1000);
		
		
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
				String saveas="/Users/liukang2002507/Desktop/simulation/LRAEM/"+"sub wealth <L="+fmt.format(aem.L1)+", R="+fmt.format(aem.R)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
		        if(picture)
		        	{
		        	String path1=picpath+"wealth_";
		        	Tools.Picture(grid1, intgamma, t, path1);
		        	String path2=picpath+"logw_";
		        	Tools.Picture(grid2, intgamma, t, path2);
		        	}
			}
			if(t==2000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/LRAEM/"+"sub wealth <L="+fmt.format(aem.L1)+", R="+fmt.format(aem.R)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
		        if(picture)
	        	{
		        	String path1=picpath+"wealth_";
		        	Tools.Picture(grid1, intgamma, t, path1);
		        	String path2=picpath+"logw_";
		        	Tools.Picture(grid2, intgamma, t, path2);
	        	}
			}
			if(t==5000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/LRAEM/"+"sub wealth <L="+fmt.format(aem.L1)+", R="+fmt.format(aem.R)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
		        if(picture)
	        	{
		        	String path1=picpath+"wealth_";
		        	Tools.Picture(grid1, intgamma, t, path1);
		        	String path2=picpath+"logw_";
		        	Tools.Picture(grid2, intgamma, t, path2);
	        	}
			}
			if(t==10000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/LRAEM/"+"sub wealth <L="+fmt.format(aem.L1)+", R="+fmt.format(aem.R)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
		        if(picture)
	        	{
		        	String path1=picpath+"wealth_";
		        	Tools.Picture(grid1, intgamma, t, path1);
		        	String path2=picpath+"logw_";
		        	Tools.Picture(grid2, intgamma, t, path2);
	        	}
			}
			if(t==50000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/LRAEM/"+"sub wealth <L="+fmt.format(aem.L1)+", R="+fmt.format(aem.R)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
		        if(picture)
	        	{
		        	String path1=picpath+"wealth_";
		        	Tools.Picture(grid1, intgamma, t, path1);
		        	String path2=picpath+"logw_";
		        	Tools.Picture(grid2, intgamma, t, path2);
	        	}
			}	
			if(t==100000)
			{
				String saveas="/Users/liukang2002507/Desktop/simulation/AEM/"+"sub wealth <L="+fmt.format(aem.L1)+", R="+fmt.format(aem.R)+", t="+bmt.format(t)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+">.txt";
		        output(aem, saveas);
		        if(picture)
	        	{
		        	String path1=picpath+"wealth_";
		        	Tools.Picture(grid1, intgamma, t, path1);
		        	String path2=picpath+"logw_";
		        	Tools.Picture(grid2, intgamma, t, path2);
	        	}
			}

		}
		
		
		//String log="/Users/liukang2002507/Desktop/simulation/AEM/log "+"<"+bmt.format(steplimit)+">.txt";
		//PrintUtil.printlnToFile(log, aem.percent, aem.growth, aem.meanwealth);
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
				saveas="/Users/liukang2002507/Desktop/simulation/LRAEM/ScanMuP data/1seed/"+"Smallrange Scan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", n= "+bmt.format(number)+">.txt";
			else
				saveas="/Users/liukang2002507/Desktop/simulation/LRAEM/ScanMuP data/1seed/"+"Scan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", n= "+bmt.format(number)+">.txt";
		}
		else
		{
			if(maxMu<0.05)
				saveas="/Users/liukang2002507/Desktop/simulation/LRAEM/ScanMuP data/nseed/"+"Smallrange Scan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", n= "+bmt.format(number)+">.txt";
			else
				saveas="/Users/liukang2002507/Desktop/simulation/LRAEM/ScanMuP data/nseed/"+"Scan <L="+fmt.format(aem.L1)+", t="+bmt.format(steplimit)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", n= "+bmt.format(number)+">.txt";

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
	
	public void Metric(AEMStructure aem, int start, int steplimit, int deltastep, int seed, double mu, double gamma)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String path="/Users/liukang2002507/Desktop/simulation/LRAEM/"+"Metric <L="+fmt.format(aem.L1)+", R="+fmt.format(aem.R)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", start="+fmt.format(start)+", dt="+fmt.format(deltastep)+">.txt";
		
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
		
		double[] totalW=new double[aem.M];
		double SumW=0;
		double SavgW=0;
		double metricW=0;
		
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
					totallnW[jj]+=Math.log(aem.wealth[jj]/aem.totalwealth);
				    
					time++;
				}
				SumR=0;
				SumW=0;
				SumlnW=0;
				
				
				for(int ii=0; ii<aem.M; ii++)
				{
					SumR+=(totalR[ii]/time);
					SumW+=(totalW[ii]/time);
					SumlnW+=(totallnW[ii]/time);
				}
				
				SavgR=SumR/aem.M;
				SavgW=SumW/aem.M;
				SavglnW=SumlnW/aem.M;
				
				double SMR=0;
				double SMW=0;
				double SMlnW=0;
				
			    for(int kk=0; kk<aem.M; kk++)
			    {
			    	SMR+=((totalR[kk]/time-SavgR)*(totalR[kk]/time-SavgR));
			    	SMW+=((totalW[kk]/time-SavgW)*(totalW[kk]/time-SavgW));
			    	SMlnW+=((totallnW[kk]/time-SavglnW)*(totallnW[kk]/time-SavglnW));
			    }
				
			    metricR=SMR/aem.M;
			    metricW=SMW/aem.M;
			    metriclnW=SMlnW/aem.M;
				
				
				//SDrank=Math.sqrt(totaldrank/aem.M);
				//SDRMI=Math.sqrt(totaldRMI/aem.M);
				//correlation=Tools.Correlation(ranki, rankf, aem.M);
			    
				PrintUtil.printlnToFile(path, t, metricR, metricW, metriclnW, aem.totalwealth, aem.sumtrading, aem.sumflow);

			}
			params.set("time", t+1);
			Job.animate();
		}
				
	}
	
	public void RankCorrelation(AEMStructure aem, int start, int steplimit, int deltastep, int seed, double mu, double gamma)
	{
		double ngrowth=aem.Ngrowth;
		Random rand= new Random(seed);
		String rankpath="/Users/liukang2002507/Desktop/simulation/LRAEM/"+"Rank Correlation <L="+fmt.format(aem.L1)+", R="+fmt.format(aem.R)+", mu="+fmt.format(mu*10000)+", p= "+fmt.format(aem.percent*1000)+", growth= "+bmt.format(aem.Ngrowth)+", gamma= "+fmt.format(gamma*1000)+", start="+fmt.format(start)+", dt="+fmt.format(deltastep)+">.txt";
		
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
	
	public void load(Control LRtrading)
	{

		LRtrading.frameTogether("Display", grid1, grid2);

		params.add("L", 50);
		params.add("R", 0);
		
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
	    
	    //Metric(AHtemp, 50000, 200000, 1000, 1, mu, 0.100);
	    
	    
	    sublinearexp(AHtemp, 101000, 9, mu, 0.000, true);
	    
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
	    
	    
	    //ScanMuP(AH, 100, 0.0000, 0.04000, 0.0010, 1.100, true, 25);
	    //richScanMuP(AH, 100, 0.0000, 0.4000, 0.0100, 0.900, true, 25);
	    //ScanP(AH, 20000, 0.02, 0.40, 0.02, 1, false, 1);
	    
	    
	    
	    //Rankrun(AHtemp, 10000, 1, mu, 0.000);
	    //RankCorrelation(AHtemp, 100000, 200100, 1000, 1, mu, 0.900);
	    
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
	
	