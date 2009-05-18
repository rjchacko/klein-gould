package rachele.ising.dim2MC.apps;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static scikit.util.Utilities.format;

import java.awt.Color;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.Random;

import rachele.ising.dim2MC.IsingLR;
import rachele.util.FileUtil;
import rachele.util.FourierTransformer;
import scikit.dataset.Accumulator;
import scikit.dataset.DatasetBuffer;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
//import scikit.jobs.params.FileValue;

/**
* 
* Monte Carlo Simulation to produce structure factor vs time averages for the 
* clumps to stripes transition.
* This has code to characterize if the stripes are hor or vertical and
* collects only the dominant direction data.
* 
*/

public class MCClumpsStripesApp extends Simulation{

		Grid grid = new Grid("Long Range Ising Model");
		Plot sftPlot = new Plot("sf_t plot perpendicular");
		Plot sftParaPlot = new Plot("sf_t plot parallel");
		Plot sftAvePlot = new Plot("sf_t plot perpendicular Ave");
		Plot sftAveParaPlot = new Plot("sf_t plot parallel Ave");
		Plot compare = new Plot ("Horizontal vertical compare");
		Plot histogramPlot = new Plot("transition distribution");
		Plot normHistPlot = new Plot("normed transition distribution");
		int dx;
		IsingLR sim;
		public FourierTransformer fft;
		double [] sFactor;
		Random random = new Random();
		Accumulator sf_kAcc;
		int accNo = 10;
		Accumulator [] para_sf_tAveAcc = new Accumulator [accNo];  // parallel to stripe direction
		Accumulator [] para_sf_tAcc = new Accumulator [accNo];
		Accumulator [] perp_sf_tAveAcc = new Accumulator [accNo];  // perpendicular to stripe direction
		Accumulator [] perp_sf_tAcc = new Accumulator [accNo];
		Accumulator [] tempHor  = new Accumulator [accNo];
		Accumulator [] tempVer  = new Accumulator [accNo];
		Histogram normHist;
		Histogram hist;
		
		int [] sfLabel = new int [accNo];
		boolean init1D;	
		
		public static void main(String[] args) {
			new Control(new MCClumpsStripesApp(), "Clumps to Stripes");
		}

		public void load(Control c) {

//			c.frame(grid);
			c.frameTogether("Plots",sftAvePlot, sftAveParaPlot, grid, sftPlot, sftParaPlot, compare,histogramPlot, normHistPlot);
			params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/lraim/stripeToClumpInvestigation/mcResults/noConservedOP/testRuns"));
//			params.add("Input 1D File",new FileValue("/Users/erdomi/data/lraim/configs1D/L128R50T0-04h0"));
			params.addm("Dynamics", new ChoiceValue("Ising Glauber","Kawasaki Glauber", "Kawasaki Metropolis",  "Ising Metropolis"));
//			params.addm("init", new ChoiceValue( "Init Clumps from random"));
			params.add("Random seed", 0);
			params.add("L", 1<<6);
			params.add("R", 23);//1<<6);
			params.add("Initial magnetization", 0.0);
			params.addm("T", 0.04);
			params.addm("J", -1.0);
			params.addm("h_i", 0.8);
			params.addm("h_f", 0.0);
			params.addm("dt", 1/(double)(1<<3));
			params.addm("maxTime", 20.0);
			params.addm("Bin Size", 1.0);
//			params.addm("sfLabel", 2);
			params.add("h", 0.0);
			params.add("time");
			params.add("magnetization");
			params.add("Lp");
			params.add("Reps");
			flags.add("Write Config");
		}	

		public void animate() {
			sftPlot.setLogScale(false, true);
			sftParaPlot.setLogScale(false, true);
			sftAvePlot.setLogScale(false, true);
			sftAveParaPlot.setLogScale(false, true);

			sftPlot.setAutoScale(true);
			sftParaPlot.setAutoScale(true);
			sftAvePlot.setAutoScale(true);
			sftAveParaPlot.setAutoScale(true);
			compare.setAutoScale(true);
			
			params.set("time", format(sim.time()));
			params.set("magnetization", format(sim.magnetization()));
			sim.setParameters(params);
			params.set("Lp", sim.L/dx);
			grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
			compare.registerLines("Horizontal", tempHor[0], Color.BLACK);
			compare.registerLines("Vertical", tempVer[0], Color.RED);
			histogramPlot.registerPoints("distribution", hist, Color.blue);
			normHistPlot.registerPoints("distribution", normHist, Color.blue);
			for (int i = 0; i < accNo; i ++){
				StringBuffer sb = new StringBuffer();sb.append("s(t) Ave "); sb.append(i);
				float colorChunk = (float)i/(float)accNo;
				Color col = Color.getHSBColor(colorChunk, 1.0f, 1.0f);
				sftPlot.registerLines(sb.toString(), perp_sf_tAcc[i], col);
				sb.append("p");	
				sftParaPlot.registerLines(sb.toString(), para_sf_tAcc[i], col);
				sb.append("_ave");				
				sftAvePlot.registerLines(sb.toString(), perp_sf_tAveAcc[i], col);
				sb.append("p");
				sftAveParaPlot.registerLines(sb.toString(), para_sf_tAveAcc[i], col);
				//sb = new StringBuffer();sb.append("s(t) "); sb.append(i);
				//sftPlot.registerLines(sb.toString(), sf_tAcc[i], col);
			}
			//sftPlot.registerLines("Blah", sf_tAcc[0], Color.BLACK);
			
			if(flags.contains("Write Config")) writeConfigToFile();
			flags.clear();
		}

		public void clear() {
		}

		public void run() {
			initialize();
			fft = new FourierTransformer((int)(sim.L/dx));
			sFactor = new double [sim.L/dx*sim.L/dx];
			double step = 0.20;
			double initTime = 30.0;
			double maxTime = params.fget("maxTime");//max time after quench time
			double quenchH = params.fget("h_f");// ext field after quench (initialized at h=0)
			double initH = params.fget("h_i");
			int repNo = 0;
			int sfLabel = findBestkR();
//			int sfLabel = params.iget("sfLabel");
			System.out.println(sfLabel);

			while (true) {
				sim.restartClock();
				initializeClumps(initTime, initH, quenchH);
				sim.restartClock();
				sFactor = fft.calculate2DSF(sim.getField(dx), false, false);
				sim.restartClock();
				int recordInt = 0;
				int recordStep = 0;
				step = 0.25;
				boolean clumpy = true;
				while (sim.time() < maxTime){
					sim.step();
					Job.animate();
					if(clumpy){
						//check to see if either direction of SF has reduced to 
						//half its original value
						if(sFactor[sfLabel] < 10.0){
							normHist.accum(sim.time());
							hist.accum(sim.time());
							clumpy = false;
							System.out.println("half");
						}else if(sFactor[sfLabel*sim.L] < 10.0){
							normHist.accum(sim.time());
							hist.accum(sim.time());
							clumpy = false;
							System.out.println("half");
						}
					}
					if (sim.time() > recordStep){
						sFactor = fft.find2DSF(sim.getField(dx), sim.L);
						collectTemp(sFactor, sfLabel);
						recordStep += step;
						recordInt +=1;
					}
				}	
				//after max time, check to see if stripes are horizontal or vertical
				boolean vertStripes;
				if(sFactor[sfLabel]>sFactor[sfLabel*sim.L/dx])vertStripes = true;
				else vertStripes = false;
				for (int i = 0; i < accNo; i ++){
					perp_sf_tAcc[i].clear();
					para_sf_tAcc[i].clear();
				}
				accumCollectedData(vertStripes);
				for (int i = 0; i < accNo; i ++){
					tempHor[i].clear();
					tempVer[i].clear();
				}
				repNo += 1;
				Job.animate();

				params.set("Reps", repNo);
				writeStSCtoFile(sfLabel, initTime, kRvalues(vertStripes, sfLabel));
			}
		}

		private void collectTemp(double [] sFactor, int sfLabelHor){
			int sfLabelVert = sfLabelHor*sim.L/dx;
			for(int i = 0; i < accNo; i ++){
				//vertical data
				sfLabel [i] = sfLabelVert + i;
				tempVer [i].accum(sim.time(), sFactor[sfLabel[i]]);
				//horizontal data
				sfLabel [i] = sfLabelHor +  i*sim.L/dx;
				tempHor [i].accum(sim.time(), sFactor[sfLabel[i]]); 
			}
		}
		
		private void accumCollectedData(boolean verticalStripes){
			System.out.println("Vertical stripes " + verticalStripes);
			for (int i = 0; i < accNo; i ++){
				DatasetBuffer dataVert = tempVer[i].copyData();
				int vertSize = dataVert.size();
				DatasetBuffer dataHor = tempHor[i].copyData();
				int horSize = dataHor.size();
				if(verticalStripes){
				// if stripes are vertical, parallel direction is vert, perp direction is hor
					for(int j = 0; j < vertSize; j++){
						para_sf_tAveAcc[i].accum(dataVert.x(j), dataVert.y(j));
						para_sf_tAcc[i].accum(dataVert.x(j), dataVert.y(j));
					}

					for(int j = 0; j < horSize; j++){
						perp_sf_tAveAcc[i].accum(dataHor.x(j), dataHor.y(j));
						perp_sf_tAcc[i].accum(dataHor.x(j), dataHor.y(j));
					}									
				}else{
					// if stripes are horizontal, parallel direction is hor, perp direction is vert
					for(int j = 0; j < vertSize; j++){
						perp_sf_tAveAcc[i].accum(dataVert.x(j), dataVert.y(j));
						perp_sf_tAcc[i].accum(dataVert.x(j), dataVert.y(j));
					}

					for(int j = 0; j < horSize; j++){
						para_sf_tAveAcc[i].accum(dataHor.x(j), dataHor.y(j));
						para_sf_tAcc[i].accum(dataHor.x(j), dataHor.y(j));
					}				
				}
			}
		}
		
//		private void collect(double [] sFactor, boolean vertStripes, int sfLabelHor){
//			int sfLabelVert = sfLabelHor*sim.L/dx;
//			if(vertStripes){
//				for (int i = 0; i < accNo; i ++)
//					sfLabel[i] = sfLabelVert + i;
//			}else{
//				for (int i = 0; i < accNo; i ++)
//					sfLabel[i] = sfLabelHor + i*sim.L/dx;
//			}
//			for (int i = 0; i < accNo; i ++){
//				//int kx = sfLabel[i]%sim.L/dx; int ky = sfLabel[i]/sim.L/dx;
//				//double kxR = 2*PI*kx*sim.R/sim.L;double kyR = 2*PI*ky*sim.R/sim.L;
//				//System.out.println("kxR = " + kxR + " kyR =  " + kyR);
//				sf_tAveAcc[i].accum(sim.time(),sFactor[sfLabel[i]]);			
//				sf_tAcc[i].accum(sim.time(),sFactor[sfLabel[i]]);
//			}
//		}
		
		private double [] kRvalues(boolean vertStripes, int sfLabelHor){
			double [] kRvalue = new double [accNo*2];
			int [] kRCoord = kRLabels(vertStripes, sfLabelHor);
			for (int i = 0; i < accNo; i++){
				kRvalue[2*i] = 2*PI*kRCoord[2*i]*sim.R/(sim.L);
				kRvalue[2*i+1] = 2*PI*kRCoord[2*i+1]*sim.R/(sim.L); 
			}	
			return kRvalue;
		}
			
			private int [] kRLabels(boolean vertStripes, int sfLabelHor){
				int [] kRCoord = new int [2*accNo];
				int sfLabelVert = sfLabelHor*sim.L/dx;
				if(vertStripes){
					for (int i = 0; i < accNo; i ++){
						int label = sfLabelVert + i;
						kRCoord[i*2] = label %sim.L/dx;
						kRCoord[i*2+1] = label /(sim.L/dx);
					}
				}else{
					for (int i = 0; i < accNo; i ++){
						int label = sfLabelHor + i*sim.L/dx;
						kRCoord[i*2] = label %sim.L/dx;
						kRCoord[i*2+1] = label /(sim.L/dx);					
					}
				}
				return kRCoord;
			}
		
		private void initializeClumps(double initTime, double initH, double finalH){
			params.set("h", initH);
			sim.randomizeField(params.fget("Initial magnetization"));	
			while(sim.time() < initTime){ 
				sim.step();
				Job.animate();		
			}
			params.set("h", finalH);
		}
		
		public int findBestkR(){
			int kRint = (int)(IsingLR.kRpeak*sim.L/(2*Math.PI*sim.R));
			double trialkR1 = 2*PI*kRint*sim.R/sim.L;
			double trialkR2 = 2*PI*(kRint+1)*sim.R/sim.L;
			if (abs(IsingLR.kRpeak-trialkR1) > abs(IsingLR.kRpeak-trialkR2)){
				kRint += 1;
			}		
			return kRint;
		}

		
		private void writeStSCtoFile(int sfInt, double initializeTime, double [] kRvalues){
			String message1 = "#Glauber Monte Carlo run: S vs t for several values of k. Stripe to clump H quench.";
			String fileName = params.sget("Data Dir") + File.separator + "e0";
			StringBuffer fileBuffer = new StringBuffer(); fileBuffer.append(fileName);
			for (int i=0; i < accNo; i ++){
				//System.out.println("start " + i);
				StringBuffer mb = new StringBuffer();
				//mb.append("# init time = "); mb.append(initializeTime);
				mb.append("# kx value = ");	mb.append(kRvalues[2*i]); mb.append(" ky value = ");
				//double krvalue = 2*sim.R*Math.PI*(sfInt+i)/sim.L;
				mb.append(kRvalues[2*i+1]);			
				fileBuffer.deleteCharAt(fileBuffer.length()-1);	fileBuffer.append(i); fileName = fileBuffer.toString();
				String message2 = mb.toString();
				FileUtil.initFile(fileName, params, message1, message2);		
				FileUtil.printAccumToFile(fileName, perp_sf_tAveAcc[i]);
			}
			String message11 = "#Glauber Monte Carlo run: S vs t for several values of k. Stripe to clump H quench.";
			String fileNamePara = params.sget("Data Dir") + File.separator + "a0";
			StringBuffer fileBufferPara = new StringBuffer(); fileBufferPara.append(fileNamePara);
			for (int i=0; i < accNo; i ++){
				//System.out.println("start " + i);
				StringBuffer mb = new StringBuffer();
				//mb.append("# init time = "); mb.append(initializeTime);
				mb.append("# kx value = ");	mb.append(kRvalues[2*i]); mb.append(" ky value = ");
				//double krvalue = 2*sim.R*Math.PI*(sfInt+i)/sim.L;
				mb.append(kRvalues[2*i+1]);			
				fileBufferPara.deleteCharAt(fileBufferPara.length()-1);	fileBufferPara.append(i); fileNamePara = fileBufferPara.toString();
				String message2 = mb.toString();
				FileUtil.initFile(fileNamePara, params, message11, message2);		
				FileUtil.printAccumToFile(fileNamePara, para_sf_tAveAcc[i]);
			}
			
			String histFile = params.sget("Data Dir") + File.separator + "h";
			String mes = "Histogram of transition times unnormalized";
			FileUtil.initFile(histFile, params, mes);
			FileUtil.printHistToFile(histFile, hist);
			
			String normHistFile = params.sget("Data Dir") + File.separator + "nh";
			mes = "Histogram of transition times normalized";
			FileUtil.initFile(normHistFile, params, mes);
			FileUtil.printHistToFile(normHistFile, normHist);
		}
		
		
		public void initialize(){
			for (int i = 0; i < accNo; i++){
				perp_sf_tAcc[i] = new Accumulator();
				perp_sf_tAveAcc[i] = new Accumulator(); perp_sf_tAveAcc[i].enableErrorBars(true);
				
				para_sf_tAcc[i] = new Accumulator();
				para_sf_tAveAcc[i] = new Accumulator(); para_sf_tAveAcc[i].enableErrorBars(true);
				tempHor[i] = new Accumulator();
				tempVer[i] = new Accumulator();
			}

			sim = new IsingLR(params);
			normHist = new Histogram(params.iget("Bin Size"));
			normHist.setNormalizing(true);
			hist = new Histogram(params.iget("Bin Size"));			
			sim.randomizeField(params.fget("Initial magnetization"));		
			dx = 1;
		
		}
		
		public void writeConfigToFile(){
			String configFileName = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/configs/config";
			FileUtil.deleteFile(configFileName);
			FileUtil.writeConfigToFile(configFileName, (sim.L/dx)*(sim.L/dx), sim.getField(dx));
		}
		
		public void readInitialConfiguration(){
			try{
				File myFile = new File("../../../research/javaData/stripeToClumpInvestigation/monteCarloData/configs/config");
				DataInputStream dis = new DataInputStream(new FileInputStream(myFile));
				int spaceIndex;
				double phiValue;
				int Lp = sim.L/dx;
				try{
					while(true){
						spaceIndex =dis.readInt();
						dis.readChar();       // throws out the tab
						phiValue = dis.readDouble();
						dis.readChar();
						sim.spins.set(spaceIndex%Lp,spaceIndex/Lp,(int)phiValue);
//						[spaceIndex] = phiValue;
					}
				} catch (EOFException e) {
				}

			} catch (FileNotFoundException e) {
				System.err.println("FileStreamsTest: " + e);
			} catch (Exception ex) {
				ex.printStackTrace();
			}
		}
			

	
}
