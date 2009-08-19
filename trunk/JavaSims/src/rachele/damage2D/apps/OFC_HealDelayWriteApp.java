package rachele.damage2D.apps;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import rachele.damage2D.OFC_DamageLattice;
import rachele.util.FileUtil;
import scikit.dataset.Accumulator;
import scikit.dataset.DatasetBuffer;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.util.Utilities;

public class OFC_HealDelayWriteApp extends Simulation{

	int dt;
	Grid grid = new Grid("Lattice");
	Grid cgGrid = new Grid("CG grid");
	Grid cgGridTimeAverage = new Grid("Time ave CG grid");
	Grid plateUpdateGrid = new Grid("Plate Update grid");
	Grid deadGrid = new Grid("Alive/Dead Lattice");

	Histogram sizeHist = new Histogram(1);
	String iMetFile;
	String sizeHistFile;
	String sizeFile;
	String damageFile;
	int maxSize;	
	int cg_dt;
	double percentDead;
	double maxPercentDead;
	double maxPercentDamage;
	
	OFC_DamageLattice ofc;
	
	Accumulator invStressMetTempAcc = new Accumulator();
	Accumulator cgInvStressMetTempAcc = new Accumulator();
	Accumulator cgInvActivityMetTempAcc = new Accumulator();
	Accumulator maxSizeTempAcc = new Accumulator();
	Accumulator invSizeActMetTempAcc = new Accumulator();
	Accumulator	cgMaxSizeActBoxTempAcc = new Accumulator();
	Accumulator timeAveForMaxCGSizeActTempAcc = new Accumulator();
	Accumulator damagePercentTempAcc = new Accumulator();
	
	public static void main(String[] args) {
		new Control(new OFC_HealDelayWriteApp(), "OFC Damage-Healing Model");
	}

	public void load(Control c) {
		c.frameTogether("Grids", grid, cgGrid, plateUpdateGrid, cgGridTimeAverage);
		c.frame(deadGrid);
		
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/testRuns"));
		params.add("Interaction", new ChoiceValue( "Circle", "Fully Connected", "Square", "Small World") );
		params.addm("Random Seed", 1);
		params.addm("CG size", 30);
		params.addm("dx", 9);
		params.addm("Coarse Grained dt (PU)", 10);
		params.addm("Equilibration Updates", 500);
		params.addm("Max PU", 1000000);
		params.addm("Data points per write", 100);
		params.addm("R", 16);// 0 -> fully connected
		params.addm("Residual Stress", 0.625);
		params.addm("Dissipation Param", 0.3);
		params.addm("Res. Max Noise", 0.125);
		params.addm("Lower Cutoff", 1);
		params.addm("Mean Max Failures", 0);
		params.addm("Failures Max Noise", 0);
		params.addm("Mean Heal Time", 1);
		params.addm("Heal Time Noise", 0);
		params.addm("Init Percent Dead", .5);
//		params.addm("Max Percent Damage", 0.05);
		params.add("L");
		params.add("Time");
		params.add("Av Size");
		params.add("Plate Updates");
		params.add("Percent dead sites");
	}

	public void animate() {
		
		grid.registerData(ofc.L, ofc.L, ofc.stress);
		cgGrid.registerData(ofc.Lp, ofc.Lp, ofc.epicenterCount);
		cgGridTimeAverage.registerData(ofc.Lp, ofc.Lp, ofc.CG_ActivityTimeAve);
		plateUpdateGrid.registerData(ofc.L, ofc.L, ofc.plateUpdateFailLocations);
//		int maxNoFails = params.iget("Mean Max Failures") + params.iget("Failures Max Noise");
		deadGrid.setScale(0,1);
		double [] deadLattice = new double [ofc.N];
		for (int i = 0; i <ofc.N; i ++){
			deadLattice[i] = (ofc.noFails[i] < 0) ? 0:1;
		}
		deadGrid.registerData(ofc.L, ofc.L, deadLattice);
		System.out.println(ofc.noFails[0] + "  " + ofc.noFails[1] + " " + ofc.noFails[2]);
		
		grid.clearDrawables();
		double radius = 1.0/(2.0*ofc.L);
		double failSite_y = ((double)(ofc.epicenterSite/ofc.L))/ofc.L + radius;
		double failSite_x = ((double)(ofc.epicenterSite%ofc.L))/ofc.L + radius;
		grid.addDrawable(
				Geom2D.circle(failSite_x, failSite_y, radius, Color.GREEN));
		
		params.set("Time", Utilities.format(ofc.cg_time));
		params.set("Plate Updates", ofc.plateUpdates);

		params.set("Percent dead sites", Utilities.format(percentDead));
		params.set("Av Size", ofc.avSize);

	}

	public void clear() {
	}

	
	public void run() {
		ofc = new OFC_DamageLattice(params);
		initFiles();
		
		dt = params.iget("Coarse Grained dt (PU)");
//		maxPercentDamage = params.fget("Max Percent Damage");
		double nextAccumTime = 0;
		int maxTime = params.iget("Max PU");
		int dataPointsPerWrite = params.iget("Data points per write");
		
		maxSize = 0;
		maxPercentDead = 0;
		int dataPointsCount = 0;
		percentDead = 0;
		
		// Equilibrate in heal mode
		System.out.println("Equlibrating");
		ofc.initEquilibrate(params.iget("Equilibration Updates"));
		while (ofc.plateUpdates < 0){
//			System.out.println(ofc.fullyConnected);
			ofc.healPreStep();
			Job.animate();
			while (ofc.nextSiteToFail >= 0){
				ofc.healStepIter();
				percentDead = ofc.noDeadSites/(double)(ofc.L*ofc.L);
				Job.animate();
			}
			ofc.checkForHealing();
		}
		System.out.println("Finished Equilibration");
		ofc.cg_time = 0;
		
		while(true){
			ofc.healPreStep();
			Job.animate();
			while (ofc.nextSiteToFail >= 0){
				ofc.healStepIter();
				Job.animate();
			}	
			ofc.checkForHealing();
			ofc.cgEpicenterSizeAccum();
			
			int size = ofc.avSize;
			sizeHist.accum(size);
			if (size > maxSize) maxSize = size;
			
			percentDead = ofc.noDeadSites/(double)(ofc.L*ofc.L);
			if(percentDead > maxPercentDead) maxPercentDead = percentDead;
			
			double stressMetric = ofc.calcStressMetric();
			
			if(ofc.plateUpdates > nextAccumTime){ //Accumulate data to be written
				Job.animate();
				
				//maxSize
				maxSizeTempAcc.accum(ofc.cg_time, maxSize);
				damagePercentTempAcc.accum(ofc.cg_time, percentDead);
//				System.out.println(percentDead);
				
				int maxEpicenter=0;
				int loc = 0;
				for(int i = 0; i < ofc.Np; i++){
					if (ofc.epicenterSize[i] > maxEpicenter){
						maxEpicenter = ofc.epicenterSize[i];
						loc = i;
					}
				}
				cgMaxSizeActBoxTempAcc.accum(ofc.cg_time, maxEpicenter);
				timeAveForMaxCGSizeActTempAcc.accum(ofc.cg_time, ofc.CG_SizeActTimeAve[loc]);
				if(dataPointsCount ==0) System.out.println("max loc = " + maxEpicenter + " " + ofc.CG_SizeActTimeAve[loc]);
				
				//stress metric
				double inverseStressMetric = 1.0/stressMetric;
				invStressMetTempAcc.accum(ofc.plateUpdates, inverseStressMetric);
				//CG stress metric
				double cgInverseStressMetric = 1.0/ofc.calcCG_stressMetric();				
				cgInvStressMetTempAcc.accum(ofc.cg_time, cgInverseStressMetric);
				//CG activity metric
				double cgInverseActivityMetric = 1.0/ofc.calcCG_activityMetric();
				cgInvActivityMetTempAcc.accum(ofc.cg_time, cgInverseActivityMetric);
				// size activity metric
				double cgInverseSizeActMet = 1.0/ofc.calcCG_sizeActMetric();
				invSizeActMetTempAcc.accum(ofc.cg_time, cgInverseSizeActMet);
				
				nextAccumTime += dt;
				dataPointsCount += 1;
				maxSize = 0;
				maxPercentDead = 0;

			}
			
//			System.out.println("no data pts = " + dataPointsCount + " cg time = " + ofc.cg_time + " next accum time = " + nextAccumTime);
			
			if(dataPointsCount > dataPointsPerWrite -1){ //write accumulated data
				FileUtil.initFile(sizeHistFile, params, "avalanch size histogram");
				FileUtil.printHistToFile(sizeHistFile, sizeHist);
				printAccumsToFile(iMetFile, invStressMetTempAcc, cgInvStressMetTempAcc, cgInvActivityMetTempAcc, invSizeActMetTempAcc);
				
				printAccumsToFile2(sizeFile, maxSizeTempAcc,cgMaxSizeActBoxTempAcc, timeAveForMaxCGSizeActTempAcc);
			
				FileUtil.printAccumToFileNoErrorBars(damageFile, damagePercentTempAcc);
				
				//restart counts
				maxSize = 0;
				dataPointsCount = 0;
				
				//clear temp accumulators
				invStressMetTempAcc.clear();
				cgInvStressMetTempAcc.clear();
				cgInvActivityMetTempAcc.clear();
				maxSizeTempAcc.clear();
				invSizeActMetTempAcc.clear();
				cgMaxSizeActBoxTempAcc.clear();
				timeAveForMaxCGSizeActTempAcc.clear();
				damagePercentTempAcc.clear();

			}

		
			if(ofc.plateUpdates%maxTime == 0) Job.signalStop();

		}
	}
	
	static void printAccumsToFile(String fileName, Accumulator acc1, Accumulator acc2, Accumulator acc3, Accumulator acc4){
		DatasetBuffer data1 = acc1.copyData();
		DatasetBuffer data2 = acc2.copyData();
		DatasetBuffer data3 = acc3.copyData();
		DatasetBuffer data4 = acc4.copyData();
		int size = data1.size();
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int i = 0; i < size; i++){
				pw.println(data1.x(i) + " " + data1.y(i) + " " + data2.x(i) + " " + data2.y(i) + " " + data3.y(i) + " " + data4.y(i)) ;
			}
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	static void printAccumsToFile2(String fileName, Accumulator acc1, Accumulator acc2, Accumulator acc3){
		DatasetBuffer data1 = acc1.copyData();
		DatasetBuffer data2 = acc2.copyData();
		DatasetBuffer data3 = acc3.copyData();
		int size = data1.size();
		try{
			File file = new File(fileName);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int i = 0; i < size; i++){
				pw.println(data1.x(i) + " " + data1.y(i) + " " + data2.y(i) + " " + data3.y(i) ) ;
			}
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	void initFiles(){
		iMetFile = params.sget("Data Dir") + File.separator + "im.txt";  // to record iverse metric data
		FileUtil.initFile(iMetFile, params, " time (plate updates), stress inverse metric, time/coarse grained time, coarse grained activity inverse metric, coarse grained stress inverse metric, cg size act inverse metric");
		sizeFile = params.sget("Data Dir") + File.separator + "s.txt";
		FileUtil.initFile(sizeFile, params, "time (plate updates),  max avalanche size, time/coarse grained time,");
		sizeHistFile = params.sget("Data Dir") + File.separator + "sh.txt";	//to record size histogram data
		damageFile = params.sget("Data Dir") + File.separator + "d.txt";	//to record percent damage
		FileUtil.initFile(damageFile, params, "time (plate updates),  max percent damage");

	}
}
