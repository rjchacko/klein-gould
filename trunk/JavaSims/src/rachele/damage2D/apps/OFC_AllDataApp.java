package rachele.damage2D.apps;

import java.io.File;
import rachele.damage2D.OFC_Lattice;
import rachele.util.FileUtil;
import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.util.Utilities;

public class OFC_AllDataApp extends Simulation{
	int cg_dt;
	
	Grid grid = new Grid("Lattice");
	Grid cgGrid = new Grid(" CG grid");
	OFC_Lattice ofc;
	Accumulator cgMetricAcc = new Accumulator();
	Accumulator sizeStore = new Accumulator();  //To store data at each plate update within a bracket
	Accumulator stressMetStore = new Accumulator();
	
	Histogram sizeHist = new Histogram(1);
	String iMetFile;
	String sizeAllFile;
	String iMetAllFile;
	String sizeHistFile;
	int fileCount = 1;
	int maxSize = 0;
	int maxUpdatesPerFile = 100000;  //Keep files small enough to be managable.
	
	public static void main(String[] args) {
		new Control(new OFC_AllDataApp(), "OFC Model");
	}
	
	public void load(Control c) {
		c.frameTogether("Grids", grid, cgGrid);
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/testRuns"));
		params.addm("Random Seed", 1);
		params.addm("CG size", 32);
		params.addm("dx", 8);
		params.addm("Coarse Grained dt", 500);
		params.addm("Equilibration Updates", 1000000);
		params.addm("Max Time", 1000000);
		params.addm("R", 16);// 0 -> fully connected
		params.addm("Residual Stress", 0.625);
		params.addm("Dissipation Param", 0.05);
		params.addm("Res. Max Noise", 0.125);
		params.addm("Lower Cutoff", 1);
		params.add("L");
		params.add("Time");
		params.add("Plate Updates");
	}
	
	public void animate() {
		grid.registerData(ofc.L, ofc.L, ofc.stress);
		cgGrid.registerData(ofc.Lp, ofc.Lp, ofc.epicenterCount);
		params.set("Time", Utilities.format(ofc.cg_time));
		params.set("Plate Updates", ofc.plateUpdates);
	}

	public void clear() {
	}

	public void run() {
		ofc = new OFC_Lattice(params);
		initFiles();
		cg_dt = params.iget("Coarse Grained dt");
		double nextRecordTime = cg_dt;
		int maxTime = params.iget("Max Time");
		int nextMakeNewFileTime = maxUpdatesPerFile;
		
		//equilibrate
		ofc.initEquilibrate(params.iget("Equilibration Updates"));
		while (ofc.plateUpdates < 0){
			ofc.equilibrate();
			Job.animate();
		}
		
		while(true){

			ofc.step();
			int size = ofc.avSize;
			sizeHist.accum(size);
			sizeStore.accum(ofc.cg_time, size);
			if (size > maxSize) maxSize = size;
			double iMet = ofc.calcInverseMetric();
			stressMetStore.accum(ofc.cg_time, iMet);

			if(ofc.cg_time > nextRecordTime){


				FileUtil.initFile(sizeHistFile, params, "avalanch size histogram");
				FileUtil.printHistToFile(sizeHistFile, sizeHist);
				double cgInverseActivityMetric = 1.0/ofc.calcCG_activityMetric();
				double cgInverseStressMetric = 1.0/ofc.calcCG_stressMetric();
				double reducedTime = ofc.cg_time/cg_dt;
				FileUtil.printlnToFile(iMetFile, ofc.cg_time, iMet, reducedTime, cgInverseActivityMetric, cgInverseStressMetric, maxSize);
				FileUtil.printAccumToFileNoErrorBars(sizeAllFile, sizeStore);
				sizeStore.clear();
				FileUtil.printAccumToFileNoErrorBars(iMetAllFile, stressMetStore);
				stressMetStore.clear();
				nextRecordTime += cg_dt;
				maxSize = 0;
				if(ofc.cg_time > nextMakeNewFileTime){
					sizeAllFile = makeNewFile(sizeAllFile);
					iMetAllFile = makeNewFile(iMetAllFile);
					System.out.println(sizeAllFile + " " + iMetAllFile);
					nextMakeNewFileTime += maxUpdatesPerFile;
					fileCount += 1;
				}

			}

			if(ofc.cg_time > maxTime) Job.signalStop();
				
			Job.animate();
		}
	}
	
	String makeNewFile(String fileName){
		StringBuffer sb = new StringBuffer();
		sb.append(fileName); 
		int length = sb.length();
		for (int i = length-1; i > length - 7; i--) sb.deleteCharAt(i);
		
		System.out.println(sb.toString());
		if (fileCount > 100) Job.signalStop();
		if (fileCount< 10){
			sb.append("0"); sb.append(fileCount); sb.append(".txt");
		}else{
			sb.append(fileCount); sb.append(".txt");
		}
		String newFileName = sb.toString();
		FileUtil.initFile(newFileName, params);
		return newFileName;
	}
	
	void initFiles(){
		iMetFile = params.sget("Data Dir") + File.separator + "im.txt";  // to record iverse metric data
		FileUtil.initFile(iMetFile, params, " time (plate updates), stress inverse metric, time/coarse grained time, coarse grained activity metric, coarse grained stress metric, size of avalanche");
		sizeAllFile = params.sget("Data Dir") + File.separator + "s00.txt";   //to record size vs time data
		FileUtil.initFile(sizeAllFile, params, "avalanch size vs time");
		iMetAllFile = params.sget("Data Dir") + File.separator + "m00.txt";
		FileUtil.initFile(iMetAllFile, params, "inverse Stress metric for all times");
		sizeHistFile = params.sget("Data Dir") + File.separator + "sh.txt";	//to record size histogram data
	}


}
