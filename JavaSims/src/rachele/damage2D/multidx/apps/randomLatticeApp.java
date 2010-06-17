package rachele.damage2D.multidx.apps;

import java.io.File;

import kip.util.Random;
import rachele.damage2D.multidx.MetricCalculator;
import rachele.util.FileUtil;
import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;

public class randomLatticeApp extends Simulation{
	
	int plateUpdates;
	Random random = new Random();
	MetricCalculator mc;
	
	int [] activity;
	int [] sizeActivity;
	int pow;
	int N;
	boolean [] alive;
	
	String imFileNd;
	String imFileWd;
	String sizeHistFile;
	String infoFile;
	String [][] cgFiles;
	
	Accumulator [][] cgMetsAcc;
	Histogram sizeHist = new Histogram(1);

	public static void main(String[] args) {
		new Control(new randomLatticeApp(), "Random Lattice App");
	}
	

	public void animate() {
		params.set("Plate Updates", plateUpdates);
	}


	public void clear() {

		
	}

	public void load(Control c) {
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/contract2/testRuns/"));
		params.addm("Random Seed", 5);
		params.addm("Size Power",7);
		params.addm("Coarse Grained dt (PU)", 16384);
		params.addm("Max PU", 1638400);
		params.addm("Data points per write", 1);
		params.add("Plate Updates");
		
	}

	public void run() {
		int dt = params.iget("Coarse Grained dt (PU)");
		pow = params.iget("Size Power");
		int maxTime = params.iget("Max PU");
		int dataPointsPerWrite = params.iget("Data points per write");
		int dataPointsCount = 0;
		
		double nextAccumTime = dt;
		
		int L = (int) Math.pow(2, pow);
		N = L*L;
		activity = new int [N];
		sizeActivity = new int [N];
		alive = new boolean [N];
		for (int i = 0; i < N; i++){
			activity[i] = 0;	
			sizeActivity[i] = 0;
			alive[i] = true;
		}
		mc = new MetricCalculator(pow, alive, 0.0);
		allocate();
		initFiles();
		
		while(true){
			plateUpdates += 1;
			int site = (int)(random.nextDouble()*N);
			activity[site] += 1;
			int size = (int)(random.nextDouble()*N);
			sizeActivity[site] += size;
			sizeHist.accum(size);

			if(plateUpdates > nextAccumTime){ //Accumulate data to be written
				accumData();
				nextAccumTime += dt;
				dataPointsCount += 1;
			}
			if(dataPointsCount > dataPointsPerWrite -1){
				writeAccumulatedData();
				dataPointsCount = 0;
			}
			maxTime = params.iget("Max PU");
			if(plateUpdates > maxTime) Job.signalStop();
		}
		
	}
	
	
	void writeAccumulatedData(){
		FileUtil.initFile(sizeHistFile, params, " Avalanch Size Histogram Data File");
		FileUtil.printHistToFile(sizeHistFile, sizeHist);

		for (int i = 0; i < pow; i++){
			for (int j = 0; j < 4; j++){
				FileUtil.printAccumToFileNoErrorBars(cgFiles[i][j], cgMetsAcc[i][j]);
			}
		}
		//clear temp accumulators
		for (int i = 0; i < pow; i++){
			for (int j = 0; j < 4; j++){
				cgMetsAcc[i][j].clear();
			}
		}
	}
	
	void accumData(){
		// returns an array of size [pow][6]
		mc.calcNewDxAveArrays(plateUpdates, activity, sizeActivity);
		double [][] cgMets = mc.calcCG_activityMetrics();
		for (int i = 0; i < pow; i++){
			for (int j = 0; j < 4; j++){
				cgMetsAcc[i][j].accum(plateUpdates, 1.0/cgMets[i][j]);
			}
		}
		for (int i = 0; i < N; i++){
			activity[i] = 0;
			sizeActivity[i] = 0;
		}
	}
	
	void allocate(){
		cgMetsAcc = new Accumulator[pow][4];
		for (int i = 0; i < pow; i++){
			for (int j = 0; j < 4; j ++){
				cgMetsAcc[i][j] = new Accumulator();
			}
		}
		cgFiles = new String [pow][4];

	}
	
	void initFiles(){
		String parentDir = params.sget("Data Dir") + File.separator;
		imFileWd = parentDir + "imWd.txt";
		FileUtil.initFile(imFileWd, params, "Stress Met not coarse grained in time, dead sites included");
		imFileNd = parentDir + "imNd.txt";
		FileUtil.initFile(imFileNd, params, "Stress Met not coarse grained in time, no dead sites included");		

		for(int i = 0; i < pow; i++){
			int dx = (int)Math.pow(2,i);
			String dirName_wd = parentDir + "wd" + File.separator + "dx" + dx + File.separator;
			FileUtil.makeDirs(dirName_wd);
			String dirName_nd = parentDir + "nd" + File.separator + "dx" + dx + File.separator;
			FileUtil.makeDirs(dirName_nd);
			String im = " Inverse Coarse-grained Metric File: ";
			String nd = " No dead sites included in metric calculation. ";
			String wd = " Dead sites included in metric calculation. ";
			String time = " time (Plate Updates), ";
			String imm = " inverse cg metric.";
			cgFiles[i][0] = dirName_nd + "m.txt"; 
			FileUtil.initFile(cgFiles[i][0], params, im + nd + time + " stress" + imm);
			cgFiles[i][1] = dirName_wd + "m.txt";
			FileUtil.initFile(cgFiles[i][1], params, im + wd + time + " stress" + imm);
			cgFiles[i][2] = dirName_nd + "a.txt";
			FileUtil.initFile(cgFiles[i][2], params, im + nd + time + " activity" + imm);
			cgFiles[i][3] = dirName_wd + "a.txt";
			FileUtil.initFile(cgFiles[i][3], params, im + wd + time + " activity" + imm);
			cgFiles[i][4] = dirName_nd + "s.txt";
			FileUtil.initFile(cgFiles[i][4], params, im + nd + time + " size-activity" + imm);
			cgFiles[i][5] = dirName_wd + "s.txt";
			FileUtil.initFile(cgFiles[i][5], params, im + wd + time + " size-activity" + imm);
		}
		sizeHistFile = parentDir + "sh.txt";
	}

}
