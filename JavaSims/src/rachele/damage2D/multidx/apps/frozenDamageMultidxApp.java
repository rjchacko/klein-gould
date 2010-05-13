package rachele.damage2D.multidx.apps;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import rachele.damage2D.multidx.FrozenDamageLattice;
import rachele.damage2D.multidx.MetricCalculator;
import rachele.util.FileUtil;
import scikit.dataset.Accumulator;
import scikit.dataset.DatasetBuffer;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.util.DoubleArray;

public class frozenDamageMultidxApp extends Simulation{

	
	int dt;
	int pow;
	int maxSize;	
	int cg_dt;
	
	String imFileNd;
	String imFileWd;
	String sizeHistFile;
	String sizeFile;
	String infoFile;
	String [][] cgFiles;
	
	double stressMetric;
	boolean includeDead;
	
	FrozenDamageLattice ofc;
	MetricCalculator mc;
	
	Histogram sizeHist = new Histogram(1);
	Accumulator maxSizeTempAcc = new Accumulator();	
	Accumulator [][] cgMetsAcc;
	Accumulator [] sMetsAcc = new Accumulator[2];
	
	int noGrids = 6;
	Grid [] deadGrid = new Grid [noGrids];
	
	public static void main(String[] args) {
		new Control(new frozenDamageMultidxApp(), "OFC Damage Model multi-dx");
	}

	public void load(Control c) {
		for(int i = 0; i < noGrids; i++){
			int dx = (int)Math.pow(2,i);
			deadGrid[i] = new Grid("Dead Blocks dx = " + dx);
		}

		c.frameTogether("Grids", deadGrid);

		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/contract2/P14-DamageCompare/L256/alpha0-04/R16/"));
		String rd = "Random";
		String br = "Random Blocks";
		String ds = "Dead Strip";
		String pr = "Place Random Dead";
		String rs = "Random + Dead Strip";
		String db = "Dead Block";
		params.add("Type of Damage", new ChoiceValue(br, ds, pr, br, rd, rs, db));
		params.add("Dead dissipation?", new ChoiceValue("Yes", "No") );
		params.addm("Random Seed", 1);
		params.addm("Size Power",9);
		params.addm("R", 16);
		params.addm("Init Percent Dead", 0.0625);
		params.addm("Dead Parameter", 8);
		params.addm("Coarse Grained dt (PU)", 10000);
		params.addm("Equilibration Updates", 100000);
		params.addm("Max PU", 1000000);
		params.addm("Data points per write", 1);
		params.addm("Residual Stress", 0.625);
		params.addm("Res. Max Noise", 0.125);
		params.addm("Dissipation Param", 0.04);
		params.addm("Damage Tolerance", 0.05);
		params.add("Av Size");
		params.add("Plate Updates");
	}

	public void animate() {
		params.set("Plate Updates", ofc.plateUpdates);
		params.set("Av Size", ofc.avSize);
	}

	public void clear() {
	}

	
	public void run() {
		infoFile = params.sget("Data Dir") + File.separator + "info.txt";
		ofc = new FrozenDamageLattice(params, infoFile);
		double tolerance = params.fget("Damage Tolerance");
		mc = new MetricCalculator(ofc.pow, ofc.aliveLattice, tolerance);
		
		pow = params.iget("Size Power");
		int noGrids = 6;
		allocate(noGrids);
		initFiles();
		drawLattices();
		
		dt = params.iget("Coarse Grained dt (PU)");
		double nextAccumTime = dt;
		int maxTime = params.iget("Max PU");
		int dataPointsPerWrite = params.iget("Data points per write");
		maxSize = 0;
		int dataPointsCount = 0;
		
		System.out.println("Equlibrating");
		FileUtil.printlnToFile(infoFile, "Equilibrating");
		ofc.initEquilibrate(params.iget("Equilibration Updates"));
		while (ofc.plateUpdates < 0){
			ofc.equilibrate();
			Job.animate();
		}
		System.out.println("Finished Equilibration");
		FileUtil.printlnToFile(infoFile, "Finished Equilibration");
		
		ofc.plateUpdates = 0;
		ofc.ofcStep();
		mc.initAveArrays(ofc.cg_time, ofc.stress, ofc.epicenterCount, ofc.epicenterSize);
		
		while(true){
			ofc.ofcStep();
			Job.animate();
			int size = ofc.avSize;
			sizeHist.accum(size);
			if (size > maxSize) {
				maxSize = size;
			}
			mc.calcNewStressAveArray(ofc.cg_time, ofc.stress);

			if(ofc.plateUpdates > nextAccumTime){ //Accumulate data to be written
				accumData();
				nextAccumTime += dt;
				dataPointsCount += 1;
			}
			if(dataPointsCount > dataPointsPerWrite -1){
				writeAccumulatedData();
				dataPointsCount = 0;
			}
			maxTime = params.iget("Max PU");
			if(ofc.plateUpdates > maxTime) Job.signalStop();
		}
	}
	
	void accumData(){
		maxSizeTempAcc.accum(ofc.cg_time, maxSize);
		maxSize = 0;
		double [] mets = mc.calcStressMets();
		sMetsAcc[0].accum(ofc.cg_time, 1.0/mets[0]);
		sMetsAcc[1].accum(ofc.cg_time, 1.0/mets[1]);
		// returns an array of size [pow][6]
		mc.calcNewDxAveArrays(ofc.cg_time, ofc.stress, ofc.epicenterCount, ofc.epicenterSize);
		double [][] cgMets = mc.calcCG_metrics();
		for (int i = 0; i < pow; i++){
			for (int j = 0; j < 6; j++){
				cgMetsAcc[i][j].accum(ofc.cg_time, 1.0/cgMets[i][j]);
			}
		}
		ofc.clearCounts();
	}
	
	void allocate(int noGrids){
		sMetsAcc[0] = new Accumulator();
		sMetsAcc[1] = new Accumulator();
		cgMetsAcc = new Accumulator[pow][6];
		for (int i = 0; i < pow; i++){
			for (int j = 0; j < 6; j ++){
				cgMetsAcc[i][j] = new Accumulator();
			}
		}
		cgFiles = new String [pow][6];

	}
	
	void drawLattices(){

		for (int j = 0; j < noGrids; j++){
			int dx = (int)Math.pow(2, j);
			int Lp = ofc.L/dx;
			int Np = Lp*Lp;
			double [] deadSites = new double [Np];
			for (int i = 0; i < Np; i++){
				if(mc.alive[mc.findAliveIndex(i,j)]) deadSites[i]=1;
				else deadSites[i]=0;
			}
			double percentAlive = DoubleArray.sum(deadSites)/(double)Np;
			deadGrid[j].setScale(0.0, 1.0);
			deadGrid[j].registerData(Lp, Lp, deadSites);
			FileUtil.printlnToFile(infoFile, "Percent alive for dx=" + dx + " is " + percentAlive);
		}
	}
	
	void writeAccumulatedData(){
		FileUtil.initFile(sizeHistFile, params, " Avalanch Size Histogram Data File");
		FileUtil.printHistToFile(sizeHistFile, sizeHist);
		FileUtil.printAccumToFileNoErrorBars(sizeFile, maxSizeTempAcc);
		maxSize = 0;

		FileUtil.printAccumToFileNoErrorBars(imFileNd, sMetsAcc[0]);
		FileUtil.printAccumToFileNoErrorBars(imFileWd, sMetsAcc[1]);
		for (int i = 0; i < pow; i++){
			for (int j = 0; j < 6; j++){
				FileUtil.printAccumToFileNoErrorBars(cgFiles[i][j], cgMetsAcc[i][j]);
			}
		}
		
		//clear temp accumulators
		sMetsAcc[0].clear();
		sMetsAcc[1].clear();
		for (int i = 0; i < pow; i++){
			for (int j = 0; j < 6; j++){
				cgMetsAcc[i][j].clear();
			}
		}
		maxSizeTempAcc.clear();
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
		String parentDir = params.sget("Data Dir") + File.separator;
//		infoFile = parentDir + "info.txt";  // to record iverse metric data
//		FileUtil.initFile(infoFile, params, "Info file for multi dt runs");
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
		sizeFile = parentDir +"ms.txt";
		FileUtil.initFile(sizeFile, params, " Max Avalanche Size Data File", " time (plate updates),  max avalanche size");
		sizeHistFile = parentDir + "sh.txt";	
		
	}
	

}
