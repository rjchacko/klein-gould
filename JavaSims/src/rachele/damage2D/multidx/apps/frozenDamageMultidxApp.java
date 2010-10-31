package rachele.damage2D.multidx.apps;

import java.io.File;
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
import scikit.util.Utilities;

public class frozenDamageMultidxApp extends Simulation{

	
	int dt;
	int pow;
	int maxSize;	
	int cg_dt;
	int [] latticeEventCount;
	int [] latticeEventSize;
	
	String imFileNd;
	String imFileWd;
	String sizeHistFile;
	String sizeFile;
	String infoFile;
	String [][] cgFiles;
	
	double stressMetric;
	double alphaDissRate;
	boolean includeDead;
	
	FrozenDamageLattice ofc;
	MetricCalculator mc;
	
	Histogram sizeHist = new Histogram(1);
	Accumulator maxSizeTempAcc = new Accumulator();	
	Accumulator [][] cgMetsAcc;
	Accumulator [] sMetsAcc = new Accumulator[2];
	Accumulator scaledTimeAcc = new Accumulator();
	Accumulator scaledLiveTimeAcc = new Accumulator();
	Accumulator alphaAcc = new Accumulator();
	
	int noGrids = 6;
	Grid [] deadGrid = new Grid [noGrids];
	Grid eventNoGrid = new Grid("No of Events");
	Grid aveSizeGrid = new Grid("Ave Size of Events");
	Grid alphaGrid = new Grid("Alpha Prime");
	Grid sizeActGrid = new Grid("Size Act");
	
	public static void main(String[] args) {
		new Control(new frozenDamageMultidxApp(), "OFC Damage Model multi-dx");
	}

	public void load(Control c) {
		for(int i = 0; i < noGrids; i++){
			int dx = (int)Math.pow(2,i);
			deadGrid[i] = new Grid("Dead Blocks dx = " + dx);
		}

		c.frameTogether("Grids", deadGrid);
		c.frameTogether("Lattice Info", alphaGrid, eventNoGrid, aveSizeGrid, sizeActGrid);

		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/contract2/testruns/"));
		String rd = "Random";
		String br = "Random Blocks";
		String ds = "Dead Strip";
		String pr = "Place Random Dead";
		String db = "Dead Block";
		String cs = "Cascade";
		String dr = "Dead Rectangle";
		String cr = "Cascade Random";
		String bl = "Dead Blocks";
		String pd = "Place Dead Blocks";
		String ac = "Alive Cascade";
		params.add("Type of Damage", new ChoiceValue(rd, ac, pd, bl, cs, br, ds, pr, cr, db, dr));
		String ca = "Constant";
		String ei = "Eights";
		String qu = "Quarters";
		String fr = "Fractal";
		String mg = "Many Gaussians";
		String da = "Dead Blocks";
		String ah = "Gaussian about half";
		String gs = "Gaussian split";
		String ga = "Gaussian";
		String fa = "Flat Random";
		String az =  "Gaussian about zero";
		params.add("Alpha Distribution", new ChoiceValue(ei, ca, qu , fr, mg, da, ah, gs, ga, fa, az));
		params.add("Dead dissipation?", new ChoiceValue("Yes", "No") );
		params.add("Boundary Conditions", new ChoiceValue( "Periodic", "Open"));
		params.addm("Random Seed", 1);
		params.addm("Size Power",6);
		params.addm("R", 8);
		params.addm("Init Percent Dead", 0.0);
		params.addm("Dead Parameter", 8);
		params.addm("Number Dead", 1);
		params.addm("Coarse Grained dt (PU)", 1);
		params.addm("Equilibration Updates", 1000);
		params.addm("Max PU",100000);
		params.addm("Data points per write", 10000);
		params.addm("Residual Stress", 0.625);
		params.addm("Res. Max Noise", 0.125);
		params.addm("Dissipation Param", 0.05);
		params.addm("Damage Tolerance", 1.0);
		params.add("Av Size");
		params.add("Plate Updates");
		params.add("Percent Damage");
		params.add("Effective Alpha");
		params.add("Error");
		

	}

	public void animate() {
		params.set("Plate Updates", ofc.plateUpdates);
		params.set("Av Size", ofc.avSize);
		
		eventNoGrid.registerData(ofc.L, ofc.L, latticeEventCount);
		sizeActGrid.registerData(ofc.L, ofc.L, latticeEventSize);
		double [] aveSize = new double [ofc.N];
		for (int i = 0; i < ofc.N; i++){
			aveSize[i] = (double)(latticeEventSize[i])/(double)(latticeEventCount[i]);
		}
		aveSizeGrid.registerData(ofc.L, ofc.L, aveSize);

	}

	public void clear() {
	}

	
	public void run() {
		infoFile = params.sget("Data Dir") + File.separator + "info.txt";
		ofc = new FrozenDamageLattice(params, infoFile);
		double tolerance = params.fget("Damage Tolerance");
		mc = new MetricCalculator(ofc.pow, ofc.aliveLattice, tolerance);
		alphaAcc.enableErrorBars(true);
		
		pow = params.iget("Size Power");
		allocate();
		initFiles();
		drawLattices();
		
		dt = params.iget("Coarse Grained dt (PU)");
		double nextAccumTime = dt;
		int maxTime = params.iget("Max PU");
		int dataPointsPerWrite = params.iget("Data points per write");
		maxSize = 0;
		int dataPointsCount = 0;
		
		System.out.println("Equlibrating");
		FileUtil.printlnToFile(infoFile, "# Equilibrating");
		ofc.initEquilibrate(params.iget("Equilibration Updates"));
		while (ofc.plateUpdates < 0){
			ofc.equilibrate();
			Job.animate();
		}
		System.out.println("Finished Equilibration");
		FileUtil.printlnToFile(infoFile, "# Finished Equilibration");
		FileUtil.printlnToFile(infoFile, "#");
		FileUtil.printlnToFile(infoFile, "# Time, Effective Alpha, Error, Dead Dissipation, Error, Alpha Dissipation, Error");
		
		ofc.plateUpdates = 0;
		ofc.ofcStep();
		mc.initAveArrays(ofc.cg_time, ofc.stress, ofc.epicenterCount, ofc.epicenterSize);
		int shCount = 1;
				
		while(true){
			ofc.ofcStep();
			Job.animate();
			
			double effAlpha = ofc.deadDiss + ofc.alphaDiss;
			alphaAcc.accum(0, effAlpha);
			alphaAcc.accum(1, ofc.deadDiss);
			alphaAcc.accum(2,ofc.alphaDiss);
			DatasetBuffer data = alphaAcc.copyData();
			params.set("Effective Alpha", Utilities.format(data.y(0)));
			params.set("Error", Utilities.format(data.errorY(0)));
			
			int size = ofc.avSize;
			sizeHist.accum(size);
			if (size > maxSize) {
				maxSize = size;
			}
			latticeEventSize[ofc.epicenterSite] += size;
			latticeEventCount[ofc.epicenterSite] += 1;
			
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
			if(ofc.plateUpdates%1000==0){
				writeShFile(shCount);
				shCount += 1;
			}
			maxTime = params.iget("Max PU");
			if(ofc.plateUpdates > maxTime) Job.signalStop();
		}
	}
	
	void writeShFile(int i){
		String parentDir = params.sget("Data Dir") + File.separator;
		String newShFile = parentDir + "sh" + i + ".txt";
		FileUtil.initFile(newShFile, params, " Avalanch Size Histogram Data File");
		FileUtil.printHistToFile(newShFile, sizeHist);
	}
	
	void accumData(){
		double scaledTime = ofc.cg_time/(double)ofc.N;
		scaledTimeAcc.accum(ofc.cg_time, scaledTime);
		double scaledLiveTime = ofc.cg_time/(double)ofc.noLiveSites;
		scaledLiveTimeAcc.accum(ofc.cg_time, scaledLiveTime);
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
	
	void allocate(){
		sMetsAcc[0] = new Accumulator();
		sMetsAcc[1] = new Accumulator();
		cgMetsAcc = new Accumulator[pow][6];
		for (int i = 0; i < pow; i++){
			for (int j = 0; j < 6; j ++){
				cgMetsAcc[i][j] = new Accumulator();
			}
		}
		cgFiles = new String [pow][6];
		
		latticeEventSize = new int [ofc.N];
		latticeEventCount = new int [ofc.N];

		
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
			FileUtil.printlnToFile(infoFile, "# Percent alive for dx=" + dx + " is " + percentAlive);
		}
		
		//print gamma plot
		int Np = ofc.L*ofc.L;
		double [] deadSites = new double [Np];
		for (int i = 0; i < Np; i++){
			if(ofc.aliveLattice[i]) deadSites[i]=1.0;
			else deadSites[i]=0;
		}
		for (int i = 0; i < Np; i++){
			if(ofc.aliveLattice[i]) deadSites[i]=1.0 - ofc.alphaP[i];
			else deadSites[i]=0;
		}
		alphaGrid.setScale(0.0, 1.0);
		alphaGrid.registerData(ofc.L, ofc.L, deadSites);


	}
	
	void writeAccumulatedData(){
		FileUtil.initFile(sizeHistFile, params, " Avalanch Size Histogram Data File");
		FileUtil.printHistToFile(sizeHistFile, sizeHist);
		FileUtil.printAccumsToFile(sizeFile, maxSizeTempAcc, scaledTimeAcc, scaledLiveTimeAcc);
		maxSize = 0;

		FileUtil.printAccumsToFile(imFileNd, sMetsAcc[0], scaledTimeAcc, scaledLiveTimeAcc);
		FileUtil.printAccumsToFile(imFileWd, sMetsAcc[1], scaledTimeAcc, scaledLiveTimeAcc);
		for (int i = 0; i < pow; i++){
			for (int j = 0; j < 6; j++){
				FileUtil.printAccumsToFile(cgFiles[i][j], cgMetsAcc[i][j], scaledTimeAcc, scaledLiveTimeAcc);
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
		scaledTimeAcc.clear();
		scaledLiveTimeAcc.clear();
		maxSizeTempAcc.clear();

		DatasetBuffer data = alphaAcc.copyData();
		FileUtil.printlnToFile(infoFile, ofc.cg_time, data.y(0), data.errorY(0), data.y(1), data.errorY(1), data.y(2), data.errorY(2));
//		double effAlpha = ofc.alphaDissRate + ofc.deadDissRate;
//		double scaledTime = ofc.cg_time/(double)ofc.N;
//		double scaledLiveTime = ofc.cg_time/(double)ofc.noLiveSites;
//		FileUtil.printlnToFile(infoFile, ofc.cg_time, effAlpha, ofc.alphaDissRate, ofc.deadDissRate, scaledTime, scaledLiveTime);
//		System.out.println("Dead Diss Rate = " + ofc.deadDissRate + " Effective Alpha = " + effAlpha +" = " + ofc.alphaDissRate + " " + ofc.deadDissRate);

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
		sizeFile = parentDir +"ms.txt";
		FileUtil.initFile(sizeFile, params, " Max Avalanche Size Data File", " time (plate updates),  max avalanche size");
		sizeHistFile = parentDir + "sh.txt";	
		
	}
	

}
