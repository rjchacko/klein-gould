package rachele.damage2D.apps;


import java.io.File;
import rachele.damage2D.multidx.FrozenDamageLattice;
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

/**
 * Code is a downsizing of frozenDamageMultidxApp
 * No metric calcs are done.
 */
public class FrozenDamageScalingOnlyApp extends Simulation{
	
	int dt;
	int pow;
	int maxSize;	
	int cg_dt;
	int siteSize;
	int trackSite;
	
	String sizeHistFile;
	String alphaHistFile;
	String sizeFile;
	String infoFile;
	String avDataFile;
	double alphaDissRate;
		
	FrozenDamageLattice ofc;
		
	Histogram sizeHist = new Histogram(1);
	
	Accumulator maxSizeTempAcc = new Accumulator();	
	Accumulator scaledTimeAcc = new Accumulator();
	Accumulator scaledLiveTimeAcc = new Accumulator();
	Accumulator alphaAcc = new Accumulator();
	Accumulator cmSeedDist = new Accumulator(1);


	Grid deadGrid = new Grid("Lattice");
	Grid gammaGrid = new Grid("Alpha Prime");

	
	int [][] sizeCount;
	int [] oneSite;
	int [][] largeFailSites;
	double [] cm = new double [3];
//	int [] eventCount;

	public static void main(String[] args) {
		new Control(new FrozenDamageScalingOnlyApp(), "OFC Damage Model Scaling Only");
	}

	public void load(Control c) {
		c.frameTogether("Data",deadGrid, gammaGrid);
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/contract2/testRuns"));
		String rd = "Random";
		String br = "RandomBlocks";
		String ds = "DeadStrip";
		String pr = "PlaceRandomDead";
		String db = "DeadBlock";
		String cs = "Cascade";
		String dr = "DeadRectangle";
		String cr = "CascadeRandom";
		String bl = "DeadBlocks";
		String pd = "PlaceDeadBlocks";
		String ac = "AliveCascade";
		params.add("Type of Damage", new ChoiceValue(cr, rd, ac, pd, bl, cs, br, ds, pr, cr, db, dr));
		params.add("Dead dissipation?", new ChoiceValue("Yes", "No") );
		params.add("Boundary Conditions", new ChoiceValue("Periodic", "Open"));
		String ca = "Constant";
		String ei = "Eights";
		String qu = "Quarters";
		String fr = "Fractal";
		String mg = "ManyGaussians";
		String da = "DeadBlocks";//Set to work for L=2**8
		String ah = "GaussianAboutHalf";
		String gs = "GaussianSplit";
		String ga = "Gaussian";
		String fa = "Flat Random";
		String az =  "Gaussian about zero";
		params.add("Alpha Distribution", new ChoiceValue(ca, mg, ca, ei, qu , fr, mg, da, ah, gs, ga, fa, az));
		params.addm("Random Seed", 0);
		params.addm("Size Power",6);
		params.addm("R", 16);
		params.addm("Init Percent Dead", .2);
		params.addm("Dead Parameter", 150);
		params.addm("Number Dead", 1024);
		params.addm("Coarse Grained dt (PU)", 1);
		params.addm("Equilibration Updates", 100);
		params.addm("Max PU",1000000);
		params.addm("Data points per write", 10000);
		params.addm("Residual Stress", 0.625);
		params.addm("Res. Max Noise", 0.125);
		params.addm("Dissipation Param", 0.0);
		params.addm("Damage Tolerance", 1.0);
		params.add("Av Size");
		params.add("Large Av Size");
		params.add("Plate Updates");
		params.add("Percent Damage");
		params.add("Effective Alpha");
		params.add("Error");
	}

	public void animate() {
		params.set("Plate Updates", ofc.plateUpdates);
		params.set("Av Size", siteSize);
		params.set("Large Av Size", ofc.avSize);
	
	}

	public void clear() {
	}

	
	public void run() {
		infoFile = params.sget("Data Dir") + File.separator + "info.txt";
		String alphaHistFile = params.sget("Data Dir") + File.separator + "ah.txt";
		ofc = new FrozenDamageLattice(params, infoFile);
		ofc.initLattice(params);
		FileUtil.initFile(alphaHistFile, params);
		alphaAcc.enableErrorBars(true);
		cmSeedDist.enableErrorBars(true);
//		trackSite = 6165;
		pow = params.iget("Size Power");
		allocate();
		initFiles();
		drawLattices();
//		FileUtil.printHistToFile(alphaHistFile, ofc.ad.alphaHist);
		
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
			siteSize=size;

			if(ofc.plateUpdates > nextAccumTime){ //Accumulate data to be written
				accumData();
				nextAccumTime += dt;
				dataPointsCount += 1;
			}
			if(dataPointsCount > dataPointsPerWrite -1){
				writeAccumulatedData();
				dataPointsCount = 0;
			}
			if(ofc.plateUpdates%100000==0){
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
		ofc.clearCounts();
	}
	
	void allocate(){
//		eventCount = new int [ofc.N];
		sizeCount = new int [10][ofc.N];	
		oneSite = new int [ofc.N];
		largeFailSites = new int [10][ofc.N];
	}
	
	void drawLattices(){

			int dx = 1;
			int Lp = ofc.L/dx;
			int Np = Lp*Lp;
			double [] deadSites = new double [Np];
			for (int i = 0; i < Np; i++){
				if(ofc.aliveLattice[i]) deadSites[i]=1.0;
				else deadSites[i]=0;
			}
			double percentAlive = DoubleArray.sum(deadSites)/(double)Np;
			deadGrid.setScale(0.0, 1.0);
			deadSites[trackSite]=0.5;
			deadGrid.registerData(Lp, Lp, deadSites);
			for (int i = 0; i < Np; i++){
				if(ofc.aliveLattice[i]) deadSites[i]=1.0 - ofc.gamma[i];
				else deadSites[i]=0;
			}
			deadSites[trackSite]=0.5;
			gammaGrid.setScale(0.0, 1.0);
			gammaGrid.registerData(Lp, Lp, deadSites);
			FileUtil.printlnToFile(infoFile, "# Percent alive for dx=" + dx + " is " + percentAlive);
			
	}
	
	void writeAccumulatedData(){
		FileUtil.initFile(sizeHistFile, params, " Avalanch Size Histogram Data File");
		FileUtil.printHistToFile(sizeHistFile, sizeHist);
		FileUtil.printAccumsToFile(sizeFile, maxSizeTempAcc, scaledTimeAcc, scaledLiveTimeAcc);
		maxSize = 0;

		
		//clear temp accumulators
		scaledTimeAcc.clear();
		scaledLiveTimeAcc.clear();
		maxSizeTempAcc.clear();

		DatasetBuffer data = alphaAcc.copyData();
		FileUtil.printlnToFile(infoFile, ofc.cg_time, data.y(0), data.errorY(0), data.y(1), data.errorY(1), data.y(2), data.errorY(2));
	}
	

	void initFiles(){
		String parentDir = params.sget("Data Dir") + File.separator;
		sizeFile = parentDir +"ms.txt";
		avDataFile = parentDir +"ad.txt";
		FileUtil.initFile(avDataFile, params, " Large avalanches:", " site index, av size");
		FileUtil.initFile(sizeFile, params, " Max Avalanche Size Data File", " time (plate updates),  max avalanche size");
		sizeHistFile = parentDir + "sh.txt";
	}
}
