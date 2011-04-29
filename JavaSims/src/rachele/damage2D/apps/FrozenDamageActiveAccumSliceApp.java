package rachele.damage2D.apps;

import java.awt.Color;
import java.io.File;
import rachele.damage2D.FrozenDamageLatticeMin;
import rachele.util.FileUtil;
import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.dataset.PointSet;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import scikit.util.DoubleArray;



/**
 * Newly minimized 1/5/2011.
 * Removed extraneous functions.
 * This version runs OFC and records active sites for various sizes and has slices with slide function.
 * Seed site accumulation and ALL active site accumulation have been removed.
 * This version specifically made to generate data for Monograph paper due Jan 31, 2011.
 * Changed accumulation of active sites to read from a list of failed sites 
 * (rather than cycle through all sites at each plate update.)
 * Hopefully much faster.
 */
public class FrozenDamageActiveAccumSliceApp extends Simulation{
	
	int dt;
	int pow;
	int maxSize;	
	int cg_dt;
	int siteSize;

	String sizeHistFile;
	String infoFile;
		
	FrozenDamageLatticeMin ofc;
	
	Accumulator maxSizeTempAcc = new Accumulator();	
	Accumulator scaledTimeAcc = new Accumulator();
	Accumulator scaledLiveTimeAcc = new Accumulator();
	Accumulator alphaAcc = new Accumulator();
	Accumulator cmSeedDist = new Accumulator(1);
	Accumulator radGyrationAcc = new Accumulator(1);

	Histogram sizeHist = new Histogram(1);
	

	Grid deadGrid = new Grid("Lattice");
	Grid gammaGrid = new Grid("Gamma");
	
	Grid activeGrid1 = new Grid("size=1");
	Grid activeGrid10 = new Grid("1<=size<10");
	Grid activeGrid100 = new Grid("10<=size<100");
	Grid activeGrid1000 = new Grid("100<=size<1000");
	Grid activeGrid10000 = new Grid("1000<=size<10000");
	Grid activeGrid100000 = new Grid("10000<=size<100000");
	Plot activeSlice1 = new Plot("size=1");
	Plot activeSlice10 = new Plot("1<=size<10");
	Plot activeSlice100 = new Plot("10<=size<100");
	Plot activeSlice1000 = new Plot("100<=size<1000");
	Plot activeSlice10000 = new Plot("1000<=size<10000");
	Plot activeSlice100000 = new Plot("10000<=size<100000");
	Plot gammaSlicePlot = new Plot("Gamma Slice");
	int [][] largeFailSites;

	public static void main(String[] args) {
		new Control(new FrozenDamageActiveAccumSliceApp(), "OFC Damage Model Scaling Only");
	}

	public void load(Control c) {
		c.frameTogether("Data",deadGrid, gammaGrid);//, allActiveGrid, allActiveSlice, seedSitesGrid);
		c.frameTogether("Active Sites for Event Size = Pi*r^2", activeGrid1, activeGrid10, activeGrid100,activeGrid1000, activeGrid10000, activeGrid100000);
		c.frameTogether("Slices", activeSlice1 ,activeSlice10 ,activeSlice100, activeSlice1000,activeSlice10000, activeSlice100000);
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
		String da = "DeadBlocks";
		String ah = "GaussianAboutHalf";
		String gs = "GaussianSplit";
		String ga = "Gaussian";
		String fa = "Flat Random";
		String az =  "Gaussian about zero";
		params.add("Alpha Distribution", new ChoiceValue(ca, mg, ca, ei, qu , fr, mg, da, ah, gs, ga, fa, az));
		params.addm("Random Seed", 0);
		params.addm("Size Power",8);
		params.addm("R", 16);
		params.addm("Init Percent Dead", .14);
		params.addm("Dead Parameter", 16);
		params.addm("Number Dead", 1024);
		params.addm("Coarse Grained dt (PU)", 1);
		params.addm("Equilibration Updates", 50000);
		params.addm("Max PU",3000000);
		params.addm("Data points per write", 10000);
		params.addm("Residual Stress", 0.625);
		params.addm("Res. Max Noise", 0.125);
		params.addm("Dissipation Param", 0.0);
		params.addm("Damage Tolerance", 1.0);
		params.add("Av Size");
		params.add("Large Av Size");
		params.add("Plate Updates");
		params.add("Percent Damage");
		params.add("Gamma");
		params.addm("Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
	}

	public void animate() {
		params.set("Plate Updates", ofc.plateUpdates);
		params.set("Av Size", siteSize);
		params.set("Large Av Size", ofc.avSize);
		double slice = params.fget("Slice");
		gammaSlicePlot.setAutoScale(true);
		gammaSlicePlot.registerLines("Gamma Slice", getSlice(slice, ofc.gamma), Color.GREEN);
		gammaGrid.clearDrawables();
		gammaGrid.addDrawable(
				Geom2D.line(0, slice, 1, slice, Color.GREEN));
		activeGrid1.clearDrawables();
		activeGrid10.clearDrawables();
		activeGrid100.clearDrawables();
		activeGrid1000.clearDrawables();
		activeGrid10000.clearDrawables();
		activeGrid100000.clearDrawables();
		activeGrid1.registerData(ofc.L, ofc.L, largeFailSites[0]);
		activeGrid10.registerData(ofc.L, ofc.L, largeFailSites[1]);
		activeGrid100.registerData(ofc.L, ofc.L, largeFailSites[2]);
		activeGrid1000.registerData(ofc.L, ofc.L, largeFailSites[3]);
		activeGrid10000.registerData(ofc.L, ofc.L, largeFailSites[4]);
		activeGrid100000.registerData(ofc.L, ofc.L, largeFailSites[5]);
		activeGrid1.addDrawable(
				Geom2D.line(0, slice, 1, slice, Color.RED));
		activeGrid10.addDrawable(
				Geom2D.line(0, slice, 1, slice, Color.RED));
		activeGrid100.addDrawable(
				Geom2D.line(0, slice, 1, slice, Color.RED));
		activeGrid1000.addDrawable(
				Geom2D.line(0, slice, 1, slice, Color.RED));
		activeGrid10000.addDrawable(
				Geom2D.line(0, slice, 1, slice, Color.RED));
		activeGrid100000.addDrawable(
				Geom2D.line(0, slice, 1, slice, Color.RED));
		activeSlice1.setAutoScale(true);
		activeSlice10.setAutoScale(true);
		activeSlice100.setAutoScale(true);
		activeSlice1000.setAutoScale(true);
		activeSlice10000.setAutoScale(true);
		activeSlice100000.setAutoScale(true);
		activeSlice1.registerPoints("Gamma Slice", getGammaSlice(slice, ofc.gamma), Color.GREEN);
		activeSlice10.registerPoints("Gamma Slice", getGammaSlice(slice, ofc.gamma), Color.GREEN);
		activeSlice100.registerPoints("Gamma Slice", getGammaSlice(slice, ofc.gamma), Color.GREEN);
		activeSlice1000.registerPoints("Gamma Slice", getGammaSlice(slice, ofc.gamma), Color.GREEN);
		activeSlice10000.registerPoints("Gamma Slice", getGammaSlice(slice, ofc.gamma), Color.GREEN);
		activeSlice100000.registerPoints("Gamma Slice", getGammaSlice(slice, ofc.gamma), Color.GREEN);		
		activeSlice1.registerPoints("size=1", getSlice(slice, largeFailSites[0]), Color.RED);
		activeSlice10.registerPoints("1<=size<10", getSlice(slice, largeFailSites[1]), Color.RED);
		activeSlice100.registerPoints("10<=size<100", getSlice(slice, largeFailSites[2]), Color.RED);
		activeSlice1000.registerPoints("100<=size<1000", getSlice(slice, largeFailSites[3]), Color.RED);
		activeSlice10000.registerPoints("1000<=size<10000", getSlice(slice, largeFailSites[4]), Color.RED);
		activeSlice100000.registerPoints("10000<=size<100000", getSlice(slice, largeFailSites[5]), Color.RED);		
	}

	public void clear() {
	}
	
	public PointSet getSlice(double horizontalSlice, double [] latticeData){
		int y = (int) (horizontalSlice * ofc.L);
		double slice[] = new double[ofc.L];
		for (int x = 0; x < ofc.L; x++) {
			slice[x] = latticeData[ofc.L*y + x];
		}
		return new PointSet(0, 1.0, slice);
	}

	public PointSet getGammaSlice(double horizontalSlice, double [] latticeData){
		int y = (int) (horizontalSlice * ofc.L);
		double slice[] = new double[ofc.L];
		for (int x = 0; x < ofc.L; x++) {
			slice[x] = 1.0-latticeData[ofc.L*y + x];
		}
		return new PointSet(0, 1.0, slice);
	}
	
	public PointSet getSlice(double horizontalSlice, int [] latticeDataInt){
		int y = (int) (horizontalSlice * ofc.L);
		double slice[] = new double[ofc.L];
		double max=0;
		for (int i=0; i<latticeDataInt.length; i++) if(latticeDataInt[i]>max) max=latticeDataInt[i];
		double [] latticeData = new double [latticeDataInt.length];
		for (int i=0; i<latticeDataInt.length; i++) latticeData[i] = (double)latticeDataInt[i]/max;
		for (int x = 0; x < ofc.L; x++) {
			slice[x] = latticeData[ofc.L*y + x];
		}
		return new PointSet(0, 1.0, slice);
	}

	
	public void run() {
		infoFile = params.sget("Data Dir") + File.separator + "info.txt";
		String alphaHistFile = params.sget("Data Dir") + File.separator + "ah.txt";
		ofc = new FrozenDamageLatticeMin(params, infoFile);
		ofc.initLattice(params);
		FileUtil.initFile(alphaHistFile, params);
		alphaAcc.enableErrorBars(true);
		cmSeedDist.enableErrorBars(true);
		radGyrationAcc.enableErrorBars(true);
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
		int shCount = 1;

		while(true){
			ofc.ofcStep();
			int size = ofc.avSize;
			sizeHist.accum(size);
			if (size > maxSize) {
				maxSize = size;
			}
			siteSize=size;

			for (int st = 0; st < size; st++){
				for (int r = 0; r<=5; r++){
					double sz = Math.pow(10, r-1);
					double sz1 = Math.pow(10,r);
					if((double)size > sz & size <= sz1) largeFailSites[r][ofc.failedSitesList[st]]+=1;
				}
			}

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
				Job.animate();
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
	}
	
	void allocate(){
		largeFailSites = new int [6][ofc.N];
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
//			deadSites[trackSite]=0.5;
			deadGrid.registerData(Lp, Lp, deadSites);
			for (int i = 0; i < Np; i++){
				if(ofc.aliveLattice[i]) deadSites[i]=1.0 - ofc.gamma[i];
				else deadSites[i]=0;
			}
//			deadSites[trackSite]=0.5;
			gammaGrid.setScale(0.0, 1.0);
			gammaGrid.registerData(Lp, Lp, deadSites);
			FileUtil.printlnToFile(infoFile, "# Percent alive for dx=" + dx + " is " + percentAlive);
			double perDam = 1.0-percentAlive;
			System.out.println("Percent Damage = " + perDam);
			params.set("Percent Damage", perDam);
			params.set("Gamma", ofc.aveGamma);
			
	}
	
	void writeAccumulatedData(){
		FileUtil.initFile(sizeHistFile, params, " Avalanch Size Histogram Data File");
		FileUtil.printHistToFile(sizeHistFile, sizeHist);
		maxSize = 0;
		//clear temp accumulators
		scaledTimeAcc.clear();
		scaledLiveTimeAcc.clear();
		maxSizeTempAcc.clear();
	}
	
	void initFiles(){
		String parentDir = params.sget("Data Dir") + File.separator;
		sizeHistFile = parentDir + "sh.txt";		
	}
}
