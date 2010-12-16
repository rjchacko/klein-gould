package rachele.damage2D.apps;

import static scikit.util.Utilities.asList;

import java.awt.Color;
import java.io.File;

import rachele.damage2D.multidx.FrozenDamageLattice;
import rachele.util.FileUtil;
import scikit.dataset.Accumulator;
import scikit.dataset.DatasetBuffer;
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
import scikit.util.Utilities;

public class FrozenDamageSeedActiveAccum extends Simulation{
	
//	private static final String RED = null;
	int dt;
	int pow;
	int maxSize;	
	int cg_dt;
	int siteSize;
	int trackSite;
	
	String sizeHistFile;
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
	Accumulator radGyrationAcc = new Accumulator(1);


	Grid deadGrid = new Grid("Lattice");
	Grid gammaGrid = new Grid("Gamma");
//	Grid sizeGrid1 = new Grid("1<=size<10");
//	Grid sizeGrid10 = new Grid("10<=size<100");
//	Grid sizeGrid100 = new Grid("100<=size<1000");
//	Grid sizeGrid1000 = new Grid("1000<=size<10000");
//	Grid sizeGrid10000= new Grid("10000<=size<100000");
//	Grid sizeGrid100000 = new Grid("100000<=size<1000000");
	
	Grid activeGrid1 = new Grid("size=1");
	Grid activeGrid10 = new Grid("1<=size<10");
	Grid activeGrid100 = new Grid("10<=size<100");
	Grid activeGrid1000 = new Grid("100<=size<1000");
	Grid activeGrid10000 = new Grid("1000<=size<10000");
	Grid activeGrid100000 = new Grid("10000<=size<100000");
	Grid accumOneSite = new Grid("Seed site (red), CM (green)");
	
	Plot activeSlice1 = new Plot("size=1");
	Plot activeSlice10 = new Plot("1<=size<10");
	Plot activeSlice100 = new Plot("10<=size<100");
	Plot activeSlice1000 = new Plot("100<=size<1000");
	Plot activeSlice10000 = new Plot("1000<=size<10000");
	Plot activeSlice100000 = new Plot("10000<=size<100000");
	Plot radGyrationPlot = new Plot("Rad of Gyration vs Size");
	
	Plot distPlot = new Plot("CM-Seed Distance");
	Plot gammaSlicePlot = new Plot("Gamma Slice");
//	Plot activeSlicesPlot = new Plot("Active Slices");
	
	int [][] sizeCount;
	int [] oneSite;
	int [][] largeFailSites;
	double [] cm = new double [4];
//	int [] eventCount;

	public static void main(String[] args) {
		new Control(new FrozenDamageSeedActiveAccum(), "OFC Damage Model Scaling Only");
	}

	public void load(Control c) {
		c.frameTogether("Data",deadGrid, gammaGrid,distPlot);
//		c.frameTogether("Seed Sites for Event Size = Pi*r^2", sizeGrid1, sizeGrid10, sizeGrid100,sizeGrid1000, sizeGrid10000, sizeGrid100000);
		c.frameTogether("Active Sites for Event Size = Pi*r^2", activeGrid1, activeGrid10, activeGrid100,activeGrid1000, activeGrid10000, activeGrid100000);
		c.frame(accumOneSite);
		c.frame(radGyrationPlot);
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
		params.addm("Size Power",7);
		params.addm("R", 16);
		params.addm("Init Percent Dead", .05);
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
		params.addm("Slice", new DoubleValue(0.5, 0, 0.9999).withSlider());
	}

	public void animate() {
		params.set("Plate Updates", ofc.plateUpdates);
		params.set("Av Size", siteSize);
		params.set("Large Av Size", ofc.avSize);
//		sizeGrid1.registerData(ofc.L, ofc.L, sizeCount[0]);
//		sizeGrid10.registerData(ofc.L, ofc.L, sizeCount[1]);
//		sizeGrid100.registerData(ofc.L, ofc.L, sizeCount[2]);
//		sizeGrid1000.registerData(ofc.L, ofc.L, sizeCount[3]);
//		sizeGrid10000.registerData(ofc.L, ofc.L, sizeCount[4]);
//		sizeGrid100000.registerData(ofc.L, ofc.L, sizeCount[5]);
		accumOneSite.setAutoScale(true);
		accumOneSite.registerData(ofc.L, ofc.L, oneSite);

//		distPlot.registerLines("Distance CM-Seed", cmSeedDist, Color.RED);

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
	
		radGyrationPlot.registerPoints("Rad of Gyration", radGyrationAcc, Color.BLUE);
		
		double x = (double)(trackSite%ofc.L)/(double)ofc.L;
		double y = (double)(trackSite/ofc.L)/(double)ofc.L;
		accumOneSite.setDrawables(asList(
				Geom2D.circle(x, y, 1.0/64.0,Color.RED),
//				Geom2D.circle(cm[0], cm[1], 1.0/64.0,Color.GREEN),
				Geom2D.circle(cm[0], cm[1], cm[3]/ofc.L,Color.GREEN)));
		
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
		ofc = new FrozenDamageLattice(params, infoFile);
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
			trackSite = ofc.epicenterSite;
//			if(ofc.avSize > 1){
				accumOneSite.clear();
				cm = findAvCM();
				radGyrationAcc.accum((double)size, cm[3]);
				cmSeedDist.accum(size, cm[2]);
//				System.out.println(size + " " + cm[2]);
//				FileUtil.printlnToFile(avDataFile, ofc.epicenterSite, size);
				//clear failed sites

				for(int st = 0; st <ofc.N; st++){
					oneSite[st]=ofc.failedSites[st];
					if(ofc.failedSites[st]>0){
						for (int r = 0; r<=5; r++){
							double sz = Math.pow(10, r-1);
							double sz1 = Math.pow(10,r);
							if((double)size > sz & size <= sz1) largeFailSites[r][st]+=ofc.failedSites[st];
						}
					}

				}

				for (int r = 0; r<=5; r++){
					double sz = Math.pow(10, r);
					double sz1 = Math.pow(10,r+1);
					if((double)size >= sz & size < sz1) sizeCount[r][ofc.epicenterSite]+=1;				
				}
//			}



//			if(ofc.epicenterSite == trackSite){
//				siteSize = ofc.avSize;
//				if(siteSize > 1){
//					for(int st = 0; st <ofc.N; st++){
//						oneSite[st]+=ofc.failedSites[st];
//					}
//				}
//			}



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
//			if(size>1000){
//				Job.animate();
//				Job.signalStop();
//			}
		}
	}
	
	double[] findAvCM(){
		int L=ofc.L;
		int sumdX=0;
		int sumdY=0;
		int sumdX2=0;
		int sumdY2=0;
		int xt = trackSite%L;
		int yt = trackSite/L;
//		System.out.println("new");
		for(int i = 0; i<ofc.N; i++){
			if(ofc.failedSites[i]>0){
				int x = i%L;
				int y = i/L;
				int dx = x-xt;
				int dy = y-yt;
				if(dx>L/2) dx -= L;
				else if(dx<-L/2) dx += L;
				if(dy>L/2) dy -= L;
				else if(dy<-L/2) dy +=L;
//				System.out.println("dx = " + dx);
//				System.out.println("dy = " + dy);
				sumdX+=dx;
				sumdY+=dy;
				sumdX2+=dx*dx;
				sumdY2+=dy*dy;
			}
		}
		double[] ret = new double [4];
		double dX = (double)sumdX/(double)siteSize;
		double dY = (double)sumdY/(double)siteSize;
		double dX2 = (double)sumdX2/(double)siteSize;
		double dY2 = (double)sumdY2/(double)siteSize;
//		double radGyration = (dX2+dY2-dX*dX-dY*dY);
		double radGyration = Math.sqrt(dX2+dY2-dX*dX-dY*dY);
		dX =  (double)((dX+xt+L)%L)/(double)L;
		dY = (double)((dY+yt+L)%L)/(double)L;
		ret[0] = dX;
		ret[1] = dY;
		double dist = Math.sqrt(dX*dX+dY*dY);
		ret[2] = dist;
		ret[3] = radGyration;
		return ret;
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
