package rachele.ofcNetworks.apps;

import java.awt.Color;
import java.io.File;

import rachele.ofcNetworks.SmallWorldLatticeFixedNoRewire;
import rachele.util.FileUtil;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;

public class SmallWorldAppFixedNoRewires extends Simulation{


	Grid latticeGrid = new Grid("Lattice");
	Grid connectionGrid = new Grid("Connections");
	SmallWorldLatticeFixedNoRewire ofc;
	Histogram sizeHist = new Histogram(1);
	Histogram failedSiteHist = new Histogram(1);
	Histogram cumSizeHist = new Histogram(1);
	String sizeHistFile;
	String failedSiteHistFile;
	String cumSizeHistFile;
	String infoFile;
	int shCount = 1;
	
	public static void main(String[] args) {
		new Control(new SmallWorldAppFixedNoRewires(), "Small World OFC Fixed No Rewires");
	}

	public void load(Control c) {
		c.frame(latticeGrid);
		c.frame(connectionGrid);
		params.add("Data Dir",new DirectoryValue("/Users/racheledominguez/data/RMC-Summer2013/testRuns/"));
		params.addm("Random Seed", 0);
		int Lstart = 128;
		params.addm("L",Lstart);
		params.addm("alpha", 0.001);
//		params.addm("R", 1);
		params.addm("Rewire Probability", 0.00586);
//		params.addm("Rewire Probability", 1.0);
		int equilStart = Lstart*Lstart*100;
		params.addm("Equilibration Updates", equilStart);
		params.add("Av Size");
		params.add("Plate Updates");
	}

	public void clear() {
		connectionGrid.clearDrawables();
	}
	
	public void animate() {
		params.set("Plate Updates", ofc.plateUpdates);
		params.set("Av Size", ofc.failCount);
		latticeGrid.setScale(0.0, 1.0);
		latticeGrid.registerData(ofc.L, ofc.L, ofc.stress);

	}
	
	public void run() {
		//construct a new lattice
		ofc = new SmallWorldLatticeFixedNoRewire(params);
		if(ofc.L<128){
			drawConnections();
		}
		Job.animate();
		initFiles();

		ofc.plateUpdates = -params.iget("Equilibration Updates");

		while(true){
			//equilibrate
			if(ofc.plateUpdates <= 0){
				ofc.ofcStep();
				Job.animate();
			}else{
				ofc.ofcStep();
				sizeHist.accum(ofc.failCount);
				failedSiteHist.accum(ofc.failedSiteCount);
				int count = ofc.failCount;
				while(count>0){
					cumSizeHist.accum(count);
					count--;
				}

				if(ofc.plateUpdates % 10000 == 0) writeAccumulatedData(); 
				if(ofc.plateUpdates%100000==0 && ofc.plateUpdates < 1000000){
					writeShFile(shCount);
				}else if (ofc.plateUpdates % 10000000 == 0){
					writeShFile(shCount);
				}
				Job.animate();
				if(ofc.plateUpdates >= 1000000000){
					Job.signalStop();
				}
			}
		}
	}

	void writeShFile(int i){
		String parentDir = params.sget("Data Dir") + File.separator;
		String newShFile = parentDir + "sh" + i + ".txt";
		String newFsFile = parentDir + "sh" + i + ".txt";
		String newCumShFile = parentDir + "ch" + i + ".txt";
		FileUtil.initFile(newShFile, params, " Avalanch Size Histogram Data File");
		FileUtil.printHistToFile(newShFile, sizeHist);
		FileUtil.initFile(newFsFile, params, " No of Failed sites Size Histogram Data File");
		FileUtil.printHistToFile(newFsFile, failedSiteHist);
		FileUtil.initFile(newCumShFile, params, " Cumulative Avalanch Size Histogram Data File");
		FileUtil.printHistToFile(newCumShFile, cumSizeHist);
		shCount += 1;
	}
	
	void writeAccumulatedData(){
		FileUtil.initFile(sizeHistFile, params, " Avalanch Size Histogram Data File");
		FileUtil.printHistToFile(sizeHistFile, sizeHist);
		FileUtil.initFile(cumSizeHistFile, params, " Cumulative Avalanch Size Histogram Data File");
		FileUtil.printHistToFile(cumSizeHistFile, cumSizeHist);
		FileUtil.initFile(failedSiteHistFile, params, " No of Failed sites Size Histogram Data File");
		FileUtil.printHistToFile(failedSiteHistFile, failedSiteHist);
	}
	
	void initFiles(){
		String parentDir = params.sget("Data Dir") + File.separator;
		sizeHistFile = parentDir + "sh.txt";
		failedSiteHistFile = parentDir + "fs.txt";
		cumSizeHistFile = parentDir + "ch.txt";
		infoFile = parentDir + "info.txt";
		FileUtil.printlnToFile(infoFile, "Percent Rewired = ", ofc.percentRewired);
		FileUtil.printlnToFile(infoFile, "Using app SmallWorldAppFixedNoRewires");
	}
	
	void drawConnections(){
		System.out.println("Drawing connections");
		connectionGrid.setScale(0.0, 1.0);
		double [] blank = new double[ofc.N];
		for(int i = 0; i < ofc.N; i++){
			blank[i] = 1.0;
		}
		connectionGrid.registerData(ofc.L, ofc.L, blank);
		for(int i = 0; i < ofc.N; i++){
			double [] position = findCenter(i);
			//connectionGrid.addDrawable(Geom2D.circle(position[0], position[1], 0.005, Color.black));
			//connectionGrid.addDrawable(Geom2D.circle(position[0], position[1], 0.002, Color.black));
			for(int j = 0; j < ofc.noNbors; j++){
				int nbor = ofc.nbor[i][j];
				if(nbor>=0){
					double [] nposition = findCenter(ofc.nbor[i][j]);
					connectionGrid.addDrawable(Geom2D.line(position[0], position[1], nposition[0], nposition[1], Color.black));
				}
			}
		}
		System.out.println("Finished Drawing connections");
		
	}

	double [] findCenter(int index){
		double L = (double)ofc.L;
		int xi = index%ofc.L;
		int yi = index/ofc.L;
		double x = ((double)xi + 0.5)/L;
		double y = ((double)yi + 0.5)/L;
		double [] ret = new double [2];
		ret[0] = x;
		ret[1] = y;
		return ret;
	}

}
