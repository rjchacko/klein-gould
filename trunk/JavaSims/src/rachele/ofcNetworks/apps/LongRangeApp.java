/**
 * LongRangeApp is a long stress transfer range OFC app
 * 
 * @author Rachele Dominguez
 * @version 1.0 July 13, 2012
 */
package rachele.ofcNetworks.apps;

import java.io.File;

import rachele.ofcNetworks.LongRangeLattice;
import rachele.util.FileUtil;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;

public class LongRangeApp extends Simulation{


	Grid latticeGrid = new Grid("Lattice");
	LongRangeLattice ofc;
	
	Histogram sizeHist = new Histogram(1);
	Histogram cumSizeHist = new Histogram(1);
	String sizeHistFile;
	String cumSizeHistFile;
	
	public static void main(String[] args) {
		new Control(new LongRangeApp(), "Long Range OFC");
	}

	public void load(Control c) {
		c.frame(latticeGrid);
		params.add("Data Dir",new DirectoryValue("/Users/racheledominguez/data/RM/OFC_networks/testRuns/"));
		params.addm("Random Seed", 40);
		int Lstart = 64;
		params.addm("L",Lstart);
		params.addm("alpha", 0.1);
		params.addm("R", 8);
		int equilStart = Lstart*Lstart*100;
		params.addm("Equilibration Updates", equilStart);
		params.add("Av Size");
		params.add("Plate Updates");
	}

	public void clear() {
		
	}
	
	public void animate() {
		params.set("Plate Updates", ofc.plateUpdates);
		params.set("Av Size", ofc.failCount);
		latticeGrid.setScale(0.0, 1.0);
		latticeGrid.registerData(ofc.L, ofc.L, ofc.stress);

	}
	
	public void run() {
		//construct a new lattice
		ofc = new LongRangeLattice(params);
		Job.animate();
		initFiles();
		int shCount = 1;
		ofc.plateUpdates = -params.iget("Equilibration Updates");

		while(true){
			//equilibrate
			if(ofc.plateUpdates <= 0){
				ofc.ofcStep();
				Job.animate();
			}else{
				ofc.ofcStep();
				sizeHist.accum(ofc.failCount);
				int count = ofc.failCount;
				while(count>0){
					cumSizeHist.accum(count);
					count--;
				}

				if(ofc.plateUpdates % 10000 == 0) writeAccumulatedData(); 
	//			if(ofc.plateUpdates%10000000==0){
				if(ofc.plateUpdates%100000==0){
					writeShFile(shCount);
					shCount += 1;
				}
				Job.animate();
				if(ofc.plateUpdates >= 10000000){
					Job.signalStop();
				}
			}
		}
	}

	void writeShFile(int i){
		String parentDir = params.sget("Data Dir") + File.separator;
		String newShFile = parentDir + "sh" + i + ".txt";
		String newCumShFile = parentDir + "sh" + i + ".txt";
		FileUtil.initFile(newShFile, params, " Avalanch Size Histogram Data File");
		FileUtil.printHistToFile(newShFile, sizeHist);
		FileUtil.initFile(newCumShFile, params, " Cumulative Avalanch Size Histogram Data File");
		FileUtil.printHistToFile(newCumShFile, cumSizeHist);
	}
	
	void writeAccumulatedData(){
		FileUtil.initFile(sizeHistFile, params, " Avalanch Size Histogram Data File");
		FileUtil.printHistToFile(sizeHistFile, sizeHist);
		FileUtil.initFile(cumSizeHistFile, params, " Cumulative Avalanch Size Histogram Data File");
		FileUtil.printHistToFile(cumSizeHistFile, cumSizeHist);
		}
	
	void initFiles(){
		String parentDir = params.sget("Data Dir") + File.separator;
		sizeHistFile = parentDir + "sh.txt";
		cumSizeHistFile = parentDir + "ch.txt";
	}



}
