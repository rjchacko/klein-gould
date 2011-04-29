package rachele.damage2D.apps;

import java.io.File;

import rachele.damage2D.multidx.FrozenDamageLattice;
import rachele.util.FileUtil;
import rachele.util.MathTools;
import rachele.util.ReadWriteTextFile;
import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.util.DoubleArray;

public class FrozenDamageIterConnectParam extends Simulation{
	
	int dt;
	int pow;
	int maxSize;	
	int cg_dt;
	
	String sizeHistFile;
	String sizeFile;
	String cpFile;
	String cpinfoFile;
	String inputFile;
	String boundaryConditions;
	double alphaDissRate;
	double [] cParam;
	int [][] nborList;
		
	FrozenDamageLattice ofc;
		
	Histogram sizeHist = new Histogram(1);
	
	Accumulator maxSizeTempAcc = new Accumulator();	
	Accumulator scaledTimeAcc = new Accumulator();
	Accumulator scaledLiveTimeAcc = new Accumulator();
	Accumulator alphaAcc = new Accumulator();

	Grid deadGrid = new Grid("Lattice");
	Grid gammaGrid = new Grid("Alpha Prime");
	Plot alphaPlot = new Plot("Alpha Histogram");
	Plot alphaPlote = new Plot("Alpha Histogram");
    Grid iter1 = new Grid("Iteration 1");
	Grid iter2 = new Grid("Iteration 2");
	Grid iter3 = new Grid("Iteration 3");
	Grid iter4 = new Grid("Iteration 4");
	Grid iter5 = new Grid("Iteration 5");
	Grid iter6 = new Grid("Iteration 6");
	Grid iter7 = new Grid("Iteration 7");
	
	
	public static void main(String[] args) {
		new Control(new FrozenDamageIterConnectParam(), "OFC Connection Param Calculator");
	}

	public void load(Control c) {

		deadGrid = new Grid("Dead Sites");
		gammaGrid = new Grid("Gamma");
		c.frameTogether("All Data",deadGrid, gammaGrid, iter1, iter2, iter3, iter4, iter5, iter6, iter7);
		params.add("Iterations", 36);
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
		params.add("Type of Damage", new ChoiceValue( rd, ac, pd, bl, cs, br, ds, pr, cr, db, dr));
//		params.add("Type of Damage", "CascadeRandom");//, rd, ac, pd, bl, cs, br, ds, pr, cr, db, dr));
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
		params.add("Alpha Distribution", new ChoiceValue(da, mg, ca, ei, qu , fr, mg, da, ah, gs, ga, fa, az));
//		params.add("Alpha Distribution", "Constant");//, mg, ca, ei, qu , fr, mg, da, ah, gs, ga, fa, az));
		params.addm("Random Seed", 0);
		params.addm("Size Power",0);
		params.addm("R", 0);
		params.addm("Init Percent Dead", 0.0);
		params.addm("Dead Parameter", 0);
		params.addm("Number Dead", 0);
		params.addm("Coarse Grained dt (PU)", 0);
		params.addm("Equilibration Updates", 0);
		params.addm("Max PU",0);
		params.addm("Data points per write", 0);
		params.addm("Residual Stress", 0.0);
		params.addm("Res. Max Noise", 0.0);
		params.addm("Dissipation Param", 0.0);
		params.addm("Damage Tolerance", 0.0);
		params.add("Av Size");
		params.add("Plate Updates");
		params.add("Percent Damage");
		params.add("Effective Alpha");
		params.add("Error");
	}

	public void animate() {
		params.set("Plate Updates", ofc.plateUpdates);
		params.set("Av Size", ofc.avSize);
	}

	public void clear() {
	}

	
	public void run() {
		cpFile = params.sget("Data Dir") + File.separator + "cpi.txt";
		cpinfoFile = params.sget("Data Dir") + File.separator + "cpi_info.txt";
		inputFile = params.sget("Data Dir") + File.separator + "info.txt";
		setParams(inputFile);
		FileUtil.initFile(cpFile, params);
		System.out.println("Init OFC");
		ofc = new FrozenDamageLattice(params, cpinfoFile);
		System.out.println("Init OFC LAttice");
		ofc.initLattice(params);
		pow = params.iget("Size Power");
		
		cParam = new double [ofc.N];
		
		int dx = 1;
		int Lp = ofc.L/dx;
		int Np = Lp*Lp;
		double [] deadSites = new double [Np];
		for (int i = 0; i < Np; i++){
			if(ofc.aliveLattice[i]) deadSites[i]=1.0;
			else deadSites[i]=0;
		}
		deadGrid.setScale(0.0, 1.0);
		deadGrid.registerData(Lp, Lp, deadSites);
		for (int i = 0; i < Np; i++){
			if(ofc.aliveLattice[i]) deadSites[i]=1.0 - ofc.gamma[i];
			else deadSites[i]=0;
		}
		gammaGrid.setScale(0.0, 1.0);
		gammaGrid.registerData(Lp, Lp, deadSites);

		//init c param
		for (int i = 0; i < ofc.N; i ++){
			cParam[i]=1.0-ofc.gamma[i];
		}
		for (int i=0;i<params.iget("Iterations");i++){
			drawLattices(i);
		}


	}
	
	void setParams(String file){
		File input = new File(file);
		params.set("Type of Damage", ReadWriteTextFile.getLastWordOfLine(input, 1));
		params.set("Dead dissipation?", ReadWriteTextFile.getLastWordOfLine(input, 2));
		params.set("Boundary Conditions", ReadWriteTextFile.getLastWordOfLine(input, 3));
		params.set("Alpha Distribution", ReadWriteTextFile.getLastWordOfLine(input, 4));//, mg, ca, ei, qu , fr, mg, da, ah, gs, ga, fa, az));
		params.set("Random Seed", Integer.parseInt(ReadWriteTextFile.getLastWordOfLine(input, 5).trim()));
		params.set("Size Power",Integer.parseInt(ReadWriteTextFile.getLastWordOfLine(input, 6).trim()));
		params.set("R", Integer.parseInt(ReadWriteTextFile.getLastWordOfLine(input, 7).trim()));
		params.set("Init Percent Dead", Double.parseDouble(ReadWriteTextFile.getLastWordOfLine(input, 8)));
		params.set("Dead Parameter", Integer.parseInt(ReadWriteTextFile.getLastWordOfLine(input, 9).trim()));
		params.set("Number Dead", Integer.parseInt(ReadWriteTextFile.getLastWordOfLine(input, 10).trim()));
		params.set("Coarse Grained dt (PU)", Integer.parseInt(ReadWriteTextFile.getLastWordOfLine(input, 11).trim()));
		params.set("Equilibration Updates", Integer.parseInt(ReadWriteTextFile.getLastWordOfLine(input, 12).trim()));
		params.set("Max PU",Integer.parseInt(ReadWriteTextFile.getLastWordOfLine(input, 13).trim()));
		params.set("Data points per write", Integer.parseInt(ReadWriteTextFile.getLastWordOfLine(input, 14).trim()));
		params.set("Residual Stress", Double.parseDouble(ReadWriteTextFile.getLastWordOfLine(input, 15)));
		params.set("Res. Max Noise", Double.parseDouble(ReadWriteTextFile.getLastWordOfLine(input, 16)));  
		params.set("Dissipation Param", Double.parseDouble(ReadWriteTextFile.getLastWordOfLine(input, 17)));
		params.set("Damage Tolerance", Double.parseDouble(ReadWriteTextFile.getLastWordOfLine(input, 18)));
	
		boundaryConditions = params.sget("Boundary Conditions");
	}


	void drawLattices(int i){

		int dx = 1;
		int Lp = ofc.L/dx;
//		int Np = Lp*Lp;
//		FileUtil.printlnToFile(cpFile, "# FileUtil.printlnToFile(cpFile, area, cave, liveAve, csum, i, cvar);");

		System.out.println("iteration = " + i);
		calcCP();
		if(i==0) iter1.registerData(Lp, Lp, cParam);
		if(i==1) iter2.registerData(Lp, Lp, cParam);
		if(i==2) iter3.registerData(Lp, Lp, cParam);
		if(i==3) iter4.registerData(Lp, Lp, cParam);
		if(i==4) iter5.registerData(Lp, Lp, cParam);
		if(i==5) iter6.registerData(Lp, Lp, cParam);
		if(i==6) iter7.registerData(Lp, Lp, cParam);
//		if(i==7) iter8.registerData(Lp, Lp, cParam);
		
		double cave = MathTools.mean(cParam);
		double cvar = MathTools.variance(cParam);
		double csum = DoubleArray.sum(cParam);
		double sum = 0;
		int count = 0;
		for(int j = 0; j < ofc.N; j++){
			sum+=cParam[j];
			if(cParam[j]!=0) count+=1;
		}
		double liveAve = sum/(double)count;
//		double area = Math.PI*Math.pow(i, 2);
		FileUtil.printlnToFile(cpFile, i, cave, liveAve, csum, cvar);
		System.out.println("iter = " + i+1 + " Ave = " + cave + " Var = " + cvar + " Sum = " + csum);



	}
	
	void calcCP(){
		double [] temp = new double [ofc.N];
		for (int i = 0; i < ofc.N; i ++){
			double sum = 0.0;
			double thresholdCt=ofc.noNbors;
			for(int j=0; j<ofc.noNbors; j++){
				double s = cParam[ofc.nborList[i][j]];
				sum+=s;
				if(s==0) thresholdCt-=1.0;
				
			}
//			temp[i]=sum*cParam[i]/thresholdCt;
			temp[i]=sum*cParam[i]/ofc.noNbors;
		}
		for (int i = 0; i < ofc.N; i ++){
			cParam[i] =  temp[i];
		}
	}


}
