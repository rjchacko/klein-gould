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
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.util.DoubleArray;

public class FrozenDamageConnectionParamApp  extends Simulation{
	
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
    Grid cp2 = new Grid("cp2");
	Grid cp4 = new Grid("cp4");
	Grid cp6 = new Grid("cp6");
	Grid cp8 = new Grid("cp8");
	Grid cp10 = new Grid("cp10");
	Grid cp12 = new Grid("cp12");
	Grid cp16 = new Grid("cp32");

	public static void main(String[] args) {
		new Control(new FrozenDamageConnectionParamApp(), "OFC Damage Model Scaling Only");
	}

	public void load(Control c) {

		deadGrid = new Grid("Dead Sites");
		gammaGrid = new Grid("Gamma");
		c.frameTogether("All Data",deadGrid, gammaGrid, cp2, cp4, cp6, cp8, cp10, cp12,cp16);
		params.add("maxRange", 16);
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/contract2/testRuns"));
//		String cr = "Cascade Random";
//		params.add("Type of Damage", new ChoiceValue("Random"));//, rd, ac, pd, bl, cs, br, ds, pr, cr, db, dr));
		params.add("Type of Damage", "Random");//, rd, ac, pd, bl, cs, br, ds, pr, cr, db, dr));
		params.add("Dead dissipation?", new ChoiceValue("Yes", "No") );
		params.add("Boundary Conditions", new ChoiceValue("Periodic", "Open"));
//		String ca = "Constant";
		params.add("Alpha Distribution", "ManyGaussians");//, mg, ca, ei, qu , fr, mg, da, ah, gs, ga, fa, az));
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
		cpFile = params.sget("Data Dir") + File.separator + "cpn.txt";
		cpinfoFile = params.sget("Data Dir") + File.separator + "cpninfo.txt";
		inputFile = params.sget("Data Dir") + File.separator + "info.txt";
		setParams(inputFile);
		FileUtil.initFile(cpFile, params);
		ofc = new FrozenDamageLattice(params, cpinfoFile);
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
//		double percentAlive = DoubleArray.sum(deadSites)/(double)Np;
		deadGrid.setScale(0.0, 1.0);
		deadGrid.registerData(Lp, Lp, deadSites);
		for (int i = 0; i < Np; i++){
			if(ofc.aliveLattice[i]) deadSites[i]=1.0 - ofc.gamma[i];
			else deadSites[i]=0;
		}
		gammaGrid.setScale(0.0, 1.0);
		gammaGrid.registerData(Lp, Lp, deadSites);
		int maxDist = params.iget("maxRange");
//		int maxNoNbors = findNoCircleNbors(maxDist);
		//		nborList = new int [ofc.N][maxNoNbors];
//		int i=1;
//		while(true){
			for(int i = 1; i <= maxDist; i++){
				drawLattices(i);

//			if(i>=maxDist) Job.signalStop();
//			i++;
		}

	}
	
	void setParams(String file){
		File input = new File(file);
		String dam = ReadWriteTextFile.getLastWordOfLine(input, 1);
		String test = ReadWriteTextFile.getLastWordOfLine(input, 4);
		System.out.println(test + "test");
		params.set("Type of Damage", dam);
		params.set("Dead dissipation?", ReadWriteTextFile.getLastWordOfLine(input, 2));
		params.set("Boundary Conditions", ReadWriteTextFile.getLastWordOfLine(input, 3));

		params.set("Alpha Distribution", ReadWriteTextFile.getLastWordOfLine(input, 4));//, mg, ca, ei, qu , fr, mg, da, ah, gs, ga, fa, az));
		params.set("Alpha Distribution", ReadWriteTextFile.getLastWordsOfLine(input, 4,3));//, mg, ca, ei, qu , fr, mg, da, ah, gs, ga, fa, az));

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
		int Np = Lp*Lp;
		FileUtil.printlnToFile(cpFile, "# FileUtil.printlnToFile(cpFile, area, cave, liveAve, csum, i, cvar);");

		System.out.println("maxRange = " + i);
		calcCP(i);
		if(i==1) cp2.registerData(Lp, Lp, cParam);
		if(i==4) cp4.registerData(Lp, Lp, cParam);
		if(i==6) cp6.registerData(Lp, Lp, cParam);
		if(i==8) cp8.registerData(Lp, Lp, cParam);
		if(i==10) cp10.registerData(Lp, Lp, cParam);
		if(i==12) cp12.registerData(Lp, Lp, cParam);
		if(i==32) cp16.registerData(Lp, Lp, cParam);

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
		double area = Math.PI*Math.pow(i, 2);
		FileUtil.printlnToFile(cpFile, area, cave, liveAve, csum, i, cvar);
		System.out.println("maxRange = " + i + " Ave = " + cave + " Var = " + cvar + " Sum = " + csum);



	}
	
	void calcCP(int maxRange){
//		int noRangeNbors = findNoCircleNbors(maxRange);
//		findCircleNbors(maxRange,noRangeNbors);
		for (int i = 0; i < 1; i ++){
			double sum = 1.0;
//			if(i==1) cp2.registerData(ofc.L, ofc.L, cParam);
//			Job.animate();
//			System.out.println("site = " + i);
//			for(int n=0; n<noRangeNbors; n++){
//				sum*=(1.0-ofc.gamma[nborList[i][n]]);//*Math.pow(1.0/(double)R+1.0-nborDistance[n], 2);///(1.0-ave);
//			}

//			for(int i=0; n<ofc.N; n++){
				int x = i%ofc.L;
				int y = i/ofc.L;
				for(int dy = -maxRange; dy <= maxRange; dy++){
					for(int dx = -maxRange; dx <= maxRange; dx++){
						double distance = Math.sqrt(dx*dx + dy*dy);
						if (distance <= maxRange){
							if(distance !=0){

								int xx = (x+dx+ofc.L)%ofc.L;
								int yy = (y+dy+ofc.L)%ofc.L;
								int nborSite = yy*ofc.L+xx;

								sum*=(1.0-ofc.gamma[nborSite]);
//								System.out.println("sum  = " + sum + " gamma " + ofc.gamma[nborSite]);
//								nborct+=1;
//								System.out.println(" no nobors " + nborct);
							}
					}
				}
			}
			cParam[i]=(1.0-ofc.gamma[i])*Math.pow(sum,(double)(1.0/(Math.PI*maxRange*maxRange)));
//			if(i%1000==0)System.out.println("c param " + i + " = " + cParam[i]);
		}
	}
//
//	
//	int findNoCircleNbors(int range){
//		int count = 0;
//		 for(int dy = -range; dy <= range; dy++){
//			 for(int dxx = -range; dxx <= range; dxx++){
//				 double distance = Math.sqrt(dxx*dxx + dy*dy);
//				 if (distance <= range){
//						 count += 1;
//				 }
//			 }
//		 }
//		 count = count -1;
//		 return count;
//	}
//	
//	void findCircleNbors(int range, int maxNo){
//		
//		if(boundaryConditions == "Periodic"){
//			int nborIndex = 0;
//			for (int s = 0; s < ofc.N; s++){
//				nborIndex = 0;
//				int x = s%ofc.L;
//				int y = s/ofc.L;
//				for(int dy = -range; dy <= range; dy++){
//					for(int dx = -range; dx <= range; dx++){
//						double distance = Math.sqrt(dx*dx + dy*dy);
//						if (distance <= range){
////							if(distance == 1.0) System.out.println("nearest neighbor = " + nborIndex + " distance = " + distance);
//							if(distance !=0){
//								int xx = (x+dx+ofc.L)%ofc.L;
//								int yy = (y+dy+ofc.L)%ofc.L;
//								int nborSite = yy*ofc.L+xx;
//								nborList[s][nborIndex] = nborSite;
//								nborIndex += 1;
//							}
//						}
//					}
//				}
//			}
////			if (nborIndex != (maxNbors)) System.out.println("Problem in findCircleNbors");
//		}else if(boundaryConditions == "Open"){
//			int nborIndex = 0;
//			for (int s = 0; s < ofc.N; s++){
//				nborIndex = 0;
//				int x = s%ofc.L;
//				int y = s/ofc.L;
//				for(int dy = -range; dy <= range; dy++){
//					for(int dx = -range; dx <= range; dx++){
//						double distance = Math.sqrt(dx*dx + dy*dy);
//						if (distance <= range){
//							int xx = (x+dx);
//							int yy = (y+dy);
//							if(xx>=0 & xx<ofc.L & yy>=0 & yy<ofc.L){
//								int nborSite = yy*ofc.L+xx;
//								nborList[s][nborIndex] = nborSite;
//							}else{
//								nborList[s][nborIndex] = -1;
//							}
//							nborIndex += 1;
//						}
//					}
//				}
//			}
//		}
//	}
//	


}
