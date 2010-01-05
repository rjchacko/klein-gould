package chris.ofcdamage.apps;

import java.awt.Color;
import java.io.File;
import java.text.DecimalFormat;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.ofcdamage.damage;
import chris.old.ofc.Damage2D;
import chris.util.ReadInUtil;

public class ergodicDiagramApp extends Simulation{

	damage model;
	
	Grid gridS = new Grid ("Stress");
	Grid gridL = new Grid ("Lives");
	
	ColorPalette palette1;
	ColorGradient cGradient;
	
	DecimalFormat fmt = new DecimalFormat("0000.00");
	
	private double ScMax;
	private Boolean pretime;
	private int pt, ptmax;
	private boolean draw;
//	private FitUtil fitter;
	private ReadInUtil imprt;
	
	public static void main(String[] args) {
		new Control(new ergodicDiagramApp(), "OFC Parameters");
	}
	
public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Random Seed",0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)",(int)(30));
		params.add("Lattice Size",1<<8);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Bordered"));
		params.add("Intitial Stess", new ChoiceValue("Random","Constant"));
		params.add("Min Lives", 2);	
		params.add("Max Lives", 10);
		params.add("Failure Stress (\u03C3_f)",2.0);
		params.add("\u03C3_f width",(double)(0));
		params.add("Residual Stress (\u03C3_r)",1.25);
		params.add("\u03C3_r width",0.5);
		params.add("Dissipation (\u03B1)",new DoubleValue(0.01,0,1));
		params.add("\u03B1 width", (double)(0));
		params.add("Equil Time", 100000);
		params.add("Trend Time", 100000);
		params.add("Animation", new ChoiceValue("Off","On"));
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("Number of Plate Updates");
		params.add("Last Avalanche Size");
		params.add("N_dead");
		
		c.frameTogether("OFC Model with Damage", gridS, gridL);
		
	}
	
	public void run() {
		
		// begin looping over the range
		// r = 1, 5, 10, 15, 20, 25, ...
		// r = 5*n ; n > 0
		
		for (int Rstep = 0 ; Rstep < 1 ; Rstep++){
			
			if(Rstep ==0){
				params.set("Interaction Radius (R)",1);
			}
			else{
				params.set("Interaction Radius (R)",10);
			}
			
			// begin looping over noise size
			// eta = 0.5 / 2^n = 1 / 2^(n+1)
			
			int etaStep = 0;
			while (true){
				params.set("\u03C3_r width",1/(Math.pow(2,etaStep)));
				// use the fact that for R++ noise -- (or so we would think)
				
				// run simulation
				// Setup model
				params.set("N_dead","Initializing");
				model = new damage(params);
				
				// Setup display
				setupDisplays();
				
				// Setup output files
				configOF(Rstep, etaStep);
				
				// Run the simulation
				
				// Equilibrate
				equil();
				// Simulate w/o Damage for Data
				ideal();
				
				// read in t and 1/Omega
				// and fit M regions
				if(isErgodic(Rstep, etaStep)) break;
				etaStep++;
				
			}
			
		}
		


		
		// compare slope from fits (print this to file)
		
		// decide whether or not to break *eta* loop
		
		// next r
		
		
		
		
		return;
	}
	
	
	@SuppressWarnings("unused")
	private boolean isErgodic(int a, int b){

		double[] slopes = new double[10];
		double[] chi2   = new double[10];		
		double tofit[][], temp[], tempx[], tempy[], totx[], toty[], mbar, dmbar; 
		int Dlength;
		
		String fin = model.getOutfile();
		String fout = model.getOutdir() + File.separator + "Slopes_"+a+"_"+b+".txt";
		
		imprt   = new ReadInUtil(fin);
		tofit   = imprt.getData(new int[]{1,8},1);
		Dlength = tofit[0].length;
		totx    = new double[Dlength];
		toty    = new double[Dlength];		
		Dlength = Dlength / 10;
//		fitter  = new FitUtil(Dlength);
//		mbar    = 0;
//		
//		for (int jj = 0 ; jj < 10 ; jj++){
//			
//			tempx = new double [Dlength];
//			tempy = new double [Dlength];
//			for (int kk = 0 ; kk < Dlength ; kk++){
//				tempx[kk] = tofit[0][kk+jj*Dlength];
//				tempy[kk] = 1/tofit[1][kk+jj*Dlength];
//				totx[kk+jj*Dlength] = tempx[kk];
//				toty[kk+jj*Dlength] = tempy[kk];
//			}
//			
//			temp = fitter.fit(tempx, tempy, 1,false);
//			slopes[jj] = temp[0];
//			chi2[jj]   = temp[4];
//			mbar += slopes[jj];
//			PrintUtil.printlnToFile(fout, jj, slopes[jj], chi2[jj]/(Dlength-2));
//		}
//		
//		fitter = new FitUtil(totx.length);
//		temp = fitter.fit(totx, toty, 1,false);
//		System.out.println(temp[0]);
//		PrintUtil.printlnToFile(fout, -7, temp[0], temp[4]/(totx.length-2));

		

		return (b>9) ? true : false;
	}
	
	
	@SuppressWarnings("unused")
	private void configOF(){
		
		Damage2D.PrintParams(model.getOutdir()+File.separator+"Params.txt",params);	
		model.writeDataHeaders();
		
		return;
	}
	
	private void configOF(int a, int b){
		
		String apend = a+"_"+b;
		
		model.setDataFile(apend);
		Damage2D.PrintParams(model.getOutdir()+File.separator+"Params_"+apend+".txt",params);	
		model.writeDataHeaders();
		
		return;
	}
	
	
	private void equil(){
	
		// Equilibrate
		pretime = true;
		pt      = 0;
		ptmax   = params.iget("Equil Time");
		
		model.setEquil(true);
		model.runClocks(false);
		while(pt < ptmax){
			model.equilibrate();
			pt++;
			Job.animate();
		}
		
		return;
	}
	
	private void ideal(){
		
		// Simulate w/o Damage for Data
		params.set("N_dead","Est. Equilib");
		model.setEquil(true);
		model.runClocks(true);		
		int pt0 = 0;
		ptmax = params.iget("Trend Time"); 
		pt    = params.iget("Trend Time");
		while(pt0 < ptmax){
			model.avalanche();
			model.takeData();
			pt++;
			pt0++;
			Job.animate();
		}
		
		return;
	}
	
	@SuppressWarnings("unused")
	private void damage(){
		
		Job.animate();
		pretime = false;
		model.setEquil(false);
		model.runClocks(true);
		while(model.avalanche()){
			model.takeData();
			Job.animate();
			if (params.sget("Record").equals("On")){
				// FIX ME!!!!!!!!!!!!
				model.takePicture(gridS, true);
				model.takePicture(gridL, true);
			}
		}
		
		return;
	}
	
	private void setupDisplays(){
		
		if(draw = (params.sget("Animation").equals("On"))){
		
			// Setup color scheme
			palette1  = new ColorPalette();
			cGradient = new ColorGradient();

			int MostLives = model.maxLives();
			Color[] Carray = new Color[]{Color.YELLOW,Color.RED,Color.GREEN,Color.BLUE,Color.GRAY};		
			palette1.setColor(0,Color.BLACK);
			for (int jj = 1 ; jj <= MostLives ; jj++){
				palette1.setColor(jj,Carray[jj%5]);
			}

			gridL.setColors(palette1);
			
		}
		
		return;
	}

	public void animate() {

		if(pretime){
			params.set("Number of Plate Updates", pt - ptmax);
			params.set("Last Avalanche Size",model.getAvlnchSize());
			params.set("N_dead", "Equilibrating");
		}
		else{
			params.set("Number of Plate Updates",model.getTime(1));
			params.set("Last Avalanche Size", model.getAvlnchSize());
			params.set("N_dead", model.getNdead());
		}
		
		if(draw){
			
			
			int L = model.getL();
			int N = L*L;
			
			int[] foo = model.getLives();

			double[] copyStress = model.getStress();
			
			for (int jj=0 ; jj < N ; jj++){
				cGradient.getColor(copyStress[jj],-2,ScMax);
			}

			gridS.setColors(cGradient);
			gridS.registerData(L,L,copyStress);
			gridL.registerData(L, L, foo);
			
		}
	
	}

	public void clear() {
		gridS.clear();
		gridL.clear();
	}


}
