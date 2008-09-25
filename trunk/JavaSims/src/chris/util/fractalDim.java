package chris.util;

import java.awt.Color;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;

import javax.imageio.ImageIO;

import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.FileValue;

public class fractalDim extends Simulation{

	private Grid gridC = new Grid("Cluster");
	private Grid gridV = new Grid("Visualization");
	private ColorPalette palette;
	
	private LatticeNeighbors LN;
	
	private int pc, cn[], dX[], dY[], counter, countMAX, L, N, cOFm[], disp[];
	private String basename, baseINname, datafile, cumfile;
	private boolean draw;
	
	
	public static void main(String[] args) {
		new Control(new fractalDim(), "Measure Fractal Dimension");
	}
	
	public void load(Control c) {
		
		params.add("Example File", new FileValue("/Users/cserino/Desktop/"));
		params.add("Number of Files", (int) 1);
		params.add("L", (int) 256);
		params.add("Animation", new ChoiceValue("Off","On"));
		params.add("Files Left");
		
		c.frameTogether("Cluster Analysis", gridC, gridV);
		
	}

	public void run() {
		
		L = params.iget("L");
		N = L*L;
		
		cn = new int[N];
		dX = new int[N];
		dY = new int[N];
		
		// setup displays
		setDisplay();
		
		// find and set the file basename
		setBaseInname();
		setBasename(params.sget("Example File"));
		cumfile = basename + "_AllFractalDim.txt";
		dataHeader(cumfile);
		
		countMAX = params.iget("Number of Files");
		counter  = countMAX;
		
		while(counter > 0){

			params.set("Files Left", counter);

			// define outfile
			datafile = basename + "_" + counter + "_fractalDim.txt";
			dataHeader();
			
			// fill the cluster arrays
			setArrays();
			Job.animate();
			
			// calculate center of mass of cluster
			cOFm = getCM();
			setcOFmDisplay();
			Job.animate();
					
			// calculate the density as a function of radius
			getDensity();
			
			// take picture of percolating cluster with center of mass
			capturePC();
			
			counter--;
		}
		
		Job.animate();
		params.set("Files Left", "Done");

		return;
	}
	
	private void capturePC(){
		
		if(!draw) return;

		String fn = basename + "_" + counter + "_fractal.png";
		
		try {
			ImageIO.write(gridV.getImage(), "png", new File(fn));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + fn);
		}
		
		
		
		return;
	}
	
	private void setBasename(String fin){

		int findslash = fin.indexOf('/');
		String ret = fin.substring(0,findslash+1);
		fin = fin.substring(findslash+1);
		findslash = fin.indexOf('/');
		while(findslash != -1){
			ret+=fin.substring(0,findslash+1);
			fin = fin.substring(findslash+1);
			findslash = fin.indexOf('/');
		}
		ret += "FractalDim/";
		DirUtil.MkDir(ret);
		basename = ret + File.separator +"Clusters";
		
		return;
	}
	
	private void getDensity(){
		
		int[] nbs;
		double mass;
		
		
		for(int jj = 1 ; jj <= (L-1)/2 ; jj++){	// L' x L' square has L' = 2*R + 1 <= L
			mass = 0;
			LN = new LatticeNeighbors(L, L,0,jj,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Square);
			nbs = LN.get(cOFm[0] + L*cOFm[1]);
			for(int kk = 0 ; kk < nbs.length ; kk++){
				mass += cn[nbs[kk]];
			}
			double tempr = 2*jj + 1;
			recordDensity(tempr, mass);
			// for display purposes
			drawBox(jj);
			Job.animate();
		}
		
		

		return;
	}
	
	private void drawBox(int R){
		
		if(!(R > 8 && draw)) return;
		
		int[] box;
		
		LN = new LatticeNeighbors(L, L,R-1,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Square);

		box = LN.get(cOFm[0] + L*cOFm[1]);
		for(int kk = 0 ; kk < box.length ; kk++){
			disp[box[kk]] = 3;
		}
		
		LN = new LatticeNeighbors(L, L,R-2,R-1,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Square);

		box = LN.get(cOFm[0] + L*cOFm[1]);
		for(int kk = 0 ; kk < box.length ; kk++){
			disp[box[kk]] = cn[box[kk]];
		}
		
		return;
	}
	
	private void dataHeader(){
		
		PrintUtil.printlnToFile(datafile, "Length", "Mass", "Density", "Log L", "Log Density");
		return;
	}
	
	private void dataHeader(String fout){
		
		PrintUtil.printlnToFile(fout, "Length", "Mass", "Density", "Log L", "Log Density");
		return;
	}
	
	private void recordDensity(double rad, double mass){
		
		double rho = mass/(rad*rad);
		
		PrintUtil.printlnToFile(datafile,rad, mass, rho, Math.log(rad), Math.log(rho));
		PrintUtil.printlnToFile(cumfile,rad, mass, rho, Math.log(rad), Math.log(rho));
		return;
	}
	
	private void setcOFmDisplay(){
		
		if(!draw) return;
		
		LN = new LatticeNeighbors(L, L,0,7,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Circle);

		int[] temp = LN.get(cOFm[0] + L*cOFm[1]);
	
		for (int jj = 0 ; jj < temp.length; jj++){
			disp[temp[jj]] = 2;
		}
		
		return;
	}
	
	private int[] getCM(){
		
		int cX = 0;
		int cY = 0;
		int ms = 0;
		int[] cv = new int[2];	// turn cluster coordinates into grid coordinates
		
		for(int jj = 0 ; jj < N ; jj++){
			ms += cn[jj];
			cX += cn[jj]*dX[jj];
			cY += cn[jj]*dY[jj];
			if(dX[jj] == 0 && dY[jj] == 0){
				cv[0] = jj%L;
				cv[1] = (int)((double)jj / (double)L);
			}
			
		}
		
		cX = cX/ms + cv[0];
		cY = cY/ms + cv[1];

		if (cX < 0) cX = L + cX;
		if (cY < 0) cY = L + cY;
		if (cX > L) cX = cX - L;
		if (cY > L) cY = cY - L;
		
		return new int[]{cX, cY};
	}
	
	private void setDisplay(){
		
		if(draw = (params.sget("Animation").equals("On"))){
		
			palette = new ColorPalette();
			
//			Color[] Carray = new Color[]{Color.RED, Color.ORANGE, Color.YELLOW, Color.GREEN, Color.BLUE,
//					 Color.GRAY, Color.PINK, Color.MAGENTA, Color.CYAN, Color.PINK}; 
			Color[] Carray = new Color[]{Color.BLUE, Color.YELLOW, Color.RED, Color.GREEN};    
			palette.setColor(0,Color.WHITE);
			palette.setColor(1,Color.BLACK);
			for (int jj = 2 ; jj <= N ; jj++){
				palette.setColor(jj,Carray[jj%4]);
			}

			gridC.setColors(palette);
			gridV.setColors(palette);
		
		}
		
		return;
	}
	
	private void setArrays(){
		
		String fin1 = baseINname + "_" + counter + ".txt";
		String fin2 = baseINname + "DX_" + counter + ".txt";
		String fin3 = baseINname + "DY_" + counter + ".txt";
		
		try {
			
			int itemp1, itemp2, itemp3;
			
			FileInputStream fis1 = new FileInputStream(fin1);
			BufferedInputStream bis1 = new BufferedInputStream(fis1);
			BufferedReader bir1 = new BufferedReader(new InputStreamReader(bis1));
			
			FileInputStream fis2 = new FileInputStream(fin2);
			BufferedInputStream bis2 = new BufferedInputStream(fis2);
			BufferedReader bir2 = new BufferedReader(new InputStreamReader(bis2));
			
			FileInputStream fis3 = new FileInputStream(fin3);
			BufferedInputStream bis3 = new BufferedInputStream(fis3);
			BufferedReader bir3 = new BufferedReader(new InputStreamReader(bis3));
			
			String stemp = bir1.readLine();
			itemp1       = stemp.indexOf('\t');
			pc           =  (int) Double.parseDouble(stemp.substring(itemp1+1));
			
			for (int jj = 0 ; jj < N ; jj++){
				itemp1 = (int) Double.parseDouble(bir1.readLine());
				itemp2 = (int) Double.parseDouble(bir2.readLine());
				itemp3 = (int) Double.parseDouble(bir3.readLine());
				if(itemp1 == pc){
					cn[jj] = 1;
					dX[jj] = itemp2;
					dY[jj] = itemp3;
				}
				else{
					cn[jj] = 0;
					dX[jj] = -2*N;
					dY[jj] = -2*N;
				}
//	
//				cn[N - jj - 1] = itemp1;
//				dX[N - jj - 1] = itemp2;
//				dY[N - jj - 1] = itemp3;
			
			}			
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
			params.set("Files Left","Error!");
		} catch (IOException e) {
			e.printStackTrace();
			params.set("Files Left","Error!");
		}
		
		disp = CopyArray.copyArray(cn, N);
		
		return;
	}
	
	private void setBaseInname(){
		
		String fin = params.sget("Example File");
		int findUS = fin.indexOf('_');
		baseINname   = fin.substring(0,findUS);
		return;
	}
	

	public void animate() {
		
		if(draw){
			gridC.registerData(L, L, cn);
			gridV.registerData(L,L,disp);
		}
		
		
		return;
	}

	public void clear() {
			
	}

	
}
