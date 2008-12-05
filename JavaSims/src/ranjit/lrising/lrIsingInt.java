package ranjit.lrising;

import java.awt.Color;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Random;

import scikit.dataset.Accumulator;
import scikit.dataset.Histogram;
import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.FileValue;

public class lrIsingInt extends Simulation {

	int spins[]=null;
	int L;
	int N;
	double T;
	double h;
	int M;
	double E;
	double J;
	Grid grid = new Grid ("Ising Lattice");
	
	Accumulator magAccum=new Accumulator();	
	Plot magPlot = new Plot("Magnetizations");
	
	Histogram interventionHist=new Histogram(1);
	Plot interventionPlot=new Plot("Intervention");
	
	Accumulator initMagAccum=new Accumulator();
	Plot initMagPlot=new Plot("Initial Magnetization");
	
	ArrayList<int[]> history=new ArrayList<int[]>();
	
	public void animate() {
		ColorPalette palette = new ColorPalette();
		palette.setColor(0,Color.BLACK);
		palette.setColor(1,Color.WHITE);

		ColorGradient smooth = new ColorGradient();
		for (int i=0 ; i<spins.length; i++){
			smooth.getColor(spins[i],-1,1);
		}
		grid.setColors(smooth);
		grid.registerData(L,L,spins);
		magPlot.registerPoints("magnetization", magAccum, Color.red);
		interventionPlot.registerBars("intervention", interventionHist, Color.red);
		initMagPlot.registerPoints("initmag", initMagAccum, Color.red);
	}

	@Override
	public void clear() {
		grid.clear();
		magAccum.clear();
		magPlot.clear();
		initMagAccum.clear();
		initMagPlot.clear();
	}

	@Override
	public void load(Control c) {
		params.add("Input Config",new FileValue("/Users/rjchacko/Simulations/JavaSims/history"));
		params.add("T",1.0);
		params.add("h",0.0);
		params.add("L",32);
		params.add("MCS",100000);
		params.add("mcs");
		params.add("M");
		params.add("Bond probability");
		
		c.frame(grid,magPlot,interventionPlot,initMagPlot);
	}

	@SuppressWarnings("unchecked")
	@Override
	public void run() {
		L=params.iget("L");
		try {
            FileInputStream fis = new FileInputStream(params.sget("Input Config"));
            ObjectInputStream ois = new ObjectInputStream(fis);
            history = (ArrayList<int[]>)ois.readObject();
            fis.close();
		}
		catch (Throwable e) {
            System.err.println("exception thrown");
		}
		
		for(int j=0;j<history.size();j++){
			
			J=1./(L*L);
			T=params.fget("T");
			h=params.fget("h");
			int totalMCS=params.iget("MCS");
			spins=history.get(j).clone();
			M=0;
			for(int b=0;b<spins.length;b++){
				M+=spins[b];
			}
			initMagAccum.accum(j, (double)M/((double)L*L));
			magAccum.clear();
			magPlot.clear();
			Job.animate();
			for(int n=0;n<100;n++){			
				Random r=new Random();
				spins=history.get(j).clone();
				M=0;
				for(int b=0;b<spins.length;b++){
					M+=spins[b];
				}
				
				E=J*M*M+h*M;
				
				int dM=-2;
				
				double m=(double)M/(double)(L*L);
				
				for(int mcs=0;mcs<totalMCS;mcs++){
					//Perform one MCS
					for(int k=0;k<L*L;k++){
						int nextSpin=r.nextInt(L*L);
						int sign=spins[nextSpin];
						int newM=M+sign*dM;
						double newEnergy=(J*newM*newM+h*newM);
						double dE=E-newEnergy;
						boolean decreaseEnergy=dE<=0;
						boolean acceptThermal=r.nextDouble()<Math.exp(-dE/T);
						if( decreaseEnergy || acceptThermal){
							spins[nextSpin]*=-1;
							M=newM;
							E=newEnergy;
						}
						Job.animate();
					}
					params.set("mcs", mcs);
					params.set("M", M);
					m=(double)M/(double)(L*L);
					magAccum.accum(mcs, m);
					if(m<0){
						interventionHist.accum(j);
						break;
					}
					Job.animate();
				}
			}
		}
	}


	public static void main(String[] args) {
		new Control(new lrIsingInt(),"Ising 2D Umbrella Sampling");
	}

}
