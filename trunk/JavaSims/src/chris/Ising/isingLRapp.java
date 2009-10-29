package chris.Ising;

import static scikit.util.Utilities.format;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

import chris.util.PrintUtil;

import kip.ising.dim2.IsingLR;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

public class isingLRapp extends Simulation{
	
	public static void main(String[] args) {
		new Control(new isingLRapp(), "Ising Model");
	}
	
	Grid grid = new Grid("Coarse Grained Display");
	Grid dirt = new Grid("Spin Level Display");
	int dx;
	double ebar[], p;
	IsingLR sim;
	private static int dlength = 100000;
	private static double eqt = 500;
	private static double rct = 5000;
	
	double data[] = new double[dlength]; 
	
	private DecimalFormat pfmt = new DecimalFormat("000");
	
	public void load(Control c) {
		c.frameTogether("Displays",grid,dirt);
		params.addm("Dynamics", new ChoiceValue("Ising Metropolis", "Ising Glauber", "Kawasaki Metropolis", "Kawasaki Glauber"));
		params.addm("Scale colors", new ChoiceValue("False", "True"));
		params.add("Random seed", 0);
		params.add("L", 1<<8);
		params.add("R", 1<<4);
		params.add("Initial magnetization", 0.);
		params.add("Dilution", new DoubleValue(0.,0.,1.));
		params.addm("T", 10.);
		params.addm("J", 1.0);
		params.addm("h", 0.0);
		params.addm("dt", 1);
		params.add("time");
		params.add("magnetization");
	}
	
	public void animate() {
		params.set("time", format(sim.time()));
		params.set("magnetization", format(sim.magnetization()));
		sim.setParameters(params);
		
		if (params.sget("Scale colors").equals("False"))
			grid.setScale(-1, 1);
		else
			grid.setAutoScale();
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
		dirt.registerData(sim.L, sim.L, sim.getSpins());
		
	}
	
	public void clear() {
		grid.clear();
		dirt.clear();
	}
	
	public void run() {
		double pnow = 0;
				
		while (pnow < 1){
	
			
			setupColors();

			sim = new IsingLR(params);

			sim.setDiluteField(params.fget("Initial magnetization"),params.fget("Dilution"));

			dx = Math.max(Integer.highestOneBit(sim.R)/8, 1);

			ebar = new double[sim.L*sim.L];

			p    = params.fget("Dilution");


			if(pnow == 0){
				PrintUtil.printlnToFile("/Users/cserino/Desktop/ising1/Params.log",params.toString());
				PrintUtil.printlnToFile("/Users/cserino/Desktop/ising1/Params.log", sim.getClass().getName());
				DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		        Date date = new Date();
		        PrintUtil.printlnToFile("/Users/cserino/Desktop/ising1/Params.log", dateFormat.format(date));
			}
			
			double lastUpdate = 0;
			while (true) {
				while (sim.time() - lastUpdate < sim.dt()) {
					sim.step();
					Job.animate();
				}
				lastUpdate = sim.time();
				if(lastUpdate > eqt) break;
				Job.animate();
			}
			
			sim.resetTime();
			lastUpdate = 0;
			while (true) {
				while (sim.time() - lastUpdate < sim.dt()) {
					sim.step();
					Job.animate();
				}
				lastUpdate = sim.time();
				if(lastUpdate > rct) break;
				data[(int)(lastUpdate-sim.dt())] = calcOmega();
				Job.animate();
			}
			
			pnow += 0.01;
			pnow = (double)(Math.round(100*pnow))/100;
			params.set("Dilution", pnow);
			printData(pnow);
			data = new double[dlength]; 
			
		}
	}
	
	private void setupColors(){
		
		dirt.setScale(-1,1);
		ColorPalette palette = new ColorPalette();
		palette.setColor(-1,Color.BLACK);
		palette.setColor(1, Color.WHITE);
		palette.setColor(0,Color.RED);
		dirt.setColors(palette);
		return;
	}
	
	public double calcOmega(){
		double ebarbar = 0;
		double omega   = 0;
		int x, y, s;
		for (int jj = 0 ; jj < sim.L*sim.L ; jj++){
			x = jj%sim.L;
			y = jj/sim.L;
			s = sim.spins.get(x,y);
			ebar[jj] += -sim.dt()*s*(sim.h+sim.J*(sim.spins.sumInRange(x,y)-s)/(4*sim.R*sim.R));
			ebarbar  += ebar[jj];
		}
		ebarbar /= ((1.-p)*sim.L*sim.L);
		for (int jj = 0 ; jj < sim.L*sim.L ; jj++)
			omega += (ebar[jj]-ebarbar)*(ebar[jj]-ebarbar);
		omega /= (sim.time()*sim.time());
		return omega;
	}
	
	public void printData(double pp){
		try{
			File file = new File("/Users/cserino/Desktop/ising1/metric_"+pfmt.format(100*(1.-pp))+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for(int jj = 0 ; jj < Math.min(rct, dlength) ; jj++){
				pw.print(jj);
				pw.print("\t");
				pw.println(data[jj]);
			}
		} catch (IOException ex){
			ex.printStackTrace();
		}
	}
	
	
}
