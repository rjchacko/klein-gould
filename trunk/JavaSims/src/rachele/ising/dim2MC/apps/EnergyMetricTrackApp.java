package rachele.ising.dim2MC.apps;

import static scikit.util.Utilities.format;


import java.io.File;
import rachele.ising.dim2MC.EnergyMetric;
import rachele.ising.dim2MC.IsingLR;
import rachele.util.FileUtil;
import rachele.util.FourierTransformer;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;

public class EnergyMetricTrackApp  extends Simulation{
	
	IsingLR sim;
	EnergyMetric em;
    FourierTransformer ft;
	
	Grid grid = new Grid("Long Range Ising Model");
	
	int dx;
	int recordStep = 0;
    boolean clearFile;
    boolean takeAverages;
    boolean writeSFFiles = true;
    String outputFileName, outputFileName0, sfFileNameH, sfFileNameV, sfFileNameDiff, energyFileName;
    int noParticles;
    
   	Accumulator eMetricInverseAcc = new Accumulator();
   	Accumulator metric0Acc = new Accumulator();
   	
	public static void main(String[] args) {
		new Control(new EnergyMetricTrackApp(), "Ising Model w/ particle tracking");
	}
	
	public void load(Control c) {
		c.frame(grid);	
		params.add("Output Directory",new DirectoryValue("/Users/erdomi/data/ergodicityTest/testruns"));
		params.addm("Dynamics", new ChoiceValue("Kawasaki Glauber", "Kawasaki Metropolis", "Glauber"));
		params.addm("init", new ChoiceValue( "Random"));
		params.add("Random seed", 0);
		params.add("L", 128);
		params.add("R", 46);
		params.add("JR", 46);
		params.add("R offset", 0);
		params.add("Initial magnetization", 0.0);
		params.addm("T", .1);
		params.addm("J", -1.0);
		params.addm("h", 0.0);
		params.addm("dt", .25);
		params.addm("sf int", 2);
		params.addm("init time", 30.0);
		params.addm("quench time", 1000);
		params.add("time");
		params.add("magnetization");
		params.add("Lp");
		flags.add("Clear");
		flags.add("Metric0");
		flags.add("t=0");
	}	
	
	public void animate() {
		params.set("time", format(sim.time()));
		params.set("magnetization", format(sim.magnetization()));
		sim.setParameters(params);
		params.set("Lp", sim.L/dx);
		
		grid.registerData(sim.L/dx, sim.L/dx, sim.getField(dx));
		
		if(flags.contains("Clear")){
			eMetricInverseAcc.clear();
			metric0Acc.clear();
		}
		if(flags.contains("t=0")){
			sett0();
			
		}
		if(flags.contains("Metric0")){
			double metric0 = em.findMetric0();
			System.out.println("Metric0 = " + metric0);
			metric0Acc.accum(sim.time(), metric0);
		}
		flags.clear();

	}
	
	void sett0(){
		sim.resetTime();
		recordStep = 0;
		em.clearSums();
		FileUtil.deleteFile(outputFileName);
		FileUtil.deleteFile(outputFileName0);
//		FileUtil.deleteFile(sfFileNameH);
//		FileUtil.deleteFile(sfFileNameV);
//		FileUtil.deleteFile(sfFileNameDiff);
		FileUtil.initFile(outputFileName, params);
		FileUtil.initFile(outputFileName0, params);
	}
	
	public void clear() {
		eMetricInverseAcc.clear();
		metric0Acc.clear();
	}
	
	public void run() {
		
		initialize();
		recordStep = 0;
		
		double [] sf = new double [sim.L*sim.L];

		
		double step = 0;
		int sfInt = params.iget("sf int");
		boolean initTime = true;
		boolean preQuench = true;
		double maxInitTime = params.fget("init time");
		double quenchTime = params.fget("quench time");
		
		while (true) {
			if(initTime){
				sim.step();
				if(sim.time()>maxInitTime){
						initTime = false;
						sett0();
				}
				Job.animate();
			}else{
				if(preQuench && sim.time()>quenchTime){
					params.set("T", 0.0);
					preQuench = false;
				}
				sim.step();

			em.calculateParticleMetric();
			if (sim.time() > recordStep){
				em.calculateSiteEnergies();
//				em.calculateSiteEnergies();
				double pt = 1.0/em.eMetric;
				double met0 = em.findParticleMetric0();
				FileUtil.printlnToFile(outputFileName, sim.time(), pt);
				FileUtil.printlnToFile(outputFileName0, sim.time(), met0);
				FileUtil.printlnToFile(energyFileName, sim.time(), em.calculateEnergy());
				if(writeSFFiles){
					sf = ft.calculate2DSF(sim.getField(1), false, false);
					double hpoint = 0.5*(sf[sfInt] +sf[(sim.L-sfInt)]);
					double vpoint = 0.5*(sf[sfInt*sim.L] +sf[(sim.L-sfInt)*sim.L]);
					double diff = (hpoint - vpoint);
					FileUtil.printlnToFile(sfFileNameH, sim.time(), hpoint);
					FileUtil.printlnToFile(sfFileNameV, sim.time(), vpoint);
					FileUtil.printlnToFile(sfFileNameDiff, sim.time(), diff);
				}
				recordStep += step;
			}
			Job.animate();
			}
		}
	}

	public void writeToFile(String fileName){
		FileUtil.initFile(fileName, params);
		
	}
	
	public void initialize(){
		sim = new IsingLR(params);
		//sim.setField(params.fget("Initial magnetization"));
		sim.randomizeField(params.fget("Initial magnetization"));		
		dx = 1;
		clearFile = true;
		outputFileName = params.sget("Output Directory") + File.separator + "T" + sim.T + ".txt";
		outputFileName0 = params.sget("Output Directory") + File.separator + "0T" + sim.T + ".txt";
		FileUtil.initFile(outputFileName, params);
		FileUtil.initFile(outputFileName0, params);		
		if(writeSFFiles){
			sfFileNameH = params.sget("Output Directory") + File.separator + "T" + sim.T + "h.txt";
			sfFileNameV = params.sget("Output Directory") + File.separator + "T" + sim.T + "v.txt";
			sfFileNameDiff = params.sget("Output Directory") + File.separator + "T" + sim.T + "d.txt";
			FileUtil.initFile(sfFileNameH, params, "#Horizontal structure factor");			
			FileUtil.initFile(sfFileNameV, params, "#Vertical structure factor");	
			FileUtil.initFile(sfFileNameDiff, params, "#Abs value of difference of hor & vert");
		}
		energyFileName = params.sget("Output Directory") + File.separator + "T" + sim.T + "e.txt";
		FileUtil.initFile(energyFileName, params, "#Energy");
			
		ft = new FourierTransformer(sim.L);
		
		sim.track = true;
		em = new EnergyMetric(sim, params);
		sim.setRangeOffset(params.iget("R offset"));
		
		//try to encourage hor stripes:  only use for M=0
		for (int x = 0; x < sim.L; x++){
			for(int y = 0; y < 20; y++)
				if (sim.spins.get(x, y)==1) sim.spins.flip(x, y);
		}
		for (int x = 0; x < sim.L; x++){
			for(int y = 20; y < 40; y++)
				if (sim.spins.get(x, y)==-1) sim.spins.flip(x, y);
		}

		sim.initTrackParticles();		
		em.initTrackEnergyMetric(sim.noParticles);
		
	}
	
	public void writeConfigToFile(){
		String configFileName = "../../../research/javaData/stripeToClumpInvestigation/monteCarloData/configs/config";
		FileUtil.deleteFile(configFileName);
		FileUtil.writeConfigToFile(configFileName, (sim.L/dx)*(sim.L/dx), sim.getField(dx));
	}
	

}
