package rachele.ising.dim2.apps;

import static java.lang.Math.exp;
import static java.lang.Math.floor;
import static java.lang.Math.pow;
import static scikit.numerics.Math2.j1;
import rachele.ising.dim2.IsingField2D;
import rachele.ising.dim2.IsingLangevin;
import rachele.ising.dim2.StructureFactor;
import rachele.util.FileUtil;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;

public class FindDisorderOrderSpinodalApp  extends Simulation{
	StructureFactor sf;
    IsingField2D ising;

    
	public static void main(String[] args) {
		new Control(new FindDisorderOrderSpinodalApp(), "Ising Field");
	}

	public void load(Control c) {
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("Interaction", new ChoiceValue("Circle", "Square"));
		params.addm("Noise", new ChoiceValue("Off","On"));
		params.addm("Dynamics?", new ChoiceValue("Langevin No M Conservation", "Langevin Conserve M"));
		params.add("Init Conditions", new ChoiceValue("Random Gaussian", 
				"Artificial Stripe 3", "Read From File", "Constant" ));
		params.addm("Approx", new ChoiceValue("Exact Stable",
				"Avoid Boundaries", "Exact SemiStable", "Exact", "Linear",  "Phi4"));
		params.addm("T", 0.05);
		params.addm("H", 0.5);
		params.add("Magnetization", 0.0);
		params.addm("dt", 1.0);
		params.addm("J", -1.0);
		params.addm("R", 1000000.0);
		params.add("L/R", 4.0);
		params.add("R/dx", 16.0);
		params.add("kR bin-width", 0.1);
		params.add("Random seed", 0);
		params.add("Max Time", 350);
		params.add("Time");
		params.add("Reps");
		params.add("Mean Phi");
		
	}
	
	public void animate() {
		ising.readParams(params);
			
	}

	public void clear() {
	}

	public void run() {
		ising = new IsingField2D(params);
		double binWidth = params.fget("kR bin-width");
		binWidth = IsingLangevin.KR_SP / floor(IsingLangevin.KR_SP/binWidth);
        sf = new StructureFactor(ising.Lp, ising.L, ising.R, binWidth, ising.dt);
		sf.setBounds(0.1, 14);	
		while (true) {
			for (double h = -.95; h < .95; h = h + .05){
				System.out.println("h = " + h);
				ising.H = h;
				boolean downDir = true;
				ising.T = .25;
				double dT = .01;
				double upLim = ising.T;
				double dwLim = 0;
				while(upLim - dwLim > .001){
					ising.randomizeField(ising.DENSITY);
					for(int i = 0; i < 100; i ++){
						//System.out.println(i);
						//params.set("Time",ising.time());
						ising.simulate();
					}
					double sf1 = sf.getSF(ising.phi);
					for (double t = 0.0; t < 2; t = t + params.fget("dt"))
						ising.simulate();	
					double sf2 = sf.getSF(ising.phi);
					double di = sf2-sf1;
					System.out.println(di);
					if(sf2-sf1 < 0 || sf2 < pow(10,-15)){
						//above spinodal
						if(downDir==false){
							upLim = ising.T;
							dT /= 10;
							downDir=true;
							//System.out.println("upLim = " + upLim + " at " + ising.T);
						}
						ising.T -= dT;	
						
					}else{
						//under spinodal
						if(downDir==true){
							dwLim = ising.T;
							dT /= 10;
							downDir=false;
							//System.out.println("dwnLim = " + dwLim + " at " + ising.T);
						}
						ising.T += dT;
					}
				}	
				System.out.println("DONE:  Spinodal Temp  for h " +h + " = "+ising.T);
				writeDataToFile(ising.H, ising.T);
			}
		}

	}
	
	public void writeDataToFile(double a, double b){
		String file = "../../../research/javaData/sfData/DOSTvHcircle";
		FileUtil.printlnToFile(file, a, b);
	}
	
	public double findMeanPhi(){
		for (double t = 0.0; t < 50.0; t = t + params.fget("dt"))
		ising.simulate();
		return ising.mean(ising.phi);
	}
	
	public double linearTheory(double kR, double density, double time){
		//double V = ising.Lp*ising.Lp;
		double D = -circlePotential(kR) - ising.T/ (1-pow(ising.H,2));
		//double sf = (exp(2*time*D)*(V + ising.T/D)-ising.T/D)/V;
		double sf = exp(2*time*D);
		return sf;
	}
	

	
	public double circlePotential(double kR){
		return 2*j1(kR)/kR;
	}

}
