package chris.tests;


import java.text.DecimalFormat;
import java.text.Format;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import chris.ofcdamage.ofc2Dfast;
import chris.util.dummyParamUtil;

public class RandomTest2 extends Simulation {
	public static void main(String[] args) {
		new Control(new RandomTest2(), "Random Clones App");
	}

	Format df  = new DecimalFormat("0.000000000000");
	Format ifr = new DecimalFormat("000");
	double r, rc, g, gc;
	int i, ic;
	ofc2Dfast model1, model2;

	public void load(Control c) {
		params.add("Seed", 0);
		params.add("Orig. Val");
		params.add("Cloned Val");
		params.add("Orig. Val (gauss)");
		params.add("Cloned Val (gauss)");
		params.add("Orig. Val (int)");
		params.add("Cloned Val (int)");
	}

	public void animate() {
		params.set("Orig. Val",df.format(r));
		params.set("Cloned Val",df.format(rc));
		params.set("Orig. Val (gauss)",df.format(g));
		params.set("Cloned Val (gauss)",df.format(gc));
		params.set("Orig. Val (int)",ifr.format(i));
		params.set("Cloned Val (int)",ifr.format(ic));
	}

	public void clear() {
	}

	public void run() {
		
		model1 = new ofc2Dfast(dummyParamUtil.ofcParams(params));
		model2 = new ofc2Dfast(dummyParamUtil.ofcParams(params));

		
		while (true) {
			for (int jj = 0 ; jj < 25 ; jj++) 
				model1.getRand().nextDouble();
			model2.setRandToClone(model1.getRand());
//			for (int jj = 0 ; jj < 25 ; jj++){
//				model1.getRand().nextDouble();
//				model1.getRand().nextGaussian();
//				model1.getRand().nextInt(100);
//			}	
			r  = model1.getRand().nextDouble();
			rc = model2.getRand().nextDouble();
			g  = model1.getRand().nextGaussian();
			gc = model2.getRand().nextGaussian();
			i  = model1.getRand().nextInt(100);
			ic = model2.getRand().nextInt(100);
			Job.animate();
		}
	}

}
