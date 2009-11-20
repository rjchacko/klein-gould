package chris.tests;


import java.text.DecimalFormat;
import java.text.Format;
import chris.util.Random;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class RandomTest extends Simulation {
	public static void main(String[] args) {
		new Control(new RandomTest(), "Random Clones App");
	}

	Format df  = new DecimalFormat("0.000000000000");
	Format ifr = new DecimalFormat("000");
	double r, rc, g, gc;
	int i, ic;
	Random rand;

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
		
		rand = new Random(params.iget("Seed"));
		
		while (true) {
			for (int jj = 0 ; jj < 25 ; jj++) 
				rand.nextDouble();
			Random randClone = rand.clone(); //= new Random();//
			r  = rand.nextDouble();
			rc = randClone.nextDouble();
			g  = rand.nextGaussian();
			gc = randClone.nextGaussian();
			i  = rand.nextInt(100);
			ic = randClone.nextInt(100);
			Job.animate();
		}
	}

}
