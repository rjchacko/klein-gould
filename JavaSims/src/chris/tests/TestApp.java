package chris.tests;


import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class TestApp extends Simulation {
	public static void main(String[] args) {
		new Control(new TestApp(), "TestApp");
	}

	double t;

	public TestApp() {
		params.add("Number", 0.1);
	}

	public void animate() {
		params.set("Number", t);
	}

	public void clear() {
	}

	public void run() {
		t = params.fget("Number");
		while (true) {
			t += 0.0001;
			Job.animate();
		}
	}

}
