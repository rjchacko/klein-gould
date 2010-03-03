package chris.tests;


import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import chris.util.FileConsumer;
import chris.util.FileDropBox;
import chris.util.FileProducer;

public class TestApp extends Simulation {
	
	double t;
	FileDropBox db;

	public static void main(String[] args) {
		new Control(new TestApp(), "TestApp").getJob().throttleAnimation(true);
	}
	
	public void load(Control c) {
		params.add("Number", 0.1);
	}

	public void animate() {
		params.set("Number", t);
	}

	public void clear() {
	}

	public void run() {
		t = params.fget("Number");
		System.out.println("Hello World!");
		while (t < 10) {
			t += 1;
			Job.animate();
		}

		db = new FileDropBox();
		(new Thread(new FileProducer(db))).start();
		(new Thread(new FileConsumer(db))).start();		

		
		wait4file();
		
//		System.out.println("In TestApp with File: \n \t" + fc.getFile());
		
		while (t < 20) {
			t += 1;
			Job.animate();
		}
		
		return;
		
	}
	
	private synchronized void wait4file(){
		
		while(!db.isDone()){
			try {
				wait(5000); 
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			System.out.println(!db.isDone());
		}
		return;
	}
	

}
