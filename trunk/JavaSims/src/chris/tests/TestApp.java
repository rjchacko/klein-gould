package chris.tests;


import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.Parameters;
import chris.util.FileDropBox;
import chris.util.PrintUtil;
import chris.util.ReadInUtil;

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
		Parameters p   = new Parameters();
		ReadInUtil riu = new ReadInUtil("/Users/cserino/Documents/Catalogue/Params_run1.log"); 
		riu.getOFCparams(p);
		PrintUtil.printlnToFile("/Users/cserino/Desktop/test.txt", p.toString());
//		Data Directory =  /Users/cserino/Documents/Catalogue
//		Data File =  run1
//		Random Seed =  0
//		Interaction Shape =  Circle
//		Interaction Radius (R) =  20
//		Lattice Size =  512
//		Boundary Condtions =  Open
//		Equil Time =  1000000
//		Sim Time =  0
//		Failure Stress (?_f) =  2.0
//		?_f width =  0.0
//		Residual Stress (?_r) =  1.0
//		?_r width =  0.025
//		Dissipation (?) =  0.1
//		? width =  0.0
//		Status =  Intializing
//
//		Launch Date & Time : 2011/03/03 15:28:35
	}
//		t = params.fget("Number");
//		System.out.println("Hello World!");
//		while (t < 10) {
//			t += 1;
//			Job.animate();
//		}
//
//		db = new FileDropBox();
//		(new Thread(new FileProducer(db))).start();
//		(new Thread(new FileConsumer(db))).start();		
//
//		
//		wait4file();
//		
////		System.out.println("In TestApp with File: \n \t" + fc.getFile());
//		
//		while (t < 20) {
//			t += 1;
//			Job.animate();
//		}
//		
//		return;
//		
//	}
//	
//	private synchronized void wait4file(){
//		
//		while(!db.isDone()){
//			try {
//				wait(5000); 
//			} catch (InterruptedException e) {
//				e.printStackTrace();
//			}
//			System.out.println(!db.isDone());
//		}
//		return;
//	}
//	

}
