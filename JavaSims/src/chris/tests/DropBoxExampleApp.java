package chris.tests;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class DropBoxExampleApp extends Simulation{

	int c;
	Drop drop;
	
	public static void main(String[] args) {
		new Control(new DropBoxExampleApp(), "Drop Box Example").getJob().throttleAnimation(true);
	}
	
	public void animate() {
		
		params.set("counter",c);
	}

	public void clear() {
		
	}

	public void load(Control c) {
		params.add("counter");
	}

	public void run() {
     drop = new Drop();
     c          = 0;
      
      while(c++ < 5)
    	  Job.animate();
            
      (new Thread(new Producer(drop))).start();
      (new Thread(new Consumer(drop))).start();	
      
      waiting();

      while(c++ < 10)
    	  Job.animate();
      
	}
	
	public synchronized void waiting(){
	      while(!drop.isDone()){
	    	  try {
	    		  wait(1000);
	    	  } catch (InterruptedException e) {
	    		  // TODO Auto-generated catch block
	    		  e.printStackTrace();
	    	  }
	    	  System.out.println(drop.isDone());
	      }
	}

}
