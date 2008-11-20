package rachele.ising.dim2.apps;

import rachele.ising.dim2.IsingArbitraryConnect;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Simulation;

public class IsingArbitraryConnectApp extends Simulation{
    IsingArbitraryConnect sim;
	Grid grid = new Grid("Phi(x)");

	public static void main(String[] args) {
		new Control(new IsingArbitraryConnectApp(), "Ising");
	}
	
	public void load(Control c) {
		c.frame(grid);
		
	}
    
	public void animate() {
		grid.registerData(sim.L, sim.L, sim.spin);
		
	}


	public void clear() {

		
	}



	public void run() {
		
        while (true) {
        	sim.step();
        }
		
	}

}
