package chris.tests;

import java.awt.Color;

import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import chris.util.LatticeNeighbors;

public class TestNeighbors extends Simulation {

	Grid grid1 = new Grid("Circle");
	Grid grid2 = new Grid("Square");
	Grid grid3 = new Grid("Diamond");
	public int Lx, Ly, N, Rmin, Rmax;
	public int s1[], s2[], s3[];	
	public LatticeNeighbors neighbors1;
	public LatticeNeighbors neighbors2;
	public LatticeNeighbors neighbors3;

	
	public static void main(String[] args) {
		new Control(new TestNeighbors(), "An App to Test LatticeNeighbors Class");
	}
	
	public void load(Control c) {

		/**
		 * Put two frames together and leave the third in a separate window.
		 */

		c.frameTogether("Plots", grid1, grid2, grid3);
		params.add("Boundary Condtions", new ChoiceValue("BORDERED","PERIODIC"));
		params.add("Cell Location" ,new ChoiceValue("CENTER","CORNER"));
		params.add("Lx", 101);
		params.add("Ly", 101);
		params.add("R_min", 0);
		params.add("R_max", 20);
		
		
	}
	
	public void animate() {
		ColorPalette palette = new ColorPalette();
		palette.setColor(0, Color.WHITE);
		palette.setColor(1,Color.BLUE);
		palette.setColor(2,Color.RED);
		palette.setColor(3,Color.GREEN);
		
		grid1.setColors(palette);
		grid2.setColors(palette);
		grid3.setColors(palette);

		grid1.registerData(Lx, Ly, s1);
		grid2.registerData(Lx, Ly, s2);
		grid3.registerData(Lx, Ly, s3);
		
	}

	public void clear() {
		
		grid1.clear();
		grid2.clear();
		grid3.clear();
		
	}

	public void run() {
		
		int ic; 
		
		Lx   = params.iget("Lx");
		Ly   = params.iget("Ly");
		N    = Lx*Ly;
		Rmin = params.iget("R_min");
		Rmax = params.iget("R_max");
		
		s1 = new int[N];
		s2 = new int[N];
		s3 = new int[N];
		
		for (int i = 0 ; i < N ; i++ ) {
			s1[i] = 0;
			s2[i] = 0;
			s3[i] = 0;
		}

		if (N%2 == 1){
			ic = (int)(N/2)+1;
			if (params.sget("Cell Location").equals("CORNER")) ic = (int)(0.8*(N-Lx));
		}
		else {
			ic = (int)(N/2)+(int)(Lx/2);
			if (params.sget("Cell Location").equals("CORNER")) ic = (int)(0.8*N+0.25*Lx);
		}
		
		if (params.sget("Boundary Condtions").equals("BORDERED")){
			neighbors1 = new LatticeNeighbors(Lx, Ly, Rmin, Rmax, LatticeNeighbors.Type.BORDERED, LatticeNeighbors.Shape.Circle);	
			neighbors2 = new LatticeNeighbors(Lx, Ly, Rmin, Rmax, LatticeNeighbors.Type.BORDERED, LatticeNeighbors.Shape.Square);	
			neighbors3 = new LatticeNeighbors(Lx, Ly, Rmin, Rmax, LatticeNeighbors.Type.BORDERED, LatticeNeighbors.Shape.Diamond);	
		}
		else{
			neighbors1 = new LatticeNeighbors(Lx, Ly, Rmin, Rmax, LatticeNeighbors.Type.PERIODIC, LatticeNeighbors.Shape.Circle);
			neighbors2 = new LatticeNeighbors(Lx, Ly, Rmin, Rmax, LatticeNeighbors.Type.PERIODIC, LatticeNeighbors.Shape.Square);
			neighbors3 = new LatticeNeighbors(Lx, Ly, Rmin, Rmax, LatticeNeighbors.Type.PERIODIC, LatticeNeighbors.Shape.Diamond);
		}		
		
		
		int[] nbs1 = neighbors1.get(ic);
		int[] nbs2 = neighbors2.get(ic);
		int[] nbs3 = neighbors3.get(ic);		
		
		for (int i = 0 ; i < nbs1.length ; i++ ) {
			s1[nbs1[i]] = 1;
			s2[nbs1[i]] = 3;
			s3[nbs1[i]] = 3;
		}
		
		for (int i = 0 ; i < nbs2.length ; i++ ) {
			s2[nbs2[i]] = 1;
		}
		
		for (int i = 0 ; i < nbs3.length ; i++ ) {
			s3[nbs3[i]] = 1;
		}

		s1[ic]=2;
		s2[ic]=2;
		s3[ic]=2;
		
		while(true){
			Job.animate();
		}
		
	}

}
