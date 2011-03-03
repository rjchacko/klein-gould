package chris.tests;

import java.awt.Color;

import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import chris.util.LatticeNeighbors;

public class TestNeighbors2 extends Simulation {

	Grid grid1 = new Grid("Circle");
	Grid grid2 = new Grid("Square");
	Grid grid3 = new Grid("Diamond");
	public int Lx, Ly, N, Rmin, Rmax;
	public int s1[], s2[], s3[];	
	public LatticeNeighbors neighbors1;
	public LatticeNeighbors neighbors2;
	public LatticeNeighbors neighbors3;

	
	public static void main(String[] args) {
		new Control(new TestNeighbors2(), "An App to Test LatticeNeighbors Class");
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
		
		int ic, st, q; 
		
		Lx   = params.iget("Lx");
		Ly   = params.iget("Ly");
		N    = Lx*Ly;
		Rmin = params.iget("R_min");
		Rmax = params.iget("R_max");
		
		s1 = new int[N];
		s2 = new int[N];
		s3 = new int[N];

		if (N%2 == 1){
			st = (int)(N/2)+1;
			if (params.sget("Cell Location").equals("CORNER")) st = (int)(2*Lx+1);
		}
		else {
			st = (int)(N/2)+(int)(Lx/2);
			if (params.sget("Cell Location").equals("CORNER")) st = (int)(2*Lx+1);
		}
		
		
		if (params.sget("Boundary Condtions").equals("BORDERED")){
			ic = (int)((1 + Lx)*(Ly/2));
			neighbors1 = new LatticeNeighbors(Lx, Ly, Rmin, Rmax, LatticeNeighbors.Type.BORDERED, LatticeNeighbors.Shape.Circle);	
			neighbors2 = new LatticeNeighbors(Lx, Ly, Rmin, Rmax, LatticeNeighbors.Type.BORDERED, LatticeNeighbors.Shape.Square);	
			neighbors3 = new LatticeNeighbors(Lx, Ly, Rmin, Rmax, LatticeNeighbors.Type.BORDERED, LatticeNeighbors.Shape.Diamond);	
		}
		else{
			ic = 0;
			neighbors1 = new LatticeNeighbors(Lx, Ly, Rmin, Rmax, LatticeNeighbors.Type.PERIODIC, LatticeNeighbors.Shape.Circle);
			neighbors2 = new LatticeNeighbors(Lx, Ly, Rmin, Rmax, LatticeNeighbors.Type.PERIODIC, LatticeNeighbors.Shape.Square);
			neighbors3 = new LatticeNeighbors(Lx, Ly, Rmin, Rmax, LatticeNeighbors.Type.PERIODIC, LatticeNeighbors.Shape.Diamond);
		}		
		
		
		int[] nbs1 = neighbors1.get(ic);
		int[] nbs2 = neighbors2.get(ic);
		int[] nbs3 = neighbors3.get(ic);		
		
		for (int i = 0 ; i < nbs1.length ; i++ ) {
			q = neighbors1.getJ(st,ic,nbs1,i);
			if(q>=0) s1[q] = 1;
		}
		
		for (int i = 0 ; i < nbs2.length ; i++ ) {
			q = neighbors2.getJ(st,ic,nbs2,i);
			if(q>=0) s2[q] = 1;
		}
		
		for (int i = 0 ; i < nbs3.length ; i++ ) {
			q = neighbors3.getJ(st,ic,nbs3,i);
			if(q>=0) s3[q] = 1;
		}

		s1[st]=2;
		s2[st]=2;
		s3[st]=2;
		
		while(true){
			Job.animate();
		}
		
	}

}
