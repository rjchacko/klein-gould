package chris.tests;

import java.awt.Color;

import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import chris.util.LatticeNeighbors;

public class TestCGnbs extends Simulation {

	Grid grid1 = new Grid("Centers");
	Grid grid2 = new Grid("CG Blocks");
	public int L, N, Lp, M;
	public int s[], cntr[];	
	public LatticeNeighbors nbs;


	
	public static void main(String[] args) {
		new Control(new TestCGnbs(), "An App to Test Coarse Graining");
	}
	
	public void load(Control c) {

		c.frameTogether("Plots", grid1, grid2);
		params.add("L", 10);
		params.add("L'", 5);
		
		return;
	}
	
	public void animate() {
		ColorPalette palette = new ColorPalette();
		palette.setColor(0,Color.WHITE);
		palette.setColor(1,Color.BLUE);
		palette.setColor(2,Color.RED);
		palette.setColor(3,Color.GREEN);
		palette.setColor(4,Color.YELLOW);
		palette.setColor(6,Color.GRAY);
		palette.setColor(7,Color.PINK);
		palette.setColor(-1,Color.BLACK);


		
		grid1.setColors(palette);
		grid2.setColors(palette);

		grid1.registerData(L, L, cntr);
		grid2.registerData(L, L, s);
		
	}

	public void clear() {
		
		grid1.clear();
		grid2.clear();
	}

	public void run() {
		
		while(true){
			L    = params.iget("L");
			Lp   = params.iget("L'");
			N    = L*L;
			M    = L/(Lp);
			s    = new int[N];
			cntr = new int[N];

			nbs = new LatticeNeighbors(L,L,0,(Lp-1)/2,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Square);
		
			for (int kk = 0 ; kk < M ; kk++){
				for (int jj = 0 ; jj < M ; jj++){
					int site = (Lp-1)/2+jj*Lp+kk*(Lp*L);
					cntr[site] = -1;
					System.out.println(site);
					int[] tn = nbs.get(site);
					for (int ll = 0 ; ll < tn.length ; ll++){
						s[tn[ll]] = (jj+M*kk)%6+1;
					}
					Job.animate();
				}
			}
			clear();
			Job.animate();
		}
		
	}

}
