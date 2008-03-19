package chris.tests;



import java.awt.Color;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.FileValue;
import scikit.jobs.params.StringValue;

public class CStest extends Simulation {
	/**
	 * so this seems to be the default starting point HOWEVER it must be created
	 * with three (3) things: (1) public void animate() (2) public void clear()
	 * (3) public void run() Thus, we create them immediately. Note that they go
	 * INSIDE the Simulation class. This can be short cutted by cmnd-click on
	 * the class name (e.g. CStest)
	 */

	// Define some variables: (we imported the Grid class)
	Grid grid1 = new Grid("Checker Board 1");
	Grid grid2 = new Grid("Checker Board 2");
	Grid grid3 = new Grid("Checker Board 3");
	public int L1, L2, L3, M;
	double T1;
	double T2;
	String word;
	public int s1[], s2[], s3[];
	public double s4[];
	public static final double PI = 4 * Math.atan(1);

	/**
	 * Now, to do anything we need a main class.
	 */

	public static void main(String[] args) {
		new Control(new CStest(), "CS First Sim Try");
	}

	/**
	 * So the main serves as a launcher. Now, we need to set up the GUI for how
	 * we want it to look when the program is opened. This is what the LOAD
	 * method is for.
	 */

	public void load(Control c) {

		/**
		 * Put two frames together and leave the third in a separate window.
		 */

		c.frameTogether("First Two Grids", grid1, grid2);
		c.frame(grid3);
		// etc.

		/**
		 * Let's define some parameters (which will in the GUI) for defining
		 * some stuff.
		 */

		// Simple Parameters (must be integers)
		params.add("L1", 1 << 8);
		params.add("L2", 1 << 3);
		// Quick way of 2^n = 1<<n (where << is a binary shift)
		params.add("L3", 1 << 4);

		// More Simple Parameters (floats / doubles)

		params.add("Float", 8.);

		// Parameters with Restrictions

		params.add("T1", new DoubleValue(7, 0, 14).withSlider());
		params.addm("T2", new DoubleValue(50, 1, 100).withSlider());

		/**
		 * params.add CANNOT be changed while app is running whereas params.addm
		 * CAN be changed while app is running
		 */

		// Parameter that cannot be initialized
		params.add("User cannot adjust");

		// Drop down menus

		params.addm("Drop Down", new ChoiceValue("Choose", "Between", "These",
				"Values"));
		params.add("Dop Down 2", new ChoiceValue("This is", "params.add",
				"above is", "params.addm"));

		// String

		params.addm("String", new StringValue("hello wolrd!"));

		// Directory

		// Directory path must be valid else it defaults to
		// /Users/$USER
		params.add("Dir", new DirectoryValue("/Users/cserino/Desktop"));

		// File

		// File path must be valid else it defaults to
		// /Users/$USER
		params.add("File", new FileValue(""));

		/**
		 * So, that's about it for menu options. At least for now.
		 */

	}

	public void animate() {
		/**
		 * In here we do all the visual updates.
		 */

		int i;

		params.set("User cannot adjust", T1);

		// assign each value a color (discrete values)

		ColorPalette palette = new ColorPalette();
		palette.setColor(1, Color.BLACK);
		palette.setColor(-1, Color.WHITE);
		palette.setColor(2, Color.RED);
		palette.setColor(-2, Color.GREEN);
		palette.setColor(3, Color.BLUE);
		palette.setColor(-3, Color.YELLOW);

		// assign each value a color (continuous values)

		ColorGradient smooth = new ColorGradient();

		for (i = 0; i < M; i++) {
			smooth.getColor(s4[i], -1., 1.);
		}

		// apply color assignments to each grid

		grid1.setColors(palette);
		grid2.setColors(palette);
		// grid3.setColors(palette);
		grid3.setColors(smooth);

		// plot the data on the grid

		grid1.registerData(L1, L1, s1);
		grid2.registerData(L1, L1, s2);
		grid3.registerData(L1, L1, s4);

		/**
		 * This if statement sets the value of T1 outside of its allowed range
		 * crashing the programs
		 * 
		 * if (T1>15){ // this should cause problems params.set("T1",T1); }
		 */
	}

	public void clear() {
		/**
		 * Do this when RESET button is clicked.
		 */
		grid1.clear();
		grid2.clear();
		grid3.clear();
	}

	public void run() {

		int i, j;

		/**
		 * Set the value of T1 by retrieving the value from the menu.
		 */
		T1 = params.fget("T1");

		/**
		 * Ends with Job.animate() which calls animate() (make sure to import
		 * Job class).
		 * 
		 * Put everything inside an infinite loop so you can use step and stop
		 * as designed.
		 */

		L1 = (int) (params.fget("L1"));
		M = L1 * L1;

		s1 = new int[M];
		s2 = new int[M];
		s3 = new int[M];
		s4 = new double[M];

		for (i = 0; i < M; i++) {
			s1[i] = 1 * (2 * (i % 2) - 1);
			s2[i] = 2 * (2 * (i % 2) - 1);
			s3[i] = 3 * (2 * (i % 2) - 1);
		}

		while (true) {
			j = (int) (M * Math.random());
			s1[j] *= (-1);
			j = (int) (M * Math.random());
			s2[j] *= -1;
			j = (int) (M * Math.random());
			s3[j] *= -1;

			T1 += 1.;

			T2 = params.fget("T2");

			for (i = 0; i < M; i++) {
				s4[i] = Math.sin((i % L1) * PI * 2 / L1 + (int) (i / L1) * 2
						* PI / L1 + T1 / T2);
			}

			Job.animate();
		}
	}

}
