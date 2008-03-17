package rachele.ising.testCode.apps;

import java.awt.Color;

import rachele.ising.testCode.ConjugateGradientMin;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
//import static scikit.util.Utilities.asList;
import static scikit.util.Utilities.frame;

public class ConjugateGradientMinApp extends Simulation{
	Grid function = new Grid ("Function");
	Accumulator acc = new Accumulator(1.0);
	ConjugateGradientMin cjMin;
	public int size = 21;
	public double conjGradOutput [] = new double [4];

	public ConjugateGradientMinApp() {
		frame(function);
		params.addm("Minimization", new ChoiceValue("Conjugate Gradient", "SteepestDecent"));
	}
	
	public static void main(String[] args) {
		new Control(new ConjugateGradientMinApp(), "Conjugate Grad Min");
	}

	public void animate() {
		//plot the function:
		double [] point = new double [2];
		double [] funcData = new double [size*size];
		int index;
		for (int x = 0; x < size; x++){
			for(int y = 0; y < size; y++){
				point[0] = x;
				point[1] = y;
				index = size*(y)+(x);
				funcData[index] = cjMin.functionMultiDim(point);
			}
		}
		function.registerData(size, size, funcData);
		acc.clear();
		
		for (int i = 0; i < cjMin.conjIterations; i++){
			function.addDrawable(Geom2D.line(cjMin.drawablePoints[i*cjMin.N]/size, cjMin.drawablePoints[i*cjMin.N+1]/size, cjMin.drawablePoints[(i+1)*cjMin.N]/size, cjMin.drawablePoints[(i+1)*cjMin.N+1]/size, Color.GREEN));
		}
	}	

	public void clear() {
	}

	public void run() {
		cjMin = new ConjugateGradientMin();
		System.out.println(cjMin.function(1.0));
		
		// Give some initial point for the conjugate grad minimization
		double [] initialPoint = new double[cjMin.N];
		//for(int i = 0; i < cjMin.N; i++)
		initialPoint[0] = size*0.0;
		initialPoint[1] = size*0.1;
		if (params.sget("Minimization") == "Conjugate Gradient"){
			cjMin.conjuageGradMin(initialPoint);
		}else{
			cjMin.steepestDecent(initialPoint);
		}
		Job.animate();
		acc.clear();
	}
	
	public Accumulator getAcc(){
		return acc;
	}
	
}