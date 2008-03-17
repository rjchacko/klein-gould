package rachele.ising.dim2;

import scikit.dataset.Accumulator;
import scikit.dataset.PointSet;

abstract public class ConjugateGradientMin {
	public int N;
	public double t;
	public double[] point;
	public double[] oldPoint;
	public double freeEnergy;
	private double f_p;
	private double [] g, h;
	public double [] xi;
	private double horizontalSlice;
	private double verticalSlice;
	private double dx;
	private int Lp;
	
	static final double tolerance = 1e-16;
	static final double EPS = 1e-16;
	
	public abstract double freeEnergyCalc(double[] point);
	public abstract double[] steepestAscentCalc(double[] point);
	private LineMin linemin;

	public boolean gotoIsing = false;
	
	public ConjugateGradientMin(double[] point, double horizontalSlice, double verticalSlice, double dx) {
		N = point.length;
		Lp =  (int)Math.sqrt((double)N);
		this.point = point;
		this.horizontalSlice = horizontalSlice;
		this.verticalSlice = verticalSlice;
		this.dx = dx;
		xi = new double [N];
		g = new double [N];
		h = new double [N];
		oldPoint = new double[N];
		for (int i = 0; i < N; i++)	oldPoint[i] = point[i];
	}
	
	public void initialize(){
		// Find the initial direction of steepest Descent
		f_p = freeEnergyCalc(point);
		xi = steepestAscentCalc(point);
    	for(int j = 0; j < N; j ++){
    		g[j] = -xi[j];
    		h[j] = g[j];
    		xi[j] = h[j];
    	}
		for (int i = 0; i < N; i ++){
			xi[i] *= -1;
		}		
	}
	
	public void step(double t) {
		linemin = new LineMin(point, xi) {
			double freeEnergyCalc(double[] point) {
				return ConjugateGradientMin.this.freeEnergyCalc(point);
			}
			double [] steepestAscentCalc(double [] point){
				return ConjugateGradientMin.this.steepestAscentCalc(point);
			}
		};
//		if(linemin.gotoLinemin == false){
//			gotoIsing = true;
//			return;
//		}

		double fret = linemin.minValue; 
    	//the following lines accept the move:
    	for (int i = 0; i < N; i ++){
    		point [i] = linemin.lineminPoint[i];
    	}

    	xi = steepestAscentCalc(point);
		// Check for doneness
    	//double sub = fret - f_p;
    	double test = tolerance*(Math.abs(fret)+ Math.abs(f_p)+ EPS);
    	//System.out.println("check CG " + sub + " " + test);
    	if(2.0*Math.abs(fret - f_p) <= test){
    		//we are done -> return
    		System.out.println("Conj Grad Min finished at time " + t + " iterations");
    		//return;
    	}//else{
    		//accept the new value of the function value
    		f_p = fret;
    		//Construct the new direction h
    		double dgg = 0.0; //  numeration of gamma scalar = varies by method
    		double gg = 0.0; // denominator of gamma scalar = g_i dot g_i
    		for(int j = 0; j < N; j ++){
    			gg += g[j]*g[j];
    			//dgg += xi[j]*xi[j];			// This statement for Fletcher-Reeves
    			dgg += (xi[j] + g[j])*xi[j];	//This statement for Polak--Ribiere
    		}
    		if(gg == 0.0){
    			System.out.println("Conj Grad Min finished gg = 0 at time " + t + " iterations");
//    			return;
//    			//if gradient is exactly zero, then we are already done
    		}	
    		double gamma = dgg/gg;
    		for(int j = 0; j < N; j++){
    			g[j] = -xi[j];
    			h[j] = g[j] +gamma*h[j];
    			xi[j] = h[j];	//This is our new direction
    		}
    	//}

		freeEnergy = linemin.minValue;
	}
	
	public Accumulator getLandscape(){
		return linemin.getLandscape();
	}
	
    public Accumulator getBracketLandscape(){
    	return linemin.getBracketLandscape();
    }
	
	public PointSet get_delHslice(){
		int y = (int) (horizontalSlice * Lp);
		double slice[] = new double[Lp];
		for (int x = 0; x < Lp; x++) {
			slice[x] = xi[Lp*y + x];
		}
		return new PointSet(0, dx, slice);
	}	

	public PointSet get_delVslice(){
		int x = (int) (verticalSlice * Lp);
		double slice[] = new double[Lp];
		for (int y = 0; y < Lp; y++) {
			slice[y] = xi[Lp*y + x];
		}
		return new PointSet(0, dx, slice);
	}
	
}
