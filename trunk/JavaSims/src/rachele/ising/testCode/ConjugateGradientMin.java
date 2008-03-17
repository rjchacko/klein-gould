package rachele.ising.testCode;

public class ConjugateGradientMin {
	public int N = 2;
	public double [] drawablePoints= new double [1000];
	public int conjIterations;
	static double GOLD = 1.618034;
	static double GOLDR = 0.61803399;
	static double GOLDC = 1-GOLDR;
	static double maxIterations = 1000;
	static double tolerance = 1e-16;
	static double EPS = 1e-16;
	public boolean checked = false;
	public double lambda;
	static double delta = 1e-10;
	public double switchD = 1;
	static double ZEPS = 0.0000000001;  //small number that protects against
	//trying to achieve fractional accuracy for a minimum that happens to be
	// exactly zero
		
	public double [] minReferencePoint = new double [N];
	public double [] minReferenceDirection = new double [N];
	public double [] xi = new double [N];
	
	public ConjugateGradientMin(){
		
	}
	
	public double function(double x) {
		double func = Math.pow(x+0.0,4)- 10*Math.pow(x,2);
		//double func = Math.pow(x+3, 2);
		return func;
	}
	/**
	* initial bracket method takes two random points, ax and bx,
	* where f(ax) > f(bx). The point ax is one end point of the 
	* bracketing procedure and bx is the midpoint.  We wish to find
	* another point cx such that f(cx) > f(bx).  These will be our 
	* first 3 input points for the actual minimization procedure. 
	*/
	public double [] initialBracket(double ax, double bx){
		double [] output = new double [3];
		double f_b = function(bx);
		double f_a = function(ax);
		double u, f_u;
		System.out.println("initial a and b points:");
		System.out.println(ax + " " + bx);
		System.out.println("initial a and b values:");
		System.out.println(f_a + " " + f_b);
		
		//Check to see if f(ax) has a higher value than f(bx)
		//If not, switch roles of ax and bx so that we can 
		//go downhill in the direction from ax to bx
		if (f_b > f_a){
			double tempA = f_a;
			f_a = f_b;
			f_b = tempA;
			tempA = ax;
			ax = bx;
			bx = tempA;
			System.out.println("switched: ax = " + ax + " bx = " + bx);			
		}
		
		//First guess for midpoint
		double cx =	bx + GOLD*(bx-ax);
		double f_c = function(cx);
		
		//repeat the following until we bracket
		int iterations = 0;
		while(f_b > f_c){
			iterations ++;
			u = cx + GOLD*(cx-bx);
			f_u = function(u);
			ax = bx;
			f_a = f_b;
			bx = cx;
			f_b = f_c;
			cx = u;
			f_c = f_u;
		}
		
		output[0] = ax;
		output[1] = bx;
		output[2] = cx;
		System.out.println("finished after " + iterations + " iterations");
		return output;
	}
	
	public double [] initialBracketMultiDim(double ax, double bx, double point [], double direction[]){
		double [] output = new double [3];
		
		for (int i = 0; i < N; i++){
			minReferencePoint[i] = point[i];
			minReferenceDirection[i] = direction[i];			
		}

		//make sure bx is positive
		
		bx = Math.abs(bx)*switchD;
		
		double f_b = f1dim(bx);
		double f_a = f1dim(ax);
		double u, f_u;
		
		//Check to see if f(ax) has a higher value than f(bx)
		//If not, make bx smaller and smaller until it does
		
		while (f_b > f_a){
			bx = .5*bx;
			f_b = f1dim(bx);
			if(bx < tolerance){
				System.out.print("bx less than tolerence");
			}
		}
		System.out.println("bx = " + bx);
		
		//If not, switch the sign of bx and check
		
//		if (f_b > f_a){
//			bx = -bx;
//			f_b = f1dim(bx);
//		}
//		
//		//Check to see if f(ax) has a higher value than f(bx)
//		//If not, switch roles of ax and bx so that we can 
//		//go downhill in the direction from ax to bx
//		if (f_b > f_a){
//			double tempA = f_a;
//			f_a = f_b;
//			f_b = tempA;
//			tempA = ax;
//			ax = bx;
//			bx = tempA;
//		}
		
		//First guess for midpoint
		double cx =	bx + GOLD*(bx-ax);
		double f_c = f1dim(cx);
		
		//repeat the following until we bracket
		int iterations = 0;
		while(f_b > f_c){
			iterations ++;
			u = cx + GOLD*(cx-bx);
			f_u = f1dim(u);
			ax = bx;
			f_a = f_b;
			bx = cx;
			f_b = f_c;
			cx = u;
			f_c = f_u;
		}
		
		output[0] = ax;
		output[1] = bx;
		output[2] = cx;
		System.out.println("init bracket finished after " + iterations + " iterations");
		return output;
	}
	
	public double [] goldenMin(double [] input){
		double [] output = new double [2];
		double x0, x1, x2, x3;
		double ax = input[0];
		double bx = input[1];
		double cx = input[2];
		System.out.println("initial points: " + ax + " " + bx + " " + cx);
		x0 = ax;
		x3 = cx;
		
		//put the new test point in the longer segment
		//if bc is the longer segment, then x2 is the 
		//new test point and goes in between b and c
		
		if(Math.abs(cx - bx) > Math.abs(bx - ax)){
			System.out.println("bc larger");
			x1 = bx;
			x2 = bx + GOLDC*(cx - bx);
		}else{
			x2 = bx;
			x1 = bx - GOLDC*(bx - ax);
		}
		
		double f_1 = function(x1);
		double f_2 = function(x2);
		
		//int iteration = 0;
		while(Math.abs(x3-x0) > tolerance*(Math.abs(x1) + Math.abs(x2))){
			if(f_2 < f_1){
				// choose new triplet as x1, x2, x3
				// shift variables x0, x1, x2 (new test point
				// which is chosen to go inside the larger segment
				x0 = x1;
				x1 = x2;
				x2 = GOLDR*x1 + GOLDC*x3; 
				// x2 goes a distance GOLDC*(x3-x1) inside the x3 x1 segment:
				// so x2 = x1 + GOLDC*(x3-x1)
				//       = x1*(1-GOLDC) + GOLDC*x3
				//       = x1*GOLDR + x3*GOLDC
				
				f_1 = f_2;
				f_2 = function(x2);
				
			}else{
				x3 = x2;
				x2 = x1;
				x1 = GOLDR*x2 + GOLDC*x0;
				// x1 = x0 + GOLDC*(x2-x0)
				//    = x0*GOLDR + x2*GOLDC
				
				f_2 = f_1;
				f_1 = function(x1);
			}
		}
		
		double xmin;
		double minValue;
		if(f_1 < f_2){
			xmin = x1;
			minValue = f_1;
		}else{
			xmin = x2;
			minValue = f_2;
		}
		

		output[0] = xmin;
		output[1] = minValue;
		return output;
	}
	
	public double goldenMinMultiDim(double [] input){
		
		double x0, x1, x2, x3;
		double ax = input[0];
		double bx = input[1];
		double cx = input[2];
		x0 = ax;
		x3 = cx;
		
		//put the new test point in the longer segment
		//if bc is the longer segment, then x2 is the 
		//new test point and goes in between b and c
		
		if(Math.abs(cx - bx) > Math.abs(bx - ax)){
			x1 = bx;
			x2 = bx + GOLDC*(cx - bx);
		}else{
			x2 = bx;
			x1 = bx - GOLDC*(bx - ax);
		}
		
		double f_1 = f1dim(x1);
		double f_2 = f1dim(x2);
		
		int iteration = 0;
		while(Math.abs(x3-x0) > tolerance*(Math.abs(x1) + Math.abs(x2))  && x1 != x2 && x2 != x3){
			iteration ++;
//			System.out.println(x1 + " " + x2 + " " + x3);
//			System.out.println(Math.abs(x3-x0) + " " + tolerance*(Math.abs(x1) + Math.abs(x2)));
			
			if(f_2 < f_1){
				// choose new triplet as x1, x2, x3
				// shift variables x0, x1, x2 (new test point
				// which is chosen to go inside the larger segment
				x0 = x1;
				x1 = x2;
				x2 = GOLDR*x1 + GOLDC*x3; 
				// x2 goes a distance GOLDC*(x3-x1) inside the x3 x1 segment:
				// so x2 = x1 + GOLDC*(x3-x1)
				//       = x1*(1-GOLDC) + GOLDC*x3
				//       = x1*GOLDR + x3*GOLDC
				
				f_1 = f_2;
				f_2 = f1dim(x2);
				
			}else{
				x3 = x2;
				x2 = x1;
				x1 = GOLDR*x2 + GOLDC*x0;
				// x1 = x0 + GOLDC*(x2-x0)
				//    = x0*GOLDR + x2*GOLDC
				
				f_2 = f_1;
				f_1 = f1dim(x1);
			}
		}
		
		System.out.println("Golden Min finished after " + iteration + " iterations");
		
		double xmin;
		double minValue;
		if(f_1 < f_2){
			xmin = x1;
			minValue = f_1;
		}else{
			xmin = x2;
			minValue = f_2;
		}
		for(int i = 0; i < N; i ++){
			minReferencePoint[i] = minReferencePoint[i] + xmin*minReferenceDirection[i];
		}
		System.out.println("golden min = " + minReferencePoint[0] + " " + minReferencePoint[1] + " lambda = " + xmin);
		minReferenceDirection = dFunctionMultiDim(minReferencePoint);
		lambda = xmin;
		return minValue;
	}
	
	//This one not working now
    public double [] brent(double [] input){
		double [] output = new double [2];
    	double xm; //midpoint between a and b
    	double e = 0.0; //This will be the distance moved on the step before last.
		double d = 0.0;
    	
		double ax = input[0];
		double bx = input[1];
		double cx = input[2];		
		
		double a = (ax < cx ? ax : cx);  //the lower bracket extreme
		double b = (ax > cx ? ax : cx);  //the higher bracket extreme
		double x = bx;  //the point with the very least function value found so far
						// or the most recent in case of a tie
		double w = bx;  //the point with the second least function value
		double v = bx; //previous value of w
		double u;  //the point at which the function was evaluated most recently
		
		double f_x = function(x);
		double f_w = f_x;
		double f_v = f_x;
		
		for(int iteration = 1; iteration <= maxIterations; iteration ++){
			xm = 0.5*(a + b);
			double tol1 = tolerance*Math.abs(x) + ZEPS;
			double tol2 = 2.0*(tol1);
			//test for done here
			if (Math.abs(x-xm) <= (tol2 - 0.5*(b-a)) ){
				//min x = x
				//minValue = f_x
				output[0] = x;
				output[1] = f_x;    	
				return output;				
			}

			if(Math.abs(e) > tol1){
				//construct a trial parabolic fit for points x, v and w
				//if the distance moved before the last was greater than the tolerence
				double r = (x-w)*(f_x-f_v);
				double q = (x-v)*(f_x-f_w);
				double p = (x-v)*q-(x-w)*r;
				q = 2.0*(q-r);
				if(q>0.0)
					p = -p;
				q = Math.abs(q);
				double eTemp = e;
				e = d;
				if(Math.abs(p) >= Math.abs(0.5*q*eTemp) || p <= q*(a-x) || p >= q*(b-x)){
					// If the parabolic fit is acceptable, just make the golden ratio step
					// into the larger of the two segments 
					e = (x >= xm ? a-x : b-x);
					d = GOLDC*e;;
				}else{
				    //else take the parabolic step
					d = p/q;  // formula for the min of parabola from 3 given points
					u = x+d;
					if (u-a < tol2 || b-u < tol2){
						d = tol1*Math.signum(xm-x);
					}
					//never take a step that is smaller than the tolerance
				}
			}else{		
				//else just make a golden ratio step
				e = (x >= xm ? a-x : b-x);
				d = GOLDC * e;
			}
			u = Math.abs(d) >= tol1 ? x+d : x + Math.signum(d)*tol1;
			// Never take a step that is smaller than tolerance:
			// There is no new info there
			double f_u = function(u);
			// Above line is the one function evaluation per iteration:
			// And now we need to decide what to do with our function evaluation
			// Housekeeping follows
			if(f_u <= f_x){
				if(u < x){
					a = x;
				}else{
					b = x;
				}
				v = w;
				w = u;
				x = u;
				f_v = f_w;
				f_w = f_x;
				f_x = f_u;
			}else{
				if(u < x){
					a = u;
				}else{
					b = u;
				}
				if(f_u <= f_w || w == x){
					v = w;
					w = u;
					f_v = f_w;
					f_w = f_u;
				}else if(f_u <= f_v || v == x || v == w){
					v = u;
					f_v = f_u;
				}
			}
		}//done with housekeeping.  back for another interation
		//some error code here-- too many iterations in brent
		
		// should never get here
		output[0] = x;
		output[1] = f_x;    	
		return output;
    }

    public double linemin(double point[], double direction[]){
    	double [] initBracket = new double [3]; 
     	// Make up two initial configurations and find an initial bracket
    	//be careful about input:  start at lambda_a = 0 amd lambda_b = positive so that it goes downhill 
    	initBracket = initialBracketMultiDim(0, 100.0, point, direction);
    	//Find the min with golden
    	double minValue = goldenMinMultiDim(initBracket);
    	return minValue;
    }
    
    public double freeEnergyCalc(double point){
    	//define 
		//freeEnergy += potential - T*entropy - H*phi[i];
		double ret = 0;
		return ret;
    }
    
    public double functionMultiDim(double point[]){
    	
    	double x = point[0] - 10.0;
    	double y = point[1] - 3.5;
    	double value = -100*x*x+ Math.pow(x,4) + y*y;
    	//double value = Math.pow(x, 2) + 10*Math.pow(y, 2);
    	//double value = 1*Math.pow(Math.sin((x-1.0)*(y-5.0)),2);
     	return value;
    }
    
    public double [] dFunctionMultiDim(double point[]){
    	double [] direction = new double [N];
    	double x = point[0] - 10.0;
    	double y = point[1] - 3.50;
    	direction[0] = -200*x + 4*x*x*x;
    	direction [1] = 2*y;
    	//direction[0] = 2*(x);
    	//direction[1] = 20*(y);
    	//direction[0] = 2*Math.sin((x-1)*(y-5))*Math.cos((x-1)*(y-5))*(y-5);
    	//direction[1] = 2*Math.sin((x-1)*(y-5))*Math.cos((x-1)*(y-5))*(x-1);
    	
    	return direction;
    }
       
    public double f1dim(double lambda){
    	double newPoint [] = new double [N];
    	for(int i = 0; i < N; i++){
    		newPoint[i] = minReferencePoint[i] + lambda*minReferenceDirection[i];
    	}
    	double ret = functionMultiDim(newPoint);
    	return ret;
    }
    
    public double [] conjuageGradMin(double point[]){
    	double g [] = new double[N];
    	double h [] = new double[N];
    	double xi [] = new double[N];
    	double fret;
    	
    	// Initializations:
    	// evaluate function and derivative at given point
    	// set all vectors (xi, g, h) equal to the direction
    	// of steepest decent at this point
    	
    	for (int i = 0; i < N; i ++){
    		drawablePoints[i] = point[i];
    	}    	
    	
    	double f_p = functionMultiDim(point);
    	xi = dFunctionMultiDim(point);
    	for(int j = 0; j < N; j ++){
    		g[j] = -xi[j];
    		h[j] = g[j];
    		xi[j] = h[j];
    	}
    	
    	System.out.println("initial xi= " + xi[0] + " " + xi[1]);    	
    	switchD = 1;
    	
    	for(int iteration = 1; iteration < maxIterations; iteration ++){
    		System.out.println("iteration = " + iteration);
    		conjIterations = iteration;
			fret = linemin(point, xi); 
			//check to see if this move is OK
			checkMove();
        	if(checked == true){
        		//the following lone accepts the move:
        		for (int i = 0; i < N; i ++){
        			//initPoint[i] = point[i];
        			point [i] = minReferencePoint[i];
        			xi[i] = minReferenceDirection[i];
        			drawablePoints[iteration*N + i] = point[i];
        		}
    		
        		// Check for doneness
        		if(2.0*Math.abs(fret - f_p) <= tolerance*(Math.abs(fret)+ Math.abs(f_p)+ EPS)){
        			//we are done -> return
        			System.out.println("Conj Grad Min finished after " + iteration + " iterations");
        			return point;
        		}else{
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
        				System.out.println("Conj Grad Min finished gg = 0 after " + iteration + " iterations");
        				return point;
        				//if gradient is exactly zero, then we are already done
        			}	
        			double gamma = dgg/gg;
        			for(int j = 0; j < N; j++){
        				g[j] = -xi[j];
        				h[j] = g[j] +gamma*h[j];
        				xi[j] = h[j];	//This is our new direction
        			}
        		}
        		switchD = 1;
    		}else{
    			for (int i = 0; i < N; i++)
    				point[i] = minReferencePoint[i];
    			switchD = -1;
    		}
    		
    	}
    System.out.println("Maximum iterations exceeded");
    return point;
    }

    public void checkMove(){
    	double trialPoint [] = new double [N];
    	double trialPointD [] = new double [N];
    	double df = 0;
    	checked = false;
  
    	//before making the move, sample some points along the way to see if you are traversing 
    	//any uphill paths
    	//if(iteration ==1){
    	for (int i = 1; i <10; i++){
    		for (int j = 0; j < N; j++){
    			trialPoint[j] = minReferencePoint[j] + lambda*(double)i*0.1*xi[j];
    			trialPointD[j] = minReferencePoint[j] + lambda*((double)i*0.1 + delta)*xi[j];
			}    				
    		//evaluate the der at the point + lambda*(double)i*.1*xi

    		df = functionMultiDim(trialPointD) - functionMultiDim(trialPoint);
    		System.out.println("df = " + df);
    		if(df > 0.0){
    			for (int j = 0; j < N; j ++){
    				minReferencePoint[j] = trialPoint[j];
    				//xi[j] = -xi[j];
    				return;
    			}
			}
		}
    	checked = true;
    }
    
    
    public double [] steepestDecent(double point []){
    	double xi [] = new double[N];
    	// Initializations:
    	// evaluate function and derivative at given point
    	// set all vectors (xi, g, h) equal to the direction
    	// of steepest decent at this point
    	
    	for (int i = 0; i < N; i ++){
    		drawablePoints[i] = point[i];
    	}    	
    	
    	double f_p = functionMultiDim(point);
    	xi = dFunctionMultiDim(point);
    	for(int i = 0; i < N; i++)
    		xi[i] = -xi[i];
    	
    	
    	for(int iteration = 1; iteration < maxIterations; iteration ++){
    		System.out.println("iteration = " + iteration);
    		conjIterations = iteration;
        	double fret = linemin(point, xi);
        	for (int i = 0; i < N; i ++){
        		point [i] = minReferencePoint[i];
        		xi[i] = -minReferenceDirection[i];
        		drawablePoints[iteration*N + i] = point[i];
        	}
    		
        	// Check for doneness
        	if(2.0*Math.abs(fret - f_p) <= tolerance*(Math.abs(fret)+ Math.abs(f_p)+ EPS)){
    			//we are done -> return
        		System.out.println("Steepest Decent finished after " + iteration + " iterations");
        		return point;
    		}else{
    			//accept the new value of the function value
    			f_p = fret;
    		}
    	}
    System.out.println("Maximum iterations exceeded");
    return point;
    }
    
}
