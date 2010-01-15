package chris.MD.TwoD;

import scikit.jobs.params.Parameters;
import chris.util.vector2d;

public abstract class InteractingSystem {
	
	// this class should be in charge of evolving the system
	// this should be a system
	// it has particles, temperature, etc

	public int N, Lx, Ly;
	public String BC;
	private particle phase[];
	private double dt, t;
	public static enum Solvers {EULER, VERLET};
	private Solvers solver;
	
	public InteractingSystem(Parameters params){
		
		IS_constructor(params);
	}
	
	public void IS_constructor(Parameters params){
	
		t = 0;
	}
	
	
	/**
	 * Abstract method that will return central potential for
	 * two particles separated by a scalar distance r.
	 * 
	 * @param r the distance between the two interacting particles
	 * @return the value of the potential between particles one and two
	 */
	public abstract double potential(double r);
	
	/**
	 * Abstract method that will return central force for
	 * two particles separated by relative coordinates x and y.
	 * 
	 * @param x the relative separation is the x_1 direction
	 * @param y the relative separation is the x_1 direction
	 * @return a vector2d that is the force on particle one due to the second
	 */
	public abstract vector2d force(double x, double y);

	public void step(){
		
		t += dt;
		
		switch (solver) {
		
			case EULER:
				
				vector2d[] accelE = new vector2d[N];

				for(int jj = 0 ; jj < N ; jj++){
					// update the position x(t+dt) <-- x(t) + v(t)*dt 
					phase[jj].r.setX(phase[jj].r.getX()+phase[jj].v.getX()*dt); 
					phase[jj].r.setY(phase[jj].r.getY()+phase[jj].v.getY()*dt); 
					
					// update the velocity v(t+dt) <-- v(t) + a(t)*dt
					phase[jj].v.setX(phase[jj].v.getX()+0.5*accelE[jj].getX()*dt);
					phase[jj].v.setY(phase[jj].v.getY()+0.5*accelE[jj].getY()*dt);
				}
			break;
			
			case VERLET:

				vector2d[] accelV = new vector2d[N];
				double dt2       = dt*dt;
				
				for(int jj = 0 ; jj < N ; jj++){
						// update the position x(t+dt) <-- x(t) + v(t)*dt + 0.5*a(t)*dt^2
						phase[jj].r.setX(phase[jj].r.getX()+phase[jj].v.getX()*dt+0.5*accelV[jj].getX()*dt2); 
						phase[jj].r.setY(phase[jj].r.getY()+phase[jj].v.getY()*dt+0.5*accelV[jj].getY()*dt2); 
						
						// calculate midpoint velocity v(t+dt/2) <-- v(t) + 0.5*a(t)*dt
						phase[jj].v.setX(phase[jj].v.getX()+0.5*accelV[jj].getX()*dt);
						phase[jj].v.setY(phase[jj].v.getY()+0.5*accelV[jj].getY()*dt);
				}
				// get the new accelerations and store them in accel
				/*
				 * 
				 * DO THIS !!!!
				 * 
				 */
				for(int jj = 0 ; jj < N ; jj++){					
					// update the velocity  v(t+dt) <-- v(t+dt/2) + 0.5*a(t+dt)*dt
					phase[jj].v.setX(phase[jj].v.getX()+0.5*accelV[jj].getX()*dt);
					phase[jj].v.setY(phase[jj].v.getY()+0.5*accelV[jj].getY()*dt);
				}
			break;
								
			default:
				throw new IllegalStateException("Solver does not exist.");
		}
		
	}
		
		
}
	
	
	

