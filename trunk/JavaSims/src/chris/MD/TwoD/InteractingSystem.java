package chris.MD.TwoD;

import chris.util.vector2d;

public abstract class InteractingSystem {
	
	// this class should be in charge of evolving the system

	public int N, Lx, Ly;
	public String BC;
	private particle phase[];
	private double dt;
	public static enum Solvers {EULER, VERLET, NOotherOPTIONS};
	private Solvers solver;
	
	
	public abstract double potential(double r);
	
	public abstract vector2d force(double r);

	public void step(){
		
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
	
	
	

