package chris.MD.TwoD;

import java.awt.Color;

import scikit.graphics.Drawable;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Gfx2D;
import scikit.jobs.params.Parameters;
import scikit.util.Bounds;
import chris.util.Random;
import chris.util.vector2d;

public abstract class InteractingSystem {
	
	// this class should be in charge of evolving the system
	// this should be a system
	// it has particles, temperature, etc

	public static enum Solvers {EULER, VERLET};
	public static enum IC {GIVEN_KE, GIVEN_LATTICE};
	public static enum BC {PERIODIC, CLOSED}
	
	public int N, Lx, Ly;
	public particle phase[];	// MAKE PRIVATE AFTER DEBUGGING
	private double dt, t;
	private Solvers solver;
	private IC initCond;
	private BC boundCond;
	private boolean haveAccel;
	private vector2d accel[];
	private Random rand;
	
	public InteractingSystem(Parameters params){
		
		IS_constructor(params);
	}
	
	public void IS_constructor(Parameters params){
	
		dt        = params.fget("dt");
		t         = 0;
		haveAccel = false;
		N         = params.iget("N");
		Lx        = params.iget("Lx");
		Ly        = params.iget("Ly");  // make more robust later
		accel     = new vector2d[N]; 
		phase     = new particle[N];
		rand      = new Random(params.iget("seed"));
		solver    = Solvers.VERLET;	// for now
		
		if(params.sget("Boundary Conditions").equals("Closed")){
			boundCond = BC.CLOSED;
		}
		else if(params.sget("Boundary Conditions").equals("Periodic")){
			boundCond = BC.PERIODIC;
		}
		
		if(params.sget("Initial Conditions").equals("Kinetic Energy")){
			initCond = IC.GIVEN_KE;
		}
		else if(params.sget("Initial Conditions").equals("Lattice")){
			initCond = IC.GIVEN_LATTICE;
		}
		
		initialise(params);
	}
	
	
	/**
	 * Abstract method that will return central potential for
	 * two particles separated by a scalar distance r.
	 * 
	 * @param r the (scalar) distance between the two interacting particles
	 * @return the value of the potential between particles one and two
	 */
	public abstract double potential(double r);
	
	/**
	 * Abstract method that will return central force for
	 * two particles separated by relative coordinates x and y.
	 * 
	 * @param r1 the vector displacement of particle 1 
	 * @param r2 the vector displacement of particle 2
	 * @return a vector2d that is the force on particle one due to the second
	 */
	public abstract vector2d force(vector2d r1, vector2d r2);
	
	public vector2d[] netForce(){
		
		for (int jj = 0 ; jj < N ; jj++){
			
		}
		
		return null;
	}

	public void step(){

		t += dt;
	
		switch (solver) {
		
			case EULER:
								
				// kk = 0, 1 , 2, ... , jj - 2, jj - 1, jj + 1, jj + 2, ... , N-2, N-1 
				for(int jj = 0 ; jj < N ; jj++){
					accel[jj] = new vector2d(); // Initializes to the zero vector
					for (int kk = 0 ; kk < jj ; kk++){
						accel[jj].add(force(phase[jj].r,phase[kk].r));
					}
					for (int kk = jj+1 ; kk < N ; kk++){
						accel[jj].add(force(phase[jj].r,phase[kk].r));
					}
					accel[jj].scale(1/phase[jj].m);
				}

					switch (boundCond){

						case PERIODIC:
							for(int jj = 0 ; jj < N ; jj++){
								// update the position x(t+dt) <-- x(t) + v(t)*dt 
								phase[jj].r.x = (phase[jj].r.x+phase[jj].v.x*dt + Lx)%Lx;
								phase[jj].r.y = (phase[jj].r.y+phase[jj].v.y*dt + Ly)%Ly;

								// update the velocity v(t+dt) <-- v(t) + a(t)*dt
								phase[jj].v.x = phase[jj].v.x+0.5*accel[jj].x*dt;
								phase[jj].v.y = phase[jj].v.y+0.5*accel[jj].y*dt;
							}
						break;

						case CLOSED:
							for(int jj = 0 ; jj < N ; jj++){
								// update the position x(t+dt) <-- x(t) + v(t)*dt 
								phase[jj].r.x = phase[jj].r.x+phase[jj].v.x*dt;
								phase[jj].r.y = phase[jj].r.y+phase[jj].v.y*dt;
																
								// update the velocity v(t+dt) <-- v(t) + a(t)*dt
								phase[jj].v.x = phase[jj].v.x+0.5*accel[jj].x*dt;
								phase[jj].v.y = phase[jj].v.y+0.5*accel[jj].y*dt;
								
								if(phase[jj].r.x < 0){
									phase[jj].r.x = -phase[jj].r.x;
									phase[jj].v.x = -phase[jj].v.x;
								}
								else if(phase[jj].r.x > Lx){
									phase[jj].r.x = 2*Lx - phase[jj].r.x; 
									phase[jj].v.x = -phase[jj].v.x;
								}
								if(phase[jj].r.y < 0){
									phase[jj].r.y = -phase[jj].r.y;
									phase[jj].v.y = -phase[jj].v.y;
								}
								else if(phase[jj].r.y > Ly){
									phase[jj].r.y = 2*Ly - phase[jj].r.y; 
									phase[jj].v.y = -phase[jj].v.y;
								}
							}
						break;

						default:

						break;
					}
			break;
			
			case VERLET:

				double dt2 = dt*dt;
				
				if(!haveAccel){
					// kk = 0, 1 , 2, ... , jj - 2, jj - 1, jj + 1, jj + 2, ... , N-2, N-1 
					for(int jj = 0 ; jj < N ; jj++){
						accel[jj] = new vector2d(); // Initializes to the zero vector
						for (int kk = 0 ; kk < jj ; kk++){
							accel[jj].add(force(phase[jj].r,phase[kk].r));
						}
						for (int kk = jj+1 ; kk < N ; kk++){
							accel[jj].add(force(phase[jj].r,phase[kk].r));
						}
						accel[jj].scale(1/phase[jj].m);
					}
				}
				
				
				switch (boundCond){

					case PERIODIC:
						for(int jj = 0 ; jj < N ; jj++){
							// update the position x(t+dt) <-- x(t) + v(t)*dt + 0.5*a(t)*dt^2
							phase[jj].r.x = (phase[jj].r.x+phase[jj].v.x*dt + 0.5*accel[jj].x*dt2 + Lx)%Lx;
							phase[jj].r.y = (phase[jj].r.y+phase[jj].v.y*dt + 0.5*accel[jj].y*dt2 + Ly)%Ly;

							// calculate midpoint velocity v(t+dt/2) <-- v(t) + 0.5*a(t)*dt
							phase[jj].v.x = phase[jj].v.x+0.5*accel[jj].x*dt;
							phase[jj].v.y = phase[jj].v.y+0.5*accel[jj].y*dt;
						}

						// get the new accelerations and store them in accel
						// kk = 0, 1 , 2, ... , jj - 2, jj - 1, jj + 1, jj + 2, ... , N-2, N-1 
						for(int jj = 0 ; jj < N ; jj++){
							accel[jj] = new vector2d(); // Initializes to the zero vector
							for (int kk = 0 ; kk < jj ; kk++){
								accel[jj].add(force(phase[jj].r,phase[kk].r));
							}
							for (int kk = jj+1 ; kk < N ; kk++){
								accel[jj].add(force(phase[jj].r,phase[kk].r));
							}
							accel[jj].scale(1/phase[jj].m);
						}
						haveAccel = true; 	// only update velocities so these accelerations 
						// can be reused at the top of this method next call

						for(int jj = 0 ; jj < N ; jj++){					
							// update the velocity  v(t+dt) <-- v(t+dt/2) + 0.5*a(t+dt)*dt
							phase[jj].v.x = phase[jj].v.x+0.5*accel[jj].x*dt;
							phase[jj].v.y = phase[jj].v.y+0.5*accel[jj].y*dt;
						}
					break;
					
					case CLOSED:
						for(int jj = 0 ; jj < N ; jj++){
							// update the position x(t+dt) <-- x(t) + v(t)*dt + 0.5*a(t)*dt^2
							phase[jj].r.x = phase[jj].r.x+phase[jj].v.x*dt+0.5*accel[jj].x*dt2;
							phase[jj].r.y = phase[jj].r.y+phase[jj].v.y*dt+0.5*accel[jj].y*dt2;
							
							// calculate midpoint velocity v(t+dt/2) <-- v(t) + 0.5*a(t)*dt
							phase[jj].v.x = phase[jj].v.x+0.5*accel[jj].x*dt;
							phase[jj].v.y = phase[jj].v.y+0.5*accel[jj].y*dt;
					}
					
					// get the new accelerations and store them in accel
					// kk = 0, 1 , 2, ... , jj - 2, jj - 1, jj + 1, jj + 2, ... , N-2, N-1 
					for(int jj = 0 ; jj < N ; jj++){
						accel[jj] = new vector2d(); // Initializes to the zero vector
						for (int kk = 0 ; kk < jj ; kk++){
							accel[jj].add(force(phase[jj].r,phase[kk].r));
						}
						for (int kk = jj+1 ; kk < N ; kk++){
							accel[jj].add(force(phase[jj].r,phase[kk].r));
						}
						accel[jj].scale(1/phase[jj].m);
					}
					haveAccel = true; 	// only update velocities so these accelerations 
										// can be reused at the top of this method next call
					
					for(int jj = 0 ; jj < N ; jj++){					
						// update the velocity  v(t+dt) <-- v(t+dt/2) + 0.5*a(t+dt)*dt
						phase[jj].v.x = phase[jj].v.x+0.5*accel[jj].x*dt;
						phase[jj].v.y = phase[jj].v.y+0.5*accel[jj].y*dt;
						if(phase[jj].r.x < 0){
							phase[jj].r.x = -phase[jj].r.x;
							phase[jj].v.x = -phase[jj].v.x;
						}
						else if(phase[jj].r.x > Lx){
							phase[jj].r.x = 2*Lx - phase[jj].r.x; 
							phase[jj].v.x = -phase[jj].v.x;
						}
						if(phase[jj].r.y < 0){
							phase[jj].r.y = -phase[jj].r.y;
							phase[jj].v.y = -phase[jj].v.y;
						}
						else if(phase[jj].r.y > Ly){
							phase[jj].r.y = 2*Ly - phase[jj].r.y; 
							phase[jj].v.y = -phase[jj].v.y;
						}
					}
					break;
						
					default:
						
					break;
					
				}
				

			break;
								
			default:
				throw new IllegalStateException("Solver does not exist.");
		}
		
	}
             	
	public void initialise(Parameters params){
		
		/*
		 * All configurations assume a single mass and a single 
		 * particle size
		 * 
		 */

		double vcmX, vcmY, vsqr, KE;

		double m  = params.fget("M");
		double RR = params.fget("R"); 
		vcmX = 0;
		vcmY = 0;
		vsqr = 0;
		
		switch (initCond) {
		
			case GIVEN_KE:

				/*
				 * 
				 * TEMP QUICK FIX
				 * 
				 */
				KE = 1;
//				KE = params.iget("KE");
				
				for(int jj = 0 ; jj < N ; jj++){
					phase[jj] = new particle();
					phase[jj].m     = m;
					phase[jj].q     = 0;
					phase[jj].s     = RR;
					phase[jj].color = Color.red;
					phase[jj].r.x   = rand.nextDouble()*Lx; 
					phase[jj].r.y   = rand.nextDouble()*Ly;	
					phase[jj].v.x   = 0.5-rand.nextDouble();
					vcmX            += phase[jj].v.x;
					phase[jj].v.y   = 0.5-rand.nextDouble();
					vcmY		    += phase[jj].v.y;
				}
				vector2d vcm = new vector2d(vcmX/N, vcmY/N);
				for(int jj = 0 ; jj < N ; jj++){
					phase[jj].v.sub(vcm);
					vsqr += phase[jj].v.length2();
				}
				double rescale = Math.sqrt(2*KE/(m*vsqr));
				for(int jj = 0 ; jj < N ; jj++){
					phase[jj].v.scale(rescale);
				}
				

			break;

			case GIVEN_LATTICE:
				throw new IllegalStateException("Configuration scheme has not been coded yet.");
				
			default:
				throw new IllegalStateException("Configuration scheme does not exist.");
			
		}
			
		return;
		
		
	}	
	
	public Drawable<Gfx2D> boundaryDw() {

		return Geom2D.rectangle(new Bounds(0., Lx, 0., Ly), Color.BLACK);
	}
	
	
	public Drawable<Gfx2D> particlesDw() {
		return new Drawable<Gfx2D>() {
			public void draw(Gfx2D g) {
				for(int jj = 0 ; jj < N ; jj++){
					g.setColor(phase[jj].color);;
					g.fillCircle(phase[jj].r.x, phase[jj].r.y, phase[jj].s);
				}
			}
			public Bounds getBounds() {
				return new Bounds(0, Lx, 0, Ly);
			}			
		};
	}
	
	public double gettime(){
		
		return t;
	}
}
	
	
	

