package chris.MD.TwoD;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

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
	public static enum BC {PERIODIC, CLOSED}
	public static enum IC {MELT, VISCOUS, COPY, READ_IN}
	
	public int N, Lx, Ly;
	private particle phase[];
	private double dt, t, Tw, E[];
	private Solvers solver;
	protected BC boundCond;
	private IC initcond;
	private vector2d accel[];
	private Random rand;
	
	public InteractingSystem(Parameters params){
		
		IS_constructor(params);
	}
	
	public void IS_constructor(Parameters params){
	
		E         = new double[2];
		dt        = params.fget("dt");
		Tw        = params.fget("T");
		t         = 0;
		N         = params.iget("N");
		if(params.containsKey("L")){
			Lx = params.iget("L");
			Ly = Lx;
		}
		else{
			Lx        = params.iget("Lx");
			Ly        = params.iget("Ly"); 
		}
		accel     = new vector2d[N]; 
		phase     = new particle[N];
		rand      = new Random(params.iget("seed"));
		
		if(params.sget("ODE Solver").equals("Velocity Verlet"))
			solver    = Solvers.VERLET;	
		else if(params.sget("ODE Solver").equals("First Order Euler"))
			solver    = Solvers.EULER;	

		if(params.sget("Initial Conditions").equals("Melt"))
			initcond = IC.MELT;
		else if(params.sget("Initial Conditions").equals("Viscous"))
			initcond = IC.VISCOUS;
		else if(params.sget("Initial Conditions").equals("Copy"))
			initcond = IC.COPY;
		else if(params.sget("Initial Conditions").equals("Read In"))
			initcond = IC.READ_IN;
		
		if(params.sget("Boundary Conditions").equals("Closed"))
			boundCond = BC.CLOSED;
		else if(params.sget("Boundary Conditions").equals("Periodic"))
			boundCond = BC.PERIODIC;
		
		initialise(params);
	}
	
	
	/**
	 * Abstract method that will return central potential for
	 * two particles separated by a scalar distance r.
	 * 
	 * @param r the (scalar) distance between the two interacting particles
	 * @return the value of the potential between particles one and two
	 */
	public abstract double potential(vector2d r1, vector2d r2);
	
	/**
	 * Abstract method that will return central force for
	 * two particles separated by relative coordinates x and y.
	 * 
	 * @param r1 the vector displacement of particle 1 
	 * @param r2 the vector displacement of particle 2
	 * @return a vector2d that is the force on particle one due to the second
	 */
	public abstract vector2d force(vector2d r1, vector2d r2);
	
	public double stepWT(){
		step();
		return t;
	}
	
	public double stepWE(){
		step();
		return E[0]+E[1];
	}

	public void step(){
		
		t    += dt;
	
		switch (solver) {
		
			case EULER:

					switch (boundCond){

						case PERIODIC:
							for(int jj = 0 ; jj < N ; jj++){
								// update the position x(t+dt) <-- x(t) + v(t)*dt 
								phase[jj].r.x = (phase[jj].r.x + phase[jj].v.x*dt + Lx)%Lx;
								phase[jj].r.y = (phase[jj].r.y + phase[jj].v.y*dt + Ly)%Ly;

								// update the velocity v(t+dt) <-- v(t) + a(t)*dt
								phase[jj].v.x = phase[jj].v.x + 0.5*accel[jj].x*dt;
								phase[jj].v.y = phase[jj].v.y + 0.5*accel[jj].y*dt;
							}
							
							getAccel();
						break;

						case CLOSED:
														
							for(int jj = 0 ; jj < N ; jj++){
								// update the position x(t+dt) <-- x(t) + v(t)*dt 
								phase[jj].r.x = phase[jj].r.x + phase[jj].v.x*dt;
								phase[jj].r.y = phase[jj].r.y + phase[jj].v.y*dt;
																
								// update the velocity v(t+dt) <-- v(t) + a(t)*dt
								phase[jj].v.x = phase[jj].v.x + 0.5*accel[jj].x*dt;
								phase[jj].v.y = phase[jj].v.y + 0.5*accel[jj].y*dt;
								
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
							
							getAccel();
						break;

						default:

						break;
					}
			break;
			
			case VERLET:

				double dt2 = dt*dt;
				
				switch (boundCond){

					case PERIODIC:
						for(int jj = 0 ; jj < N ; jj++){
							// update the position x(t+dt) <-- x(t) + v(t)*dt + 0.5*a(t)*dt^2
							phase[jj].r.x = (phase[jj].r.x + phase[jj].v.x*dt + 0.5*accel[jj].x*dt2 + Lx)%Lx;
							phase[jj].r.y = (phase[jj].r.y + phase[jj].v.y*dt + 0.5*accel[jj].y*dt2 + Ly)%Ly;

							// calculate midpoint velocity v(t+dt/2) <-- v(t) + 0.5*a(t)*dt
							phase[jj].v.x = phase[jj].v.x + 0.5*accel[jj].x*dt;
							phase[jj].v.y = phase[jj].v.y + 0.5*accel[jj].y*dt;
						}

						getAccel();
						
						for(int jj = 0 ; jj < N ; jj++){					
							// update the velocity  v(t+dt) <-- v(t+dt/2) + 0.5*a(t+dt)*dt
							phase[jj].v.x = phase[jj].v.x + 0.5*accel[jj].x*dt;
							phase[jj].v.y = phase[jj].v.y + 0.5*accel[jj].y*dt;
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

						getAccel();

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
		
		return;
	}
	
	private void getAccel(){

		vector2d tmp = new vector2d();
		
		E[0] = 0;
		E[1] = 0;
		
		for(int jj = 0 ; jj < N ; jj++){
			accel[jj] = new vector2d(); // Initializes to the zero vector
		}		
		
		for(int jj = 0 ; jj < N-1 ; jj++){
			E[0] += 0.5*phase[jj].v.length2();
			for(int kk = jj+1 ; kk < N ; kk++){
				tmp = force(phase[jj].r,phase[kk].r);
				accel[jj].add(tmp);
				accel[kk].sub(tmp);
				E[1] += potential(phase[jj].r,phase[kk].r);
			}
		}
		return;
	}
             	
	public void initialise(Parameters params){

		switch(initcond){
		
		case MELT:
		
			double pcmX, pcmY, KE0; 
			int Nx, Ny;

			double m  = params.fget("M");
			double RR = params.fget("R"); 

			if(Lx == Ly){
				Nx = (int)(Math.sqrt(N));
				if(Nx*Nx != N)
					Nx++;		
				Ny = Nx;
			}
			else{
				Nx = (int)(Math.sqrt((double)(Lx)*N/Ly));
				Ny = (int)(Math.floor(N/Nx));
				Ny += N - Nx*Ny;
			}

			pcmX = 0;
			pcmY = 0;
			KE0 = 0;

			for(int jj = 0 ; jj < N ; jj++){
				phase[jj] = new particle();
				phase[jj].q     = 0;
				phase[jj].s     = RR;
				phase[jj].color = Color.red;
				phase[jj].r.x   = (jj%Nx)*Lx/Nx+(1.01)*RR;
				phase[jj].r.y   = ((int)(jj/Nx))*Ly/Ny+(1.01)*RR;
				phase[jj].v.x   = 0.5-rand.nextDouble();
				pcmX            += m*phase[jj].v.x;
				phase[jj].v.y   = 0.5-rand.nextDouble();
				pcmY		    += m*phase[jj].v.y;
			}
			vector2d pcm = new vector2d(pcmX/N, pcmY/N);
			for(int jj = 0 ; jj < N ; jj++){
				phase[jj].v.sub(pcm);;
				KE0 += 0.5*phase[jj].v.length2();
			}
			double rescale = Math.sqrt((N*Tw)/(KE0));
			for(int jj = 0 ; jj < N ; jj++){
				phase[jj].v.scale(rescale);
			}
		
		break;
		
		case VISCOUS:
			
		case COPY:
			
		case READ_IN:
			
			throw new IllegalStateException("Not coded yet.");

		default:

			throw new IllegalStateException("Initialization scheme does not exist.");
		
		}
		
		getAccel();		
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
	
	public double getEsys(){
		
		return E[0] + E[1];
	}

	
	public void getEsys(double[] E){
		
		E[0] = this.E[0];
		E[1] = this.E[1];
		return;
	}
	
	public void changeT(double T){
		
		double KE0 = 0;
		for(int jj = 0 ; jj < N ; jj++){
			KE0 += 0.5*phase[jj].v.length2();
		}
		double rescale = Math.sqrt((N*T)/(KE0));
		for(int jj = 0 ; jj < N ; jj++){
			phase[jj].v.scale(rescale);
		}
		return;
	}
	
	public void setStepSize(double dt){
		
		this.dt = dt;
		return;
	}
	
	public double getTemp(){
		double T = 0;
		
		for(int jj = 0 ; jj < N ; jj++)
			T += 0.5*phase[jj].v.length2();

		return T/N;
	}
	
	public void printPhase(String fout){
		
		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.println("ID \t x \t y \t v_x \t v_y");
			pw.println();
			for (int jj = 0 ; jj < N ; jj++){
				pw.print(jj);
				pw.print("\t");
				pw.print(phase[jj].r.x);
				pw.print("\t");
				pw.print(phase[jj].r.y);
				pw.print("\t");
				pw.print(phase[jj].v.x);
				pw.print("\t");
				pw.print(phase[jj].v.y);
				pw.println();
			}
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
	return;
	}
}
	
/*
 * 	periodic boundary condition debugging data
 */
	
//phase[0] = new particle();
//phase[0].m     = m;
//phase[0].q     = 0;
//phase[0].s     = RR;
//phase[0].color = Color.red;
//phase[0].r.x   = 45.30531741317776;
//phase[0].r.y   = 13.696366621853066;
//phase[0].v.x   = -0.016131747215712464;
//phase[0].v.y   = 0.03230476536155055;
//
//phase[1] = new particle();
//phase[1].m     = m;
//phase[1].q     = 0;
//phase[1].s     = RR;
//phase[1].color = Color.red;
//phase[1].r.x   = 29.045508522457197;
//phase[1].r.y   = 42.13277977321583;
//phase[1].v.x   = -0.028468880134714613;
//phase[1].v.y   = -0.0024605691974465614;
//
//phase[2] = new particle();
//phase[2].m     = m;
//phase[2].q     = 0;
//phase[2].s     = RR;
//phase[2].color = Color.red;
//phase[2].r.x   = 14.576730072738897;
//phase[2].r.y   = 0.45182210819004354;
//phase[2].v.x   = 0.037416089146865096;
//phase[2].v.y   = -0.037073919915994055;
//
//phase[3] = new particle();
//phase[3].m     = m;
//phase[3].q     = 0;
//phase[3].s     = RR;
//phase[3].color = Color.red;
//phase[3].r.x   = 16.05978527112019;
//phase[3].r.y   = 48.850946642337206;
//phase[3].v.x   = -0.07665533054105472;
//phase[3].v.y   = 0.09827745829828256;

