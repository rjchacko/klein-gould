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
import chris.util.ReadInUtil;
import chris.util.vector2d;

public abstract class InteractingSystem {

	public static enum Solvers {EULER, VERLET};
	public static enum BC {PERIODIC, CLOSED}
	public static enum IC {MELT, VISCOUS, COPY, READ_IN, DEBUG}
	
	public int N;
	public double Lx, Ly;
	private particle phase[];
	private double dt, t, Tw, E[];
	private Solvers solver;
	protected BC boundCond;
	private IC initcond;
	private vector2d accel[], rij[][]; // rij[i][j] is r_i - r_j
	private Random rand;
	
	public InteractingSystem(Parameters params){
		
		IS_constructor(params);
	}
	
	public void IS_constructor(Parameters params){
	
		E  = new double[2];
		dt = params.fget("dt");
		Tw = params.fget("T");
		t  = 0;
		N  = params.iget("N");
		if(params.containsKey("L")){
			Lx = params.fget("L");
			Ly = Lx;
		}
		else{
			Lx = params.fget("Lx");
			Ly = params.fget("Ly"); 
		}
		rand = new Random(params.iget("seed"));
		
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
		else if(params.sget("Initial Conditions").equals("Debug"))
			initcond = IC.DEBUG;
		
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
	public abstract double potential(vector2d deltaR);
	
	/**
	 * Abstract method that will return central force for
	 * two particles separated by relative coordinates x and y.
	 * 
	 * @param r1 the vector displacement of particle 1 
	 * @param r2 the vector displacement of particle 2
	 * @return a vector2d that is the force on particle one due to the second
	 */
	public abstract vector2d force(vector2d deltaR);
	
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
				// get distance between jj and kk
				rij[jj][kk] = dr(phase[jj].r,phase[kk].r);
				// simple vector maths for the opposite relative distance
				rij[kk][jj] = rij[jj][kk].negate();
				// get the force
				tmp = force(rij[jj][kk]);
				accel[jj].add(tmp);
				// use Newton III
				accel[kk].sub(tmp);
				E[1] += potential(rij[jj][kk]);
			}
		}
		return;
	}
             	
	public void initialise(Parameters params){

		switch(initcond){
		
		case MELT:
		
			double pcmX, pcmY, KE0; 
			int Nx, Ny;

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

			pcmX  = 0;
			pcmY  = 0;
			KE0   = 0;
			phase = new particle[N];
			accel = new vector2d[N]; 
			
			for(int jj = 0 ; jj < N ; jj++){
				phase[jj] = new particle();
				phase[jj].q     = 0;
				phase[jj].s     = RR;
				phase[jj].color = Color.red;
				phase[jj].r.x   = (jj%Nx)*Lx/Nx+(1.01)*RR;
				phase[jj].r.y   = ((int)(jj/Nx))*Ly/Ny+(1.01)*RR;
				phase[jj].v.x   = 0.5-rand.nextDouble();
				pcmX            += phase[jj].v.x;
				phase[jj].v.y   = 0.5-rand.nextDouble();
				pcmY		    += phase[jj].v.y;
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
		
		case READ_IN:
			
			ReadInUtil ri = new ReadInUtil("/Users/cserino/Documents/BU/Research/MD/IC.txt");
			double ic[][] = ri.getDataBeforeString(new int[]{2,3,4,5}, 2,"-------- Parameters --------");
			N = ic[0].length;
			phase = new particle[N];
			accel = new vector2d[N]; 
			params.set("N", N);
			
			for(int jj = 0 ; jj < N ; jj++){
				phase[jj] = new particle();
				phase[jj].q     = 0;					// FIX THIS
				phase[jj].s     = params.fget("R"); 	// AND THIS
				phase[jj].color = Color.red;			// AND THIS
				phase[jj].r.x   = ic[0][jj];
				phase[jj].r.y   = ic[1][jj];
				phase[jj].v.x   = ic[2][jj];
				phase[jj].v.y   = ic[3][jj];
			}

		break;
		
		case VISCOUS:
			
		case DEBUG:

			params.set("L",15);
			params.set("N",225);
			params.set("M",1);
			params.set("T",0);
			
			N  = params.iget("N");
			Lx = params.fget("L");
			Ly = Lx;
			
			phase = new particle[N];
			accel = new vector2d[N]; 

			for(int kk = 0 ; kk < Ly ; kk++){
				for(int jj = 0 ; jj < Lx ; jj++){
					phase[jj+kk*(int)Ly] = new particle();
					phase[jj+kk*(int)Ly].q     = 0;
					phase[jj+kk*(int)Ly].s     = 0.01;
					phase[jj+kk*(int)Ly].color = Color.red;
					phase[jj+kk*(int)Ly].v.x   = 0;
					phase[jj+kk*(int)Ly].v.y   = 0;
					phase[jj+kk*(int)Ly].r.y   = jj/(double)(Lx);
					if(jj%2 == 0)
						phase[jj+kk*(int)Ly].r.x = kk/(double)(Ly);
					else
						phase[jj+kk*(int)Ly].r.x = (2*kk+1)/(2.*Ly);
				}
			}
//			printPhase("/Users/cserino/Desktop/debugConfig.txt");
			
			params.set("L",1);
			Lx = params.iget("L");
			Ly = Lx;

		break;
			
		case COPY:
			
			throw new IllegalStateException("Not coded yet.");
						
		default:

			throw new IllegalStateException("Initialization scheme does not exist.");
		}
		
		rij = new vector2d[N][N];
		for(int jj = 0 ; jj  < N ; jj++){ // the off diagonal terms are computed in getAccel()
			rij[jj][jj] = new vector2d();
			
		}
		getAccel();		
		
		if(initcond.equals(IC.READ_IN)){
			params.set("T",E[0]/N);
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
	
	public void printPhase(String fout, Parameters params){
		
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
			pw.println("-------- Parameters --------");
			pw.println(params.toString());
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
	return;
	}
	
	public double structureFactor(int nx, int ny, double Lx, double Ly){
	
		double dot;		
		vector2d k = new vector2d(Math.PI*nx/Lx, Math.PI*ny/Ly);
		double rp, ip, sf;  
		rp = ip = sf = 0;
		
		for(int jj = 0 ; jj < N ; jj++){
			for(int kk = 0 ; kk < N ; kk++){
				dot = vector2d.dot(k, rij[jj][kk]);
				rp += Math.cos(dot);
				ip += Math.sin(dot);
			}
			sf += (rp*rp+ip*ip);
		}
		return sf/N;
	}
	
	public vector2d dr(vector2d r1, vector2d r2){
		
		switch (boundCond){
		
			case PERIODIC: 
				
				double dx = r1.x-r2.x;
				double dy = r1.y-r2.y;
				if(dx*dx > 0.25*Lx*Lx)
					dx = -Math.signum(dx)*(Lx-Math.abs(dx));
				if(dy*dy > 0.25*Ly*Ly)
					dy = -Math.signum(dy)*(Ly-Math.abs(dy));
				return new vector2d(dx, dy);
		
			case CLOSED:
			
				return new vector2d(r1.x-r2.x, r1.y-r2.y);

			default:
				throw new IllegalStateException("Boundary condition does not exist.");
		}
	}
}
