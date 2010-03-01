package chris.MD.TwoD;

import java.awt.Color;

import chris.util.vector2d;

public class particle {
	
	public particle(){
		
		r     = new vector2d();
		v     = new vector2d();
		q     = 0;
		s     = 0;
		color = Color.WHITE;
	}
	
	
	/**
	 * the position of the particle
	 */
	public vector2d r;
	
	/**
	 * the velocity of the particle
	 */
	public vector2d v;
	
	/**
	 * the electric charge of the particle
	 */
	public double q;
	
	/**
	 * the radius of the particle
	 */
	public double s;

	/**
	 * the net force on the particle
	 */
//	public vector2d f;
//	
//	/**
//	 * the potential energy of the particle
//	 */
//	public double U;
//	
//	/**
//	 * the kinetic energy of the particle
//	 */
//	public double T;
	
	/**
	 * for graphics
	 */
	public Color color;

	
}
