package chris.MD.TwoD;

import scikit.jobs.params.Parameters;
import chris.util.vector2d;

public class TRS extends InteractingSystem {
	
	/*
	 * 	NOTE: V0 and F0 take the </i> square </i> of the inter-particle
	 *        distance as the input whereas V1, V2, V3, F1, F2, and F3
	 *        take the inter-particle distance as the input.
	 * 
	 */

	static final double A0   = 1;
	static final double A2   = 1;
	static final double A3   = 1;
	static final double c    = 1;
	static final double d2   = 1;
	static final double s2   = 1;
	static final double d3   = 1;
	static final double s3   = 1;
	static final double dlt2 = d2-s2;
	static final double dlt3 = d3-s3;
	static final double sgm2 = d2+s2;
	static final double sgm3 = d3+s3;
	static final double cm14 = 1./((c-1.)*(c-1.)*(c-1.)*(c-1.));
	static final double s24  = 1./(s2*s2*s2*s2);
	static final double s34  = 1./(s3*s3*s3*s3);

	public TRS(Parameters params) {
		
		super(params);
	}

	public vector2d force(vector2d deltaR) {

		double dr2 = deltaR.length2();
		if(dr2 >= sgm3*sgm3)  // if true, all potentials are zero.
			return new vector2d(); // the zero vector
		
		double r = Math.sqrt(deltaR.length2());
		deltaR.scale(-1./r); // make a unit vector pointing in opposite direction
		
		return vector2d.scale(deltaR,F0(r*r)+F1(r)+F2(r)+F3(r));
	}

	public double potential(vector2d deltaR) {

		double dr2 = deltaR.length2();
		if(dr2 >= sgm3*sgm3)  // if true, all potentials are zero.
			return 0;
		
		double r = Math.sqrt(dr2);
		return V0(dr2) + V1(r) + V2(r) + V3(r);
	}
	
	private double V0(double dr2){
		
		if(dr2 < 1)
			return 0;
		
		double r6  = dr2*dr2*dr2;
		double r12 = r6*r6; 
		return A0*(1./r12 - 2./r6 + 1);
	}
	
	private double V1(double r){

		if( r >= c)
			return 0;
		
		return ((r-1)*(r-1)*(r+1-2*c)*(r+1-2*c)*cm14) - 1;
	}
	
	private double V2(double r){
		
		if((dlt2 >= r) || (r <= sgm2))
			return 0;
		
		return -A2*(r-sgm2)*(r-sgm2)*(r-dlt2)*(r-dlt2)*s24;
	}
	
	private double V3(double r){
		
		if((dlt3 >= r) || (r <= sgm3))
			return 0;
		
		return A3*(r-sgm3)*(r-sgm3)*(r-dlt3)*(r-dlt3)*s34;
	}
	
	/*
	 * 	The F_i are magnitudes and, as such, are Grad V (</i> not - Grad V)
	 * 
	 */
	
	public double F0(double r2) {
		
		return 0;
	}
	
	public double F1(double r) {
		
		return 0;
	}
	
	public double F2(double r) {
		
		return 0;
	}
	
	public double F3(double r) {
		
		return 0;
	}
	
//	 /**
//     * Returns the vector sum: sum_{i=1...4} vector_i.
//     * 
//     * @param v1 the first vector
//     * @param v2 the second vector
//     * @param v3 the first vector
//     * @param v4 the second vector
//     * @return the vector sum_{i=1...4} v_i.
//     */
//    public static final vector2d add4(vector2d v1, vector2d v2, vector2d v3, vector2d v4)
//    {
//    	
//    	return new vector2d(v1.x+v2.x+v3.x+v4.x , v1.y+v2.y+v3.y+v4.y);
//    }


}
