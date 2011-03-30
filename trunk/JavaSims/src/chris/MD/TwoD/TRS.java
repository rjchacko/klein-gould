package chris.MD.TwoD;

import scikit.jobs.params.Parameters;
import chris.util.vector2d;

public class TRS extends InteractingSystem {
	
	/*
	 * 	From:
	 * 		  M. F. Laguna and E. A. Jagla </i>J. Stat. Mech.</i> P09002, (2009)
	 * 
	 */
	
	/*
	 * 	NOTE: V0 takes the </i> square </i> of the inter-particle
	 *        distance as the input whereas V1, V2, V3, F0, F1, F2, 
	 *        and F3 take the inter-particle distance as the input.
	 * 
	 */
	
	/*
	 * The set of parameters P = {A0, A2, A3, c, d2, s2, d3, s3}.
	 * P1 and P2 which drive a triangularÐrhombohedral transition (TÐR) and a 
	 * triangularÐsquare transition (TÐS) respectively:
	 * 
	 * P1	=	{A0, 0.003, 0.01, 1.722, 0.98, 0.04, 1.74, 0.2}	
	 * 			with variable A0 and Ac0 = 0.067.
	 * 
	 * P2	=	{0.024, A2, 0.01, 1.730, 0.98, 0.1, 1.74, 0.2} 
	 * 			with variable A2 and Ac2 = 0.022.
	 * 
	 * Note:  Ac1/2 is a "transition value"
	 * 
	 */
	
	double A0, A2;
	
//	// Studying the T-S transition
//	static final double A3   = 0.01;
//	static final double c    = 1.722;
//	static final double d2   = 0.98;
//	static final double s2   = 0.04;
//	static final double d3   = 1.74;
//	static final double s3   = 0.2;
//	static final double dlt2 = d2-s2;
//	static final double dlt3 = d3-s3;
//	static final double sgm2 = d2+s2;
//	static final double sgm3 = d3+s3;
//	static final double cm14 = 1./((c-1.)*(c-1.)*(c-1.)*(c-1.));
//	static final double s24  = 1./(s2*s2*s2*s2);
//	static final double s34  = 1./(s3*s3*s3*s3);

	// Studying the T-S transition
	static final double A3   = 0.01;
	static final double c    = 1.730;
	static final double d2   = 0.98;
	static final double s2   = 0.1;
	static final double d3   = 1.74;
	static final double s3   = 0.2;
	static final double dlt2 = d2-s2;
	static final double dlt3 = d3-s3;
	static final double sgm2 = d2+s2;
	static final double sgm3 = d3+s3;
	static final double cm14 = 1./((c-1.)*(c-1.)*(c-1.)*(c-1.));
	static final double s24  = 1./(s2*s2*s2*s2);
	static final double s34  = 1./(s3*s3*s3*s3);

	public TRS(Parameters params) {
		
		super(params);
		
		// for now, pick T-S transition 
		A0 = 0.024;
		//A2 = 0.03; //(A2 > Ac2) stable phase is triangle
			       // quench to A2 = 0.013 < Ac2 (square is stable)
		A2 = 0.013;
		return;
	}

	public vector2d force(vector2d deltaR) {

		double dr2 = deltaR.length2();
		if(dr2 >= sgm3*sgm3)  // if true, all potentials are zero.
			return new vector2d(); // the zero vector
		
		double r = Math.sqrt(deltaR.length2());		
		return vector2d.scale(vector2d.scale(deltaR,-1./r),F0(r)+F1(r)+F2(r)+F3(r));
	}

	public double potential(vector2d deltaR) {

		double dr2 = deltaR.length2();
		if(dr2 >= sgm3*sgm3)  // if true, all potentials are zero.
			return 0;
		
		double r = Math.sqrt(dr2);
		return V0(dr2) + V1(r) + V2(r) + V3(r);
	}
	
	private double V0(double dr2){
		
		if(dr2 >= 1)
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
		
		if((dlt2 <= r) || (r >= sgm2))
			return 0;
		
		return -A2*(r-sgm2)*(r-sgm2)*(r-dlt2)*(r-dlt2)*s24;
	}
	
	private double V3(double r){
		
		if((dlt3 <= r) || (r >= sgm3))
			return 0;
		
		return A3*(r-sgm3)*(r-sgm3)*(r-dlt3)*(r-dlt3)*s34;
	}
	
	/*
	 * 	The F_i are magnitudes and, as such, are Grad V (</i> not - Grad V)
	 * 
	 */
	
	public double F0(double r) {
		
		if(r*r >= 1)
			return 0;
		
		return (-12*V0(r*r)-A0)/r;
	}
	
	public double F1(double r) {
		
		if( r >= c)
			return 0;
		
		return 2*(V1(r)+1)*(1./(r-1) + 1./(r-2*c+1));
	}
	
	public double F2(double r) {
		
		return 2*V2(r)*(1./(r-sgm2) + 1./(r-dlt2));
	}
	
	public double F3(double r) {
		
		return 2*V3(r)*(1./(r-sgm3) + 1./(r-dlt3));
	}
	
	public void quenchA0(double A0){

		this.A0 = A0;
		return;
	}
	
	public void quenchA2(double A2){

		this.A2 = A2;
		return;
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
