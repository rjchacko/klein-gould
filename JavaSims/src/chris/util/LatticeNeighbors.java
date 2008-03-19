package chris.util;


import static java.lang.Math.ceil;
import static java.lang.Math.min;
import static scikit.numerics.Math2.sqr;


/**
 * Utility class to calculate lattice neighbor lists.
 */

public class LatticeNeighbors {
	/**
	 * PERIODIC -- Indices on the edges have fewer neighbors.
	 * BORDERED -- Indices on the edges wrap around.
	 */
	public static enum Type {BORDERED, PERIODIC};
	
	public static enum Shape {Circle, Square, Diamond};

	public LatticeNeighbors(int Nx, int Ny, double r_lo, double r_hi, Type type, Shape shape) {
		this.Nx = Nx;
		this.Ny = Ny;
		this.r_lo = r_lo;
		this.r_hi = r_hi;
		this.type = type; 
		this.shape = shape;
		int d = (int)r_hi*2 + 1;
		list = new int[d*d];
	}
	
	/**
	 * Returns a neighbor list array.
	 * For each lattice index i, neighbor[i] is an array containing all of i's neighbor
	 * indices. Neighbor indices are defined to be those whose distance d from i is
	 * in the interval [r_lo, r_hi).
	 */
	public int[] get(int i) {
		calculateNeighborList(i);
		int[] ret = new int[num];
		System.arraycopy(list, 0, ret, 0, num);
		return ret;
	}

	/**
	 * Returns a static neighbor list array, terminated by index -1.
	 * For each lattice index i, neighbor[i] is an array containing all of i's neighbor
	 * indices.
	 * The returned array will be reused for each call to unsafeGet()
	 */
	public int[] unsafeGet(int i) {
		calculateNeighborList(i);
		list[num] = -1;
		return list;
	}

	// -------------------------------------------------------------------------------------
	// Protected methods
	// 

	int Nx, Ny;    // lattice size
	double r_lo, r_hi; // neighbor list range
	Type type;
	Shape shape;

	// These variables are set be calculateNeighborList()
	int[] list;    // most recently computed neighbor list
	int num;       // number of elements

	void calculateNeighborList(int i) {
		int _r = (int)ceil(r_hi);
		int ix = i % Ny;
		int iy = i / Nx;
		num = 0;

		
		switch (shape) {
			
		case Circle:
		
				switch (type) {
				case BORDERED:
					for (int jy = Math.max(0, iy-_r); jy <= Math.min(Ny-1, iy+_r); jy++) {
						for (int jx = Math.max(0, ix-_r); jx <= Math.min(Nx-1, ix+_r); jx++) {
							double d2 = sqr(jx-ix) + sqr(jy-iy);
							if (sqr(r_lo) <= d2 && d2 <= sqr(r_hi))
								list[num++] = jy*Nx + jx;
						}
					}
					break;

				case PERIODIC:
					// take advantage of rounding down when Ny or Nx is even.  this way
					// we avoid double counting (-Ny/2 == +Ny/2).
					// for odd Ny or Nx, this ambiguity does not exist, but the subtraction
					// by one will not have an effect.
					int ylo = iy - min(_r, (Ny-1)/2);
					int yhi = iy + min(_r, Ny/2);
					int xlo = ix - min(_r, (Nx-1)/2);
					int xhi = ix + min(_r, Nx/2);
			
					for (int jy = ylo; jy <= yhi; jy++) {
						for (int jx = xlo; jx <= xhi; jx++) {
							double d2 = sqr(jx-ix) + sqr(jy-iy);
							if (sqr(r_lo) <= d2 && d2 <= sqr(r_hi)) {
								int _jy = (jy + Ny) % Ny;
								int _jx = (jx + Nx) % Nx;
								list[num++] = _jy*Nx + _jx;
							}
						}
					}
					break;
				}
				break;
			
		case Square:
			
			switch (type) {
			case BORDERED:
				for (int jx = Math.max(0,ix-_r); jx <= Math.min(Nx-1,ix+_r); jx++){
					for (int jy = Math.max(0,iy-_r); jy <= Math.min(Ny-1,iy+_r); jy++){
						
						if( (int)(r_lo) <= Math.abs(jx-ix) || (int)(r_lo) <= Math.abs(jy-iy) ) list[num++] = jy*Nx + jx;
					}
				}

				break;
				
			case PERIODIC:
				
				int ylo = iy - min(_r, (Ny-1)/2);
				int yhi = iy + min(_r, Ny/2);
				int xlo = ix - min(_r, (Nx-1)/2);
				int xhi = ix + min(_r, Nx/2);
				
				for (int jy = ylo; jy <= yhi; jy++) {
					for (int jx = xlo; jx <= xhi; jx++) {
						int _jy = (jy + Ny) % Ny;
						int _jx = (jx + Nx) % Nx;
						
						if( (int)(r_lo) <= Math.abs(_jx-ix) || (int)(r_lo) <= Math.abs(_jy-iy) ) list[num++] = _jy*Nx + _jx;			
					}
				}
				
				break;
				
			}
			
			break;
			
		case Diamond:
			
			int rmin = (int)(r_lo);
			
			switch (type) {
			case BORDERED:
				for (int jx = Math.max(0,ix-_r); jx <= Math.min(Nx-1,ix+_r); jx++){
					for (int jy = Math.max(0,iy-_r); jy <= Math.min(Ny-1,iy+_r); jy++){
				
						if ( rmin <= (Math.abs(jx-ix) + Math.abs(jy-iy)) && _r >= (Math.abs(jx-ix) + Math.abs(jy-iy)) ) list[num++] = jy*Nx + jx;
						
					}
				}
				
				break;
				
			case PERIODIC:
				
				
				int ylo = iy - min(_r, (Ny-1)/2);
				int yhi = iy + min(_r, Ny/2);
				int xlo = ix - min(_r, (Nx-1)/2);
				int xhi = ix + min(_r, Nx/2);
				
				for (int jy = ylo; jy <= yhi; jy++) {
					for (int jx = xlo; jx <= xhi; jx++) {
						int _jy = (jy + Ny) % Ny;
						int _jx = (jx + Nx) % Nx;
												
						if ( rmin <= (Math.abs(jx-ix) + Math.abs(jy-iy)) && _r >= (Math.abs(jx-ix) + Math.abs(jy-iy)) ) list[num++] = _jy*Nx + _jx;			
					}
				}
				
				break;
			}
			
			break;
			
		}
	}
}
