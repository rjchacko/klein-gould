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
	
	public static enum Shape {Circle, Square, Diamond, Ellipse, All};

	public LatticeNeighbors(int Nx, int Ny, double r_lo, double r_hi, Type type, Shape shape) {
		this.Nx = Nx;
		this.Ny = Ny;
		this.r_lo = r_lo;
		this.r_hi = r_hi;
		this.type = type; 
		this.shape = shape;
		int d = (int)r_hi*2 + 1;
		if(!(shape.equals(Shape.All))){
			list = new int[d*d];
		}
		else{
			list = null;
		}
	}
	
	public LatticeNeighbors(int Nx, int Ny, double eps, double r_lo, double r_hi, Type type, Shape shape) {
		this.Nx = Nx;
		this.Ny = Ny;
		this.eps = eps;
		this.r_lo = r_lo;
		this.r_hi = r_hi;
		this.type = type; 
		this.shape = shape;
		int d = (int)r_hi*2 + 1;
		list = new int[d*d];
		eps2 = eps*eps;
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
	 * Returns a neighbor list array.
	 * For each lattice index i, neighbor[i] is an array containing all of i's neighbor
	 * indices. Neighbor indices are defined to be those whose distance d from i is
	 * in the interval [r_lo, r_hi).
	 * 
	 * This is a faster method than above but requires at least one call
	 * to the method above prior to calling this one. It simply translates
	 * all the neighbor sites by the vector which points from the site with
	 * known neighbors to the site with neighbors to be calculated.
	 * 
	 * For case = BORDERED the site i0 *must* be positioned on the lattice such
	 * that its set of neighbors is the same set that would be returned had the 
	 * method been called with PERIODIC boundary conditions. Namely, i0 must be
	 * sufficiently centered on the lattice such that there are no edge effects.
	 * 
	 * For case = PERIODIC the site i0 *must* be taken as i0 = 0 or else you will
	 * get an array index error resulting from JAVA's % operator which is not the
	 * true modulo operator.
	 * 
	 */
	public int[] get(int i, int i0, int[] nbs0) {
		
		int Nsites = nbs0.length;
		int[] temp = new int[Nsites];
		
		int dx = i%Nx - i0%Nx;
		int dy = (int)(i/Ny) - (int)(i0/Ny);
		
		if(shape == LatticeNeighbors.Shape.All){
			return nbs0;
		}
		
		switch (type) {
		
		case BORDERED:
			
			int counter = 0;
			
			for (int jj = 0 ; jj < Nsites ; jj++){
				int xn = nbs0[jj]%Nx;
				int yn = (int)(nbs0[jj]/Ny);
				
				int xx = xn + dx;
				int yy = yn + dy;
					
				if (xx >= 0 && xx < Nx && yy >= 0 && yy < Ny) temp[counter++] = xx + yy*Ny;
	
			}
			
			int[] ret = new int[counter];
			
			for (int jj = 0 ; jj < counter ; jj++){
				ret[jj] = temp[jj];
			}
			
			return ret;
			
		case PERIODIC:
					
			int[] retP = new int[Nsites];

			for (int jj = 0 ; jj < Nsites ; jj++){
				int xn = nbs0[jj]%Nx;
				int yn = (int)(nbs0[jj]/Ny);
				retP[jj] = (xn+dx)%Nx + ((yn + dy)%Ny)*Ny;
			}

			return retP;
			
		default:
			
			return null;
			
		}
		
	}
	
	/**
	 * Returns the Jth element of the neighbor list array
	 * using the method of translation.
	 * 
	 **/
	
	
	public int getJ(int i, int i0, int[] nbs0, int J){
		
		if(shape == LatticeNeighbors.Shape.All) return (J == i) ? -1 : J;
		
		int xn, yn;
		int dx = i%Nx - i0%Nx;
		int dy = (int)(i/Ny) - (int)(i0/Ny);
		
		switch (type) {
	
		case BORDERED:
			
			xn = nbs0[J]%Nx;
			yn = (int)(nbs0[J]/Ny);
			int xx = xn + dx;
			int yy = yn + dy;
				
			return (xx >= 0 && xx < Nx && yy >= 0 && yy < Ny) ?  xx + yy*Ny : -1;
	
		case PERIODIC:
			
			xn = nbs0[J]%Nx;
			yn = (int)(nbs0[J]/Ny);
			
			return (xn+dx)%Nx + ((yn + dy)%Ny)*Ny;	
			
		default:


		}
		
		return -77; // only get here via default
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
	double eps, eps2; //eccentricity and its square (respectively) of the ellipse

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
							if((jy*Nx + jx)==i) continue;
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
							if((jy*Nx + jx)==i) continue;
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
			

		case Ellipse:
		
				switch (type) {
				case BORDERED:
					for (int jy = Math.max(0, iy-_r); jy <= Math.min(Ny-1, iy+_r); jy++) {
						for (int jx = Math.max(0, ix-_r); jx <= Math.min(Nx-1, ix+_r); jx++) {
							if((jy*Nx + jx)==i) continue;
							double d2 = sqr(jx-ix) + sqr(jy-iy)/(1-eps2);
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
							if((jy*Nx + jx)==i) continue;
							double d2 = sqr(jx-ix) + sqr(jy-iy)/(1-eps2);
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
						if((jy*Nx + jx)==i) continue;
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
						if((jy*Nx + jx)==i) continue;
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
						if((jy*Nx + jx)==i) continue;
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
						if((jy*Nx + jx)==i) continue;
						int _jy = (jy + Ny) % Ny;
						int _jx = (jx + Nx) % Nx;
												
						if ( rmin <= (Math.abs(jx-ix) + Math.abs(jy-iy)) && _r >= (Math.abs(jx-ix) + Math.abs(jy-iy)) ) list[num++] = _jy*Nx + _jx;			
					}
				}
				
				break;
			}
			
			break;
			
		case All:

			for (int jj = 0 ; jj < Nx*Ny ; jj++){
				if(jj==i) continue;
				list[num++] = jj;
			}
			
			break;
			
		}
	}
}
