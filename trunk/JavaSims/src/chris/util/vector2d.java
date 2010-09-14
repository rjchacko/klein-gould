package chris.util;


/**
 * A 2-element vector that is represented by double-precision floating point 
 * <code>(x<sub>1</sub>,x<sub>2</sub>)</code> coordinates.  The vector is 
 * assumed to live in <code>E<sup>2</sup></code>. If this value represents a 
 * normal, then it should be normalized.
 * 
 * Based on K. Barros's Tuple3d from SciKit
 *
 */
public class vector2d {

	/**
     * The x_1 coordinate.
     */
    public	double	x;

    /**
     * The x_2 coordinate.
     */
    public	double	y;


    /*	====================================================================================
     * 		CONSTRUCTORS
     * 	====================================================================================
     */    
    
    /**
     * Constructs and initializes a vector2d from the specified xy coordinates.
     * @param x the x coordinate
     * @param y the y coordinate
     */
    public vector2d(double x, double y)
    {
	this.x = x;
	this.y = y;
    }

    /**
     * Constructs and initializes a vector2d from the array of length 2.
     * @param v the array of length 2 containing xy in order
     */
    public vector2d(double[] v)
    {
	this.x = v[0];
	this.y = v[1];
    }

    /**
     * Constructs and initializes a vector2d from the specified vector2d.
     * @param v the vector2d containing the initialization x y data
     */
    public vector2d(vector2d v)
    {
	this.x = v.x;
	this.y = v.y;
    }
    
    /**
     * Constructs and initializes a vector2d to (0,0).
     */
    public vector2d()
    {
	this.x = (double) 0.0;
	this.y = (double) 0.0;
    }

    
    /*	====================================================================================
     * 		STATIC METHODS
     * 	====================================================================================
     */
    
    /**
     * Returns the vector sum of vector1 and vector2. 
     * 
     * @param v1 the first vector
     * @param v2 the second vector
     * @return the vector v1 + v2
     */
    public static final vector2d add(vector2d v1, vector2d v2)
    {
    	
    	return new vector2d(v1.x+v2.x,v1.y+v2.y);
    }
    
    /**
     * Returns the vector difference between vector1 and vector2. 
     * 
     * @param v1 the first vector
     * @param v2 the second vector
     * @return the vector v1 - v2
     */
    public static final vector2d sub(vector2d v1, vector2d v2)
    {
    	
    	return new vector2d(v1.x-v2.x,v1.y-v2.y);
    }
    
    
    /**
     * Sets the value of this vector2d to the negation of vector2d v.
     * @param v the source tuple
     */
    public static final vector2d negate(vector2d v)
    {
    	return(new vector2d(-v.x,-v.y));
    }

    
    /**
     * Returns the vector after scalar multiplication. 
     * 
     * @param v the first vector
     * @param s the scalar
     * @return the vector s * v1
     */
    
    public static final vector2d scale(vector2d v, double s)
    {
    	
    	return new vector2d(s*v.x,s*v.y);
    }
    
    
    /**
     * Returns a string that contains the values of this Tuple3d.
     * The form is (x,y).\
     * @param v the vector
     * @return the String representation
     */  
    public static final String toString(vector2d v) {
        return "(" + v.x + ", " + v.y + ")";
    }



   /**
     * Returns true if the vectors are equal, element by element.
     * Returns false otherwise.
     * @param v1 the first vector
     * @param v2 the second vector
     * @return  true or false
     */  
    public static final boolean equals(vector2d v1, vector2d v2)
    {
      try {
        return(v1.x == v2.x && v1.y == v2.y);
      }
      catch (NullPointerException e2) {return false;}
    }
    
	/**
	 * Returns the dot product (using a Euclidean metric) of vector v1 and v2.
	 * @param v1 the first vector
	 * @param v2 the second vector
	 * @return the dot product of v1 and v2
	 */
	public static final double dot(vector2d v1, vector2d v2){
		
		return v1.x*v2.x + v1.y*v2.y;
	}
    /*	====================================================================================
     * 		PUBLIC SET / GET METHODS 
     * 	====================================================================================
     */
    
    
    /**
     * Sets the value of this vector2d to the specified xy coordinates.
     * @param x the x coordinate
     * @param y the y coordinate
     */
    public final void set(double x, double y)
    {
	this.x = x;
	this.y = y;
    }

    /**
     * Sets the value of this vector2d to the value of the xy coordinates
     * located in the array of length 2.
     * @param v the array of length 2 containing xy in order
     */
    public final void set(double[] v)
    {
	this.x = v[0];
	this.y = v[1];
    }

    /**
     * Sets the value of this tuple to the value of tuple t1.
     * @param v the tuple to be copied
     */
    public final void set(vector2d v)
    {
	this.x = v.x;
	this.y = v.y;
    }

   /**
     * Copies the x,y coordinates of this tuple into the array t
     * of length 2.
     * @param v  the target array 
     */
    public final void get(double[] v)
    {
        v[0] = this.x;
        v[1] = this.y;
    }


   /**
     * Copies the x,y, coordinates of this tuple into the tuple t.
     * @param v the Tuple3d object into which the values of this object are copied
     */
    public final void get(vector2d v)
    {
        v.x = this.x;
        v.y = this.y;
    }


    /**
	 * Get the <i>x_1</i> coordinate.
	 * 
	 * @return the <i>x_1</i> coordinate.
	 */
	public final double getX() {
		return x;
	}


	/**
	 * Set the <i>x_1</i> coordinate.
	 * 
	 * @param x  value to <i>x_1</i> coordinate.
	 */
	public final void setX(double x) {
		this.x = x;
	}

    /**
	 * Get the <i>x_2</i> coordinate.
	 * 
	 * @return  the <i>x_2</i> coordinate.
	 */
	public final double getY() {
		return y;
	}


	/**
	 * Set the <i>x_2</i> coordinate.
	 * 
	 * @param y value to <i>x_2</i> coordinate.
	 */
	public final void setY(double y) {
		this.y = y;
	}
	
    
	 /*	====================================================================================
     * 		PUBLIC MATHEMATICAL MANIPULATION METHODS 
     * 	====================================================================================
     */
    
    
    /**  
     * Sets the value of this tuple to the sum of itself and v.
     * @param v the other vector2d
     */  
    public final void add(vector2d v)
    { 
        this.x += v.x;
        this.y += v.y;
    }

    /**
     * Sets the value of this vector2d to the difference of vector2ds
     * v1 and v2 (this = v1 - v2).
     * @param v1 the first tuple
     * @param v2 the second tuple
     */
    
 
    /**  
     * Sets the value of this vector2d to the difference
     * of itself and v (this = this - v).
     * @param v the other tuple
     */  
    public final void sub(vector2d v)
    { 
        this.x -= v.x;
        this.y -= v.y;
    }


    /**
     * Returns the negated vector2d.
     * @return -v
     */
    public final vector2d negate(){
    
    	return new vector2d(-this.x,-this.y);
    }


    /**
     * Sets the value of this vector2d to the scalar multiplication
     * of tuple vector2d.
     * @param s the scalar value
     * @param v the source tuple
     */
    public final void scale(double s, vector2d v)
    {
	this.x = s*v.x;
	this.y = s*v.y;
    }


    /**
     * Sets the value of this vector2d to the scalar multiplication
     * of itself.
     * @param s the scalar value
     */
    public final void scale(double s)
    {
        this.x *= s;
        this.y *= s;
    }


   /**
     * Returns true if all of the data members of Tuple3d t1 are
     * equal to the corresponding data members in this Tuple3d.
     * @param t1  the tuple with which the comparison is made
     * @return  true or false
     */  
    public boolean equals(vector2d v)
    {
      try {
        return(this.x == v.x && this.y == v.y);
      }
      catch (NullPointerException e2) {return false;}
    }


	/**
	 * Returns the dot product (using a Euclidean metric) of this vector and vector v1.
	 * @param v the other vector
	 * @return the dot product of this and v
	 */
	public final double dot(vector2d v){
		
		return x*v.x + y*v.y;
	}
	
	/**
	 * Return the only non-zero coordinate of the cross product 
	 * of the vector2d with v (in that order, i.e. this cross v).
	 * @return the <i>x_3</i> component of this cross v
	 */
	public final double cross(vector2d v){
		return this.x*v.y - this.y*v.x;
	}
	
	
	/**
	 * Returns the angle in the plane of this and v between this and the vector2d v.
	 * @param the other vector
	 * @return the angle between this and v in radians
	 */
	public final double angle(vector2d v){
		return Math.acos(this.dot(v) / Math.sqrt(this.length2()*v.length2()));
	}
	
	
	/**
	 * Return the dot product of the vector2d with itself.
	 * @return the Euclidean norm
	 */
	public final double length2(){
		return this.dot(this);
	}

	
    /**
     * Normalizes this vector in place.
     */
	public final void normalize(){
		double N = Math.sqrt(this.length2());
		x /= N;
		y /= N;
	}
	
    /**
     * Sets the value of this vector2d to unit vector created from v
     * @param v the vector2d
     */
    public final void norm(vector2d v){
    	double N = Math.sqrt(v.length2());
    	this.x = v.x/N;
    	this.y = v.y/N;
    }
    
    
    /**
     * 
     * Rotates this by an angle <i>theta</i> about the vector perpendicular to 
     * the plane in which the vector lives with positive angles corresponding
     * to rotation in accordance with right-handed coordinated systems.
     * 
     * @param theta the rotation angle in radians
     */
    public final void rotate(double theta){
    	double cpx = this.x;
    	double cpy = this.y;
    	this.x =  Math.cos(theta)*cpx + Math.sin(theta)*cpy;
    	this.y = -Math.sin(theta)*cpx + Math.cos(theta)*cpy;
    }
    
	 /*	====================================================================================
     * 		PUBOIC I/O METHODS 
     * 	====================================================================================
     */
    
    /**
     * Returns a string that contains the values of this Tuple3d.
     * The form is (x,y).
     * @return the String representation
     */  
    public String toString() {
        return "(" + this.x + ", " + this.y + ")";
    }
    
    
}
