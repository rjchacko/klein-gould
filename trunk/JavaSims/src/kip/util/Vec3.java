package kip.util;

import static java.lang.Math.*;

// TODO replace with scikit.vecmath.Vector3d
public class Vec3 implements Cloneable {
	public double x, y, z;
	public Vec3(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public Vec3 clone(Vec3 that) {
		try {
			return (Vec3)super.clone();
		} catch (Exception e) {
			return null;
		}
	}
	
	public double dot(Vec3 that) {
		return x*that.x + y*that.y + z*that.z;
	}
	
	public double norm2() {
		return this.dot(this);
	}
	
	public double norm() {
		return Math.sqrt(norm2());
	}
	
	public Vec3 normalize() {
		if (norm2() == 0)
			throw new IllegalArgumentException("Can't normalize zero vector.");
		else
			return scale(1/norm());
	}
	
	public Vec3 projectOnto(Vec3 that) {
		// v1*v2 v2 / v2^2
		return that.scale(this.dot(that) / that.norm2());
	}
	
	public Vec3 scale(double a) {
		return new Vec3(x*a, y*a, z*a);
	}
	
	public Vec3 plus(Vec3 that) {
		return new Vec3(x+that.x, y+that.y, z+that.z);
	}
	
	public Vec3 minus(Vec3 that) {
		return new Vec3(x-that.x, y-that.y, z-that.z);		
	}
	
	public Vec3 cross(Vec3 v) {
		return new Vec3(y*v.z - z*v.y, -(x*v.z - z*v.x), x*v.y - y*v.x);
	}
	
	public Vec3 rotateZ(double a) {
		double c = cos(a);
		double s = sin(a);
		return new Vec3(c*x-s*y, s*x+c*y, z);
	}
	
	public String toString() {
		return "("+x+","+y+","+z+")";
	}
}
