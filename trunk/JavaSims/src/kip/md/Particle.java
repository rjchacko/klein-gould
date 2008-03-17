package kip.md;

import kip.util.Vec3;
import scikit.util.Point;


public class Particle extends Point {
	public double vx = 0, vy = 0, vz = 0;
	public ParticleTag tag;
	
	public void force(Particle that, Vec3 f) {
		f.x = f.y = f.z = 0;
	}
	public void force(Vec3 f) {
		f.x = f.y = f.z = 0;
	}
	public double potential(Particle that) {
		return 0;
	}
	public double potential() {
		return 0;
	}
}
