package kip.md.apps.shaker;

import kip.md.Particle;
import kip.md.ParticleContext;

public interface AbstractTrajectory {
	public double startTime();
	public double endTime();
	public Particle[] get(double t);	
	public ParticleContext getContext();
}
