package rachele.damage2D;

import kip.util.Random;

public class AbstractCG_OFC {

	public int L, N, R, Lp, Np, dx, dt, noNbors;
	boolean fullyConnected;
	Random random = new Random();
	
	int findCG_site(int s){
		int x = s % L;
		int y = s / L;
		int xp = x / dx;
		int yp = y / dx;
		int cgSite = yp *Lp + xp;
		return cgSite;
	}
	
	
	
}
