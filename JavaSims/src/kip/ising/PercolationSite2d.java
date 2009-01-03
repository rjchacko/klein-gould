package kip.ising;

public class PercolationSite2d extends NewmanZiff {
	int lx, ly;
	
	public PercolationSite2d(int lx, int ly) {
		super(lx*ly);
		this.lx = lx;
		this.ly = ly;
	}
	
	public void occupyAndBondSites(int a[], int dir) {
		assert (a.length == lx*ly);
		for (int i = 0; i < lx*ly; i++)
			if (a[i] == dir)
				occupyAndBondSite(i);
	}
	
	public void occupyAndBondSite(int i) {
		super.occupySite(i);
		
		int x = i % lx;
		int y = i / lx;			
		int xp = (x+1)%lx;
		int xm = (x-1+lx)%lx;
		int yp = (y+1)%ly;
		int ym = (y-1+ly)%ly;
		
		tryBond(i, yp*lx+x, 0, 1);
		tryBond(i, ym*lx+x, 0, -1);
		tryBond(i, y*lx+xp, 1, 0);
		tryBond(i, y*lx+xm, -1, 0);
	}
	
	private void tryBond(int i, int j, int dx, int dy) {
		if (isOccupied(j))
			addBond(i, j, dx, dy);
	}
}
