package chris.util;

public class FitUtil {
	
	private double Xdata[], Ydata[], sigma;
	private int Ndata;
	
	public FitUtil(int Ndatapts){
		
		Ndata = Ndatapts;
		
		return;
	}
	
	public double[] fit(double[] XX, double[] YY, double error, boolean ssrBool){
		
		
		if( XX.length != Ndata || YY.length != Ndata ) return null;
		
		Xdata = XX;
		Ydata = YY;
		sigma = error;
		
		if(Ydata.length != Ndata) return null;
		
		double m, dm, b, db, xs, x, y, xy, chis, ssr, ssb, foo;
		
		m    = 0;
		dm   = 0;
		b    = 0;
		db   = 0;
		xs   = 0;
		x    = 0;
		y    = 0;
		xy   = 0;
		chis = 0;
		ssr  = 0;
		ssb  = 0;
		foo  = 0;
		
		for (int jj = 0 ; jj < Ndata ; jj++){
			xs += Xdata[jj]*Xdata[jj];
			x  += Xdata[jj];
			y  += Ydata[jj];
			xy += Xdata[jj]*Ydata[jj];
		}
		
		b  = (xs*y-x*xy)/(Ndata*xs-x*x);  	// Melissinos pg 450 (NB equation for a* and b* are mixed up but the
		db = xs/(Ndata*xs-x*x);				// equation for the errors in a* and b* are correct)
		db = sigma*Math.sqrt(db);
		m  = (Ndata*xy-x*y)/(Ndata*xs-x*x); 
		dm = Ndata/(Ndata*xs-x*x);
		dm = sigma*Math.sqrt(dm);
		
		
		for (int jj = 0 ; jj < Ndata ; jj++){
			foo = Ydata[jj] - (m*Xdata[jj]+b);
			chis += foo*foo;
			//PrintUtil.printlnToFile("/Users/cserino/Desktop/DeBug.txt",jj,foo*foo);
		}
		
		ssb = chis / Ndata;
		chis = chis / (sigma*sigma);

		if(ssrBool){
			for (int jj = 0 ; jj < Ndata ; jj++){
				foo = ssb - (Ydata[jj] - (m*Xdata[jj]+b))*(Ydata[jj] - (m*Xdata[jj]+b));
				ssr +=  foo*foo;
			}

			return new double[]{m, b, dm, db, chis, ssr};
		}
		else{
			return new double[]{m, b, dm, db, chis};
		}
	}
	
}
