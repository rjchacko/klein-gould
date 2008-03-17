package kip.md.apps.shaker;

import static java.lang.Math.PI;
import static java.lang.Math.exp;
import static java.lang.Math.sqrt;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;

import scikit.dataset.Accumulator;
import scikit.dataset.DynamicArray;
import scikit.dataset.Function;
import scikit.dataset.Histogram;
import scikit.util.FileUtil;
import scikit.util.Terminal;

class Particle {
	double radius, id;
	public DynamicArray t = new DynamicArray();
	public DynamicArray x = new DynamicArray();
	public DynamicArray y = new DynamicArray();
	
	public void clear() {
		t.clear();
		x.clear();
		y.clear();
	}
	
	public void acc(double x, double y, double t) {
		this.x.append(x);
		this.y.append(y);
		this.t.append(t);
	}
	
	public double dist2(int i1, int i2) {
		double dx = x(i1) - x(i2);
		double dy = y(i1) - y(i2);
		return dx*dx+dy*dy;
	}
	
	public int size() {return t.size();}
	public double x(int i) {return x.get(i); }
	public double y(int i) {return y.get(i); }
	
}

class Data {
	int maxCnt = Integer.MAX_VALUE;
	int cnt;
	FloatBuffer fb;
	
	public Data(String fname) {
		try {
			FileChannel channel = new FileInputStream(fname).getChannel();
			MappedByteBuffer bb = channel.map(FileChannel.MapMode.READ_ONLY, 0, channel.size());
			bb.order(ByteOrder.LITTLE_ENDIAN);
			
			IntBuffer ib = bb.asIntBuffer();
		    ib.get(); // unknown meaning
		    int dim = ib.get(); // spatial dimensions
		    int n = ib.get(); // columns of data
		    int m = ib.get(); // rows of data
		    int prec = ib.get(); // precision (4 for float, 8 for double)
		    int size = ib.get(); // m*n
		    for (int i = 0; i < 6; i++)
		    	ib.get(); // these should be zero
		    assert(dim == 2);
		    assert(prec == 4);
		    assert(size == m*n);
		    
		    bb.position(4*ib.position());
		    
		    cnt = 0;
		    fb = bb.asFloatBuffer();
		    
		} catch (Exception e) {
			System.out.println(e);
		}
	}
	
	public void reset() {
		cnt = 0;
		fb.rewind();
	}
	
	public boolean hasRemaining() {
		return cnt < maxCnt && fb.hasRemaining();
	}
	
	public Particle nextParticle() {
		if (!hasRemaining())
			return null;
		
		Particle ret = new Particle();
		double x, y, radius, time, id;
    	x = fb.get();
    	y = fb.get();
    	fb.get(); // brightness
    	radius = fb.get();
    	time = fb.get();
    	id = fb.get();
    	ret.id = id;
    	ret.radius = radius;
   		ret.acc(x, y, time);
   		
    	while (fb.hasRemaining()) {
    		fb.mark();
        	x = fb.get();
        	y = fb.get();
        	fb.get(); // brightness
        	radius = fb.get();
        	time = fb.get();
        	id = fb.get();
        	if (id == ret.id)
        		ret.acc(x, y, time);
        	else {
        		fb.reset();
        		break;
        	}
    	}
    	cnt++;
		return ret;
	}
}

class Alpha {
	public Accumulator x2;
	public Accumulator x4;
	public Accumulator alpha;
}

class Commands {
	Terminal term;
	int minFrames = 500;
	int frameJump = 100;
	
	public Commands(Terminal term) {
		this.term = term;
	}
	
	public Data loadData(String fname) {
		term.println("Opening " + fname);
		return new Data(fname);
	}
	
	public Data loadData() {
		try {
			String fname = FileUtil.loadDialog(term.getConsole(), "Open Shaker Data");
			term.println("Opening " + fname);
			return new Data(fname);
		} catch(IOException e) {
			term.println(e);
			return null;
		} 
	}
	
	public Histogram trajectoryDistribution(Data data) {
		data.reset();
		Histogram len = new Histogram(100);
		while (data.hasRemaining()) {
			Particle p = data.nextParticle();
			len.accum(p.t.size());
		}
		term.println(data.cnt + " trajectories");
		return len;
	}
	
	public Alpha analyze(Data data) {
		int nsteps = 500;
		int[] steps = new int[nsteps];
		for (int i = 0; i < nsteps; i++) {
			steps[i] = (int)exp(i*0.1);
		}
		
		data.reset();
		double bw = 1;
		Accumulator x2 = new Accumulator(bw);
		Accumulator x4 = new Accumulator(bw);
		while (data.hasRemaining()) {
			Particle p = data.nextParticle();
			if (p.t.size() < minFrames) continue;
			
			for (int i = 0; i < p.size(); i += frameJump) {
				for (int s = 0; s < nsteps; s++) {
					int j = i+steps[s];
					if (j >= p.size()) break;
					double d2 = p.dist2(i, j);
					x2.accum(j-i, d2);
					x4.accum(j-i, d2*d2);
				}
			}
		}
		
		Accumulator alpha = new Accumulator(bw);		
		for (double t : x2.keys()) {
			double d = 2;
			double d2 = x2.eval(t);
			double d4 = x4.eval(t);
			alpha.accum(t, (1/(1.+2./d)) * (d4/(d2*d2)) - 1.0);
		}

		Alpha ret = new Alpha();
		ret.x2 = x2;
		ret.x4 = x4;
		ret.alpha = alpha;
		return ret;
	}
	
	public Histogram vanHove(Data data, int tstar) {
		Histogram ret = new Histogram(0.05);
		ret.setNormalizing(true);
		data.reset();
		while (data.hasRemaining()) {
			Particle p = data.nextParticle();
			if (p.t.size() < minFrames) continue;
			
			for (int i = 0; i < p.size(); i += frameJump) {
				int j = i + tstar;
				if (j >= p.size()) break;
				double r = sqrt(p.dist2(i, j));
				ret.accum(r, 2*PI*r);
			}
		}
		return ret;
	}
	
	public Function vanHoveTheory(int tstar, double diffusion) {
		final double sig2 = diffusion*tstar;
		return new Function(0, 5*sqrt(sig2)) {
			public double eval(double r) {
				return (2*PI*r)*(1/(sig2*2*PI))*exp(-(r*r)/(2*sig2));
			}
		};
	}
}

public class ShakerApp extends Terminal {
	public static void main(String[] args) {
		ShakerApp term = new ShakerApp();
		term.help = "Suggested command sequence:\n"+
			"\tdata = loadData();\n"+
			"\t// plot distribution of particle tracking duration\n"+
			"\tplot(trajectoryDistribution(data));\n"+
			"\t// neglect particles tracked for less than, e.g., 500 frames\n"+
			"\tdata.minFrames = 500; // default 500\n"+
			"\t// smaller frameJump means more averages and slower to process\n"+
			"\tdata.frameJump = 100; // default 100\n"+
			"\ta = analyze(data);\n"+
			"\tplot(a.x2);\n"+
			"\tplot(a.alpha);\n"+
			"\t// do van-hove calculation with t* = 1000\n"+
			"\tplot(vanHove(data, 1000), \"Experiment\");\n"+
			"\t// compare to gaussian with t* = 1000, diffusion coef = 0.003\n"+
			"\treplot(vanHoveTheory(1000, 0.003), \"Theory\");\n"+
			"(Right click in the plot windows to enable log scale and save data)";
		term.importObject(new Commands(term));
		term.runApplication();
	}
}
