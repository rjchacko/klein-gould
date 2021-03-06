package chris.util;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;


public class movieUtil {
	
	private static DecimalFormat fmt = new DecimalFormat("000000");


	/**
	 * Scales the image up by s as to make larger images for movies
	 * 
	 * @param data the array storing the data
	 * @param Lx the number of horizontal entries in the grid
	 * @param Ly the number of vertical entries in the grid
	 * @param s scaling factor to produce a larger image
	 * @param g a grid on which to create the image of data.
	 * @param pth a string which is the path to the location to save the image
	 * @param picID a number to uniquely label the picture
	 */
	public static void saveImage(int[] data, int Lx, int Ly, int s, Grid g, String pth, int picID){
		
		int sLx = s*Lx;
		int sLy = s*Ly;

		int[] sd = new int[sLx*sLy];
		for(int jj = 0 ; jj < sLy; jj++){
			for(int kk = 0 ; kk < sLx; kk++){
				sd[kk+jj*sLx] = data[(int)(kk/s)+(int)(jj/s)*Lx];
			}
		}
		
		g.registerData(sLx,sLy,sd);
		String SaveAs = pth + File.separator + g.getTitle()+fmt.format(picID)+".png";
		try {
			ImageIO.write(g.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		return;
	}
	
	public static void saveLegend(int min, int max, ColorPalette cp, String pth){
		
		max++; // Include max in legend
		
		int[] legend = new int[max-min];
		Grid g = new Grid("Legend");
		
		for(int jj = 0 ; jj < max-min ; jj++){
			legend[jj] = jj;
		}
		g.setColors(cp);
		
		saveImage(legend, 1, max-min, 60, g, pth, 0);
		return;
	}
}
	
	
