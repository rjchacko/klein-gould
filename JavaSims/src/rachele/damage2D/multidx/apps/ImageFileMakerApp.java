package rachele.damage2D.multidx.apps;

import java.io.*;
//import java.util.*;
//import java.io.IOException;


import rachele.util.*;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;

public class ImageFileMakerApp extends Simulation{

	String parentDir;
	int L;
	String R;
	int percentDead;
	String damageType;
	
	public static void main(String[] args) {
		new Control(new ImageFileMakerApp(), "Image File Maker");
	}

	public void animate() {
	}

	public void clear() {
	}

	public void load(Control c) {
		params.add("Data Dir",new DirectoryValue("/Users/erdomi/data/damage/contract2/testRuns"));
		params.add("Image Type", new ChoiceValue(".jpg", ".eps", ".pdf"));
	}

	public void run() {
		parentDir = params.sget("Data Dir") + File.separator;
		p("Parent Directory = ", parentDir);
		readParams();
		makeShFiles();

		p("Done");
	}
	
	void makeShFiles(){
		//dd (dead dissipation files)
		String ddbrshFile = "dd" + File.separator + "br" + File.separator + "sh.txt";
//		String ddrdshFile = "dd" + File.separator + "rd" + File.separator + "sh.txt";
//		String dddsshFile = "dd" + File.separator + "ds" + File.separator + "sh.txt";
//		String dddbshFile = "dd" + File.separator + "db" + File.separator + "sh.txt";

		//nd (no dissipation files)
		String ndbrshFile = "nd" + File.separator + "br" + File.separator + "sh.txt";
//		String ndrdshFile = "nd" + File.separator + "rd" + File.separator + "sh.txt";
//		String nddsshFile = "nd" + File.separator + "ds" + File.separator + "sh.txt";
//		String nddbshFile = "nd" + File.separator + "db" + File.separator + "sh.txt";
		
		//make new pyxplot files
		String BRshFile = parentDir + "pyxplotShBR.txt";
		FileUtil.initFile(BRshFile, params);
		makeDissipationCompareFiles(BRshFile, "output", ddbrshFile, ndbrshFile);
		FileUtil.printlnToFile(BRshFile, "");
		
	}
	
	void makeDissipationCompareFiles(String newFile, String outputFile, String ddFile, String ndFile){
		setOutput(newFile, outputFile);
		sp(newFile);
		FileUtil.printlnToFile(newFile, "set multiplot");
		FileUtil.printlnToFile(newFile, "set grid x1y1x2");
		FileUtil.printlnToFile(newFile, "set width 10");
		FileUtil.printlnToFile(newFile, "set key below");
		FileUtil.printlnToFile(newFile, "set keycolumns 2");
		FileUtil.printlnToFile(newFile, "set log y");
		FileUtil.printlnToFile(newFile, "set log x");
		sp(newFile);
		makeSimpleTitle(newFile);
		FileUtil.printlnToFile(newFile, "text", "'Damage = " + damageType +"'", "at 3,3.2");
		sp(newFile);
		FileUtil.printlnToFile(newFile, "set xlabel 'Size'");
		FileUtil.printlnToFile(newFile, "set xlabel 'Frequency'");
		sp(newFile);
		FileUtil.printlnToFile(newFile, "set origin 2,0");
		sp(newFile);
		String plotLine = "plot '" + ddFile + "' t 'Dead Site Dissipation', '" + ndFile + "' t 'No Dead Dissipation'";
		FileUtil.printlnToFile(newFile, plotLine);
	}
	
	void makeSimpleTitle(String fileName){
		FileUtil.printlnToFile(fileName, "# Make title");
		String titleLine = "set title '$L=" + L + "$, $R=" + R + "$, $" + percentDead +  " $ % damaged sites '"; 
		FileUtil.printlnToFile(fileName, titleLine);
		p("Title line = " ,titleLine);
	}
	
	void setOutput(String fileName, String outFile){
		String imageType = params.sget("Image Type");
		if(imageType==".jpg")
			FileUtil.printlnToFile(fileName, "set term jpg color");
		else if(imageType==".eps")
			FileUtil.printlnToFile(fileName, "set term postscript color");
		else if(imageType==".pdf")
			FileUtil.printlnToFile(fileName, "set term pdf color");
		String outFileName = outFile + imageType;
		FileUtil.printlnToFile(fileName, "set output", outFileName);
	}
		
	
	void readParams(){
		String infoFile = parentDir + "dd" + File.separator + "br" + File.separator + "info.txt";
		
		String sPow = ReadFile.readWord(infoFile, 4, 3);
		int pow = Integer.parseInt(sPow);
		L = (int)Math.pow(2, pow);
		p("L = " , L);
		
		R = ReadFile.readWord(infoFile, 5, 2);
		p("R = ", R);
		
		//find percent damaged sites
		String sDamagedSites = ReadFile.readWord(infoFile, 22, 5);
		p("no damaged sites string = ", sDamagedSites);
		double noDamagedSites = Double.parseDouble(sDamagedSites);
		p("no damaged sites double = ", noDamagedSites);
		double size = (double)L*(double)L;
		double  dPercentDead = 100.0*noDamagedSites/size;
		if((dPercentDead - (double)((int)dPercentDead))>0.5){
			dPercentDead += 1;
			percentDead = (int)dPercentDead;
		}else{
			percentDead = (int)dPercentDead;
		}
		
		//find damage type
		p("Percent damaged sites = ", percentDead);
		damageType = ReadFile.readWord(infoFile, 1, 4) + " " +ReadFile.readWord(infoFile, 1, 5);
		p("Damage Type = ", damageType);
		
	}

	public static String fileToString(String fileName) {
        File file = new File(fileName);
        StringBuilder contents = new StringBuilder();
        BufferedReader input = null;
        try {
            input = new BufferedReader(new FileReader(file));
            String line = null;
            while((line = input.readLine()) != null)
                contents.append(line);
            return contents.toString();
        } catch (Exception e) {
            throw new RuntimeException(e);
        } finally {
            try {
                if(input != null) input.close();
            } catch(Exception e) {
            }
        }
    }

	
	void sp(String fileName){
		FileUtil.printlnToFile(fileName, " ");
	}
	
	
	void p(String t){
		System.out.println(t);
	}
	void p(String t, String t2){
		System.out.println(t + t2);
	}
	void p(String t, Double d){
		System.out.println(t + d);
	}
	void p(String t, int i){
		System.out.println(t + i);
	}	
	void p(String t, boolean i){
		System.out.println(t + i);
	}
}
