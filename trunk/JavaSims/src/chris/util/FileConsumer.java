package chris.util;

import java.io.File;


public class FileConsumer implements Runnable{
	private FileDropBox drop;
	File fin;

	public FileConsumer(FileDropBox drop) {
		this.drop = drop;
	}

	public void run() {
		
		for (File file = drop.take(); ! file.equals(null);
		file = drop.take()) {
			System.out.format("FILE RECEIVED: %s%n", file);
			fin = file;
		}
		System.out.println("DONE");
	}
	
	public File getFile(){

		return fin;
	}
	
	public FileDropBox getdb(){
		
		return drop;
	}
}


