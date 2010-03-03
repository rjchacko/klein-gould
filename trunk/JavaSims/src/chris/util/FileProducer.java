package chris.util;


public class FileProducer implements Runnable{
	
	AddFileUtil f;

    public FileProducer(FileDropBox db) {
    
        f = new AddFileUtil(db);
    }

    public void run() {
    	
    	
    	f.createAndShowGUI();
    }
    
    public AddFileUtil getAF(){
    	
    	return f;
    }
}