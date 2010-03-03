package chris.tests;

import chris.util.AddFileUtil;
import chris.util.FileDropBox;
import chris.util.Random;

public class Producer implements Runnable{
    private Drop drop;

    public Producer(Drop drop) {
        this.drop = drop;
        
    }

    public void run() {
    	AddFileUtil f = new AddFileUtil(new FileDropBox());
    	f.createAndShowGUI();
    	
        String importantInfo[] = {
            "Mares eat oats",
            "Does eat oats",
            "Little lambs eat ivy",
            "A kid will eat ivy too"
        };
        Random random = new Random();

        for (int i = 0; i < importantInfo.length; i++) {
            drop.put(importantInfo[i]);
            try {
                Thread.sleep(random.nextInt(5000));
            } catch (InterruptedException e) {}
        }
        System.out.println("DONE");
        drop.alertDone();
    }
}