package chris.util;

import java.io.File;

public class FileDropBox {
    //True if consumer should wait for producer to send message, false
    //if producer should wait for consumer to retrieve message.
    private boolean empty = true;
    private boolean done  = false;
    private File file;

    public synchronized File take() {
        //Wait until message is available.
        while (empty) {
            try {
                wait();
            } catch (InterruptedException e) {}
        }
        //Toggle status.
        empty = true;
        //Notify producer that status has changed.
        notifyAll();
        return file;
    }

    public synchronized void put(File file) {
        //Wait until message has been retrieved.
        while (!empty) {
            try { 
                wait();
            } catch (InterruptedException e) {}
        }
        //Toggle status.
        empty = false;
        //Store message.
        this.file = file;
        //Notify consumer that status has changed.
        notifyAll();
    }
    
    public void alertDone(){
    	
    	done = true;
    }
    
    public boolean isDone(){
    	
    	return done;
    }
}