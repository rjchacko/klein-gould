package chris.util;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;


public class AddFileUtil extends JPanel implements ActionListener{


	private static final long serialVersionUID = 1L;

	JButton addButton, doneButton;
	JFileChooser fc;
	FileDropBox db;

	public AddFileUtil(FileDropBox db){

		this.db = db;
		
		fc = new JFileChooser();
		addButton = new JButton("Add File");
		addButton.addActionListener(this);

		//Create the save button.  We use the image from the JLF
		//Graphics Repository (but we extracted it from the jar).
		doneButton = new JButton("Done");
		doneButton.addActionListener(this);

		//For layout purposes, put the buttons in a separate panel
		JPanel buttonPanel = new JPanel(); //use FlowLayout
		buttonPanel.add(addButton);
		buttonPanel.add(doneButton);

		//Add the buttons and the log to this panel.
		add(buttonPanel, BorderLayout.PAGE_START);
	}

	public void actionPerformed(ActionEvent e) {

		int returnVal;

		//Handle open button action.
		if (e.getSource() == addButton) {
			returnVal = fc.showOpenDialog(AddFileUtil.this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc.getSelectedFile();
				//This is where a real application would open the file.
				System.out.println(file);
				db.put(file);
			}

		}
		else if (e.getSource() == doneButton) {
//			returnVal = fc.showSaveDialog(AddFileUtil.this);
//			if (returnVal == JFileChooser.APPROVE_OPTION) {
//				File file = fc.getSelectedFile();
//				//This is where a real application would save the file.
//			}
			// store relevant data to pass
			db.put(null);
			
			// destroy GUI
			System.exit(0);
		}

	}

    public void createAndShowGUI() {
        //Create and set up the window.
        JFrame frame = new JFrame("JWSFileChooserDemo");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        //Add content to the window.
        frame.add(new AddFileUtil(db));

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }
    
    public FileDropBox getdb(){
    	
    	return db;
    }
}