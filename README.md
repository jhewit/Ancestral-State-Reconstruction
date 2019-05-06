# Ancestral-State-Reconstruction
This is a Python script for creating an Ancestry State Reconstruction. It's a small part of an overall project for the Bioinformatics (CSS383) course at the University of Washington in Bothell.

# MANUAL:
**Only applies to CSS383_Project_1_ASR.py script**

To use the program, you'll need the Python packages listed below:
1. ETE3
2. xlrd

To build from Source:
1. Ensure the following files are downloaded and in the same directory:
    a) CSS383_Project_1_ASR.py
    b) fish_anadromy.xlsx
    c) RAxML_bestTree.result

2. Open CommandLine/PowerShell, Terminal or Linux Terminal

3. Type "python CSS383_Project_1_ASR.py" in the command line

4. Press Enter/Return to enter the Main Menu

5. Press 1 and Enter/Return to Build the Tree (from here, you can view the original tree
before running maximum parsimony by skipping to Step )

6. Press 2 and Enter/Return to import the Look-Up file containing the fish information
and character states. You will need to manually enter the file path for the working
directory, followed by "fish_anadromy.xlsx" and press Enter/Return.

7. Press 3 and Enter/Return to execute the maximum parsimony algorithm.

8. Now you can view the results of the ancestral state reconstruction by either:
    a) Pressing 4 and Enter/Return to get general information about the result in text form.
       This will tell you the taxa being used, and the number of character state changes.
    b) Pressing 5 and Enter/Return to see a visualization of the tree in the terminal window,
       and another visualization that displays the tree and its branch lengths.
