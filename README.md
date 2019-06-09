# Ancestral-State-Reconstruction
This is a Python script for creating an Ancestry State Reconstruction. It's a small part of an overall project for the Bioinformatics (CSS383) course at the University of Washington in Bothell.

# MANUAL:
**Only applies to CSS383_Project_1_ASR.py & CSS383_Project_2_ASR.py scripts**

To use the program, you'll need the Python packages listed below:
1. ETE3
2. xlrd
3. matplotlib
4. Seaborn

**Only Steps 1 - 8 are applicable to CSS383_Project_1_ASR.py**

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
directory, followed by "fish_file(1).xlsx" for Project 1, or "file_file.xlsx"
for Project 2. Press Enter/Return.

7. Press 3 and Enter/Return to execute the maximum parsimony algorithm.

8. Now you can view the results of the ancestral state reconstruction by either:

    a) Pressing 4 and Enter/Return to get general information about the result in text form.
       This will tell you the taxa being used, and the number of character state changes.

    b) Pressing 5 and Enter/Return to see a visualization of the tree in the terminal window,
       and another visualization that displays the tree and its branch lengths.

9. If you're interested in running the hypothesis test via the Markov Chain Monte Carlo method:

    a) Press 6 and Enter/Return to run the default number of Monte Carlo simulations (1000).

    b) Press 7 and Enter/Return to view the histogram produced by the results of the simulations.

    c) Press 8 and Enter/Return to view the P-Value of the hypothesis test.
