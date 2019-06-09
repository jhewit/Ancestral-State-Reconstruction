#------------------------CSS383_Project_2_ASR.py--------------------------------
# Author: Johnathan Hewit
# Created: 4-28-2019
# Modified: 6-9-2019
#-------------------------------------------------------------------------------
# Purpose: This script is designed to work within the context of a class project
#          for CSS383 (Bioinformatics) at the University of Washington. It
#          creates an ancestral state reconstruction for determining when the
#          aquaporin gene in various anadromous and non-anadromous fish was
#          evolved in hopes of declaring some significance of the aquaporin
#          gene and anadromy.
#-------------------------------------------------------------------------------

from ete3 import Tree
from ete3 import EvolTree
import xlrd
import random
import matplotlib.pyplot as plt
import seaborn as sns

class ASRTree:
    #Attributes
    __tree = None #Actual tree
    __sim_tree = None #Simulation tree
    ____transition_prob_anad = None
    __transition_prob_aqp3 = None
    __sim_effect_sizes = [] #List containing simulation effect sizes
    __p_value_count = 0 #Number of times an effect size is simulated => actual
    __effect_size = 0 #Actual effect size of model
    __num_of_branches = __num_anad = __num_aqp3 = __num_anad_and_aqp3 = __num_taxa = __p_value = 0
    __anadromy_lookup = dict() #Dictionary matching FASTA file names (key) to a list of taxa names and character states
    SCIENTIFIC_INDEX = 0
    COMMON_INDEX = 1
    ANAD_INDEX = 2
    AQP3_INDEX = 3
    EPSILON = 0.00000000000000000001 #Number being added to anadromy/aqp3 variables to avoid division by 0 in effect size

#Public Methods

    #--------------------------constructor--------------------------------------
    # Description: Constructs ASTree and sets default value for tree, and creates
    #              the 2D list for transition rate matrix, setting initial
    #              values to 0.
    #---------------------------------------------------------------------------
    def __init__(self):
        self.__tree = None
        self.____transition_prob_anad = [[0.0 for x in range(2)] for y in range(2)]
        self.__transition_prob_aqp3 = [[0.0 for x in range(2)] for y in range(2)]
    #end constructor

    #-----------------------------build_tree------------------------------------
    # Description: Builds phylogenetic tree from newick tree file in RAxML result.
    #---------------------------------------------------------------------------
    def build_tree(self, path):
        rax_file = open(path, "r")
        if rax_file.mode == "r":
            contents = rax_file.read()
            self.__tree = Tree(contents)
            print("\nRAxML tree imported successully.")
        else:
            print("\nRAxML tree failed to import successfully. Please check the file path and try again.")
    #end build_tree

    #-----------------------run_max_parsimony-----------------------------------
    # Description: Calls private functions for Fitch's algorithm of maximum
    #              parsimony.
    #---------------------------------------------------------------------------
    def run_max_parsimony(self): #Calls private functions for Fitch's algorithm of maximum parsimony
        if self.__tree is None:
            print("\n****************Error****************\nTree has not been imported. Please run build_tree method first.")
        else:
            self.__tree.resolve_polytomy() #Transform tree to bifurcating - does nothing if already bifurcating
            self.__down_pass()
            self.__up_pass()
            self.__clean_tree()
            self.__find_char_states()
            self.__find_transition_prob()
            self.__effect_size = self.calc_effect_size(self.__num_anad + self.EPSILON,\
            self.__num_aqp3 + self.EPSILON, self.__num_anad_and_aqp3 + self.EPSILON)
    #end run_max_parsimony

    #-----------------------------get_num_taxa----------------------------------
    # Description: Returns number of taxa.
    #---------------------------------------------------------------------------
    def get_num_taxa(self):
        return self.__num_taxa
    #end get_num_taxa

    #-----------------------------get_p_value-----------------------------------
    # Description: Returns the P-Value of the hypothesis test.
    #---------------------------------------------------------------------------
    def get_p_value(self):
        return self.__p_value
    #end get_p_value

    #--------------------------import_lookup------------------------------------
    # Description: Imports the look-up file for assigning character state
    #              changes and taxa names.
    #---------------------------------------------------------------------------
    def import_lookup(self, path): #Imports the look-up file for assigning character state changes and taxa names
        import_file = xlrd.open_workbook(path)
        file = import_file.sheet_by_index(0)
        values = list() #Local list for holding cell row information

        for row in range(1, file.nrows): #Nested loops to cover entire spreadsheet
            for col in range(file.ncols): #Creates a list of the scientific names, common names and character states for each fish in file
                if col == 0:
                    file_name = file.cell_value(row, col)
                    values.append(file_name)
                elif col == 1:
                    scientific_name = file.cell_value(row, col)
                    values.append(scientific_name)
                elif col == 2:
                    common_name = file.cell_value(row, col)
                    values.append(common_name)
                elif col == 3:
                    anadromous = int(file.cell_value(row, col))
                    values.append(anadromous)
                else:
                    aqp3 = int(file.cell_value(row, col))
                    values.append(aqp3)
                    self.__anadromy_lookup[values[0]] = values[1:]
                    values.clear()

        __num_taxa = len(self.__anadromy_lookup)
    #end import_lookup

    #----------------------------show_tree--------------------------------------
    # Description: Displays tree in console and opens an external window to
    #              interact with tree and see branch length.
    #---------------------------------------------------------------------------
    def show_tree(self):
        print(self.__tree.get_ascii(attributes=["name", "anadromy", "aqp3"], show_internal=True))
        self.__tree.show()
    #end show_tree

    #----------------------------to_string--------------------------------------
    # Description: Prints to console number of taxa and their names, as well as
    #              the number of character state changes.
    #---------------------------------------------------------------------------
    def to_string(self):
        if self.__tree == None or self.__effect_size == 0:
            return "\n****************Error****************\nTree not constructed,\
             or maximum parsimony not yet run. Please run methods and try again."

        count = 0
        asr_info = "\n\t\tTaxa\n"
        for key in self.__anadromy_lookup:
            count += 1
            asr_info += str(count) + ": " + self.__anadromy_lookup[key][self.SCIENTIFIC_INDEX]
            asr_info += " (" + self.__anadromy_lookup[key][self.COMMON_INDEX] + ")\n"
        asr_info += "\nAnadromy Character State Changes: " + str(self.__num_anad)
        asr_info += "\nAQP3 Character State Changes: " + str(self.__num_aqp3)
        return asr_info
    #end to_string

    #------------------------calc_effect_size-----------------------------------
    # Description: Public method that calculates the effect size of the ASRTree.
    #---------------------------------------------------------------------------
    def calc_effect_size(self, numOfAnad, numOfAqp3, numAnadAndAqp3):
        effect_size = ((numAnadAndAqp3/self.__num_of_branches)/((numOfAnad/self.__num_of_branches)*(numOfAqp3/self.__num_of_branches)))
        return effect_size
    #end calc_effect_size

    #-------------------------monte_carlo_sim-----------------------------------
    # Description: Public method to run n number of Monte Carlo simulations
    #              in order to test the hypothesis. Each simulation checks
    #              the ancestral node in the tree, then refers to the transition
    #              rate matrix for the probability of getting the same or a
    #              different character state.
    #---------------------------------------------------------------------------
    def monte_carlo_sim(self, num_sims):
        #Checks if there already is a simulation tree to avoid unncessary copies
        self.__p_value_count = 0 #Initialize back to 0
        self.__sim_effect_sizes.clear() #Initialize back to empty
        if self.__sim_tree is None:
            self.__sim_tree = self.__tree.copy()
        for sim in range(num_sims):
            #Set values of each count back to the EPSILON value to avoid
            #division by 0 in the effect size
            aqp3_count = self.EPSILON
            anad_count = self.EPSILON
            anad_aqp3_count = self.EPSILON
            for node in self.__sim_tree.traverse("preorder"):
                rand_num_1 = random.randint(0, 1001)
                rand_num_2 = random.randint(0, 1001)
                if not node.is_root():
                    #Check each ancestor's character state, and roll a random
                    #number against the probability of going from that state to
                    #the same or a different state based on transition matrix
                    #and assign that character state. Tally all gains
                    if node.up.anadromy == 1:
                        if (self.____transition_prob_anad[1][0]*1000) > rand_num_1:
                            node.add_feature("anadromy", 0)
                        else:
                            node.add_feature("anadromy", 1)
                            anad_count += 1
                    else:
                        if (self.____transition_prob_anad[0][1]*1000) < rand_num_1:
                            node.add_feature("anadromy", 0)
                        else:
                            node.add_feature("anadromy", 1)
                            anad_count += 1
                    if node.up.aqp3 == 1:
                        if (self.__transition_prob_aqp3[1][0]*1000) > rand_num_2:
                            node.add_feature("aqp3", 0)
                        else:
                            node.add_feature("aqp3", 1)
                            aqp3_count += 1
                    else:
                        if (self.__transition_prob_aqp3[0][1]*1000) < rand_num_2:
                            node.add_feature("aqp3", 0)
                        else:
                            node.add_feature("aqp3", 1)
                            aqp3_count += 1
                    if node.anadromy == 1 and node.aqp3 == 1:
                        anad_aqp3_count += 1
            #Calculate the effect size and store the results.
            eff_size = self.calc_effect_size(anad_count, aqp3_count, anad_aqp3_count)
            self.__sim_effect_sizes.append(eff_size)
            if eff_size >= self.__effect_size:
                self.__p_value_count += 1
        self.__p_value = (self.__p_value_count/num_sims) #Calculate and store p-value
    #end monte_carlo_sim

    #--------------------------plot_histogram-----------------------------------
    # Description: Public method to plot the histogram for testing the null
    #              hypothesis.
    #---------------------------------------------------------------------------
    def plot_histogram(self):
        plt.style.use('seaborn')
        _ = plt.hist(self.__sim_effect_sizes, bins=100)
        plt.axvline(self.__effect_size, color = 'k', linestyle = 'dashed', linewidth=1)
        plt.text(self.__effect_size + .05, 200, '   Actual Effect Size:{:.3f}'.format(self.__effect_size))
        plt.xlabel('Effect Size')
        plt.ylabel('Effect Frequency')
        plt.title('Monte Carlo Simulation Distribution')
        plt.show()
    #end plot_histogram

    #--------------------__find_transition_prob---------------------------------
    # Description: Private method that determines the transition probability
    #              of each character trait change.
    #---------------------------------------------------------------------------
    def __find_transition_prob(self):
        #Establish counter variables and traverse tree
        zero_to_one_anad = zero_to_zero_anad = one_to_zero_anad = one_to_one_anad = 0.0
        zero_to_one_aqp3 = zero_to_zero_aqp3 = one_to_zero_aqp3 = one_to_one_aqp3 = 0.0
        for node in self.__tree.traverse("postorder"):
            if not node.is_root():
                #Find Anadromy transitions
                if (node.up.anadromy is 0 and node.anadromy is 0):
                    zero_to_zero_anad += 1
                elif (node.up.anadromy is 0 and node.anadromy is 1):
                    zero_to_one_anad += 1
                elif (node.up.anadromy is 1 and node.anadromy is 0):
                    one_to_zero_anad += 1
                else:
                    one_to_one_anad += 1
                #Find AQP3 transitions
                if (node.up.aqp3 is 0 and node.aqp3 is 0):
                    zero_to_zero_aqp3 += 1
                elif (node.up.aqp3 is 0 and node.aqp3 is 1):
                    zero_to_one_aqp3 += 1
                elif (node.up.aqp3 is 1 and node.aqp3 is 0):
                    one_to_zero_aqp3 += 1
                else:
                    one_to_one_aqp3 += 1

        #Insert the probability into the appropriate matrix
        self.____transition_prob_anad[0][0] = (zero_to_zero_anad/self.__num_of_branches)
        self.____transition_prob_anad[0][1] = (zero_to_one_anad/self.__num_of_branches)
        self.____transition_prob_anad[1][1] = (one_to_one_anad/self.__num_of_branches)
        self.____transition_prob_anad[1][0] = (one_to_zero_anad/self.__num_of_branches)

        self.__transition_prob_aqp3[0][0] = (zero_to_zero_aqp3/self.__num_of_branches)
        self.__transition_prob_aqp3[0][1] = (zero_to_one_aqp3/self.__num_of_branches)
        self.__transition_prob_aqp3[1][1] = (one_to_one_aqp3/self.__num_of_branches)
        self.__transition_prob_aqp3[1][0] = (one_to_zero_aqp3/self.__num_of_branches)
    #end findTransitionProb

#Private Methods
    #---------------------------__down_pass-------------------------------------
    # Description: Private method to perform down-pass to assign character state
    #              to tips and internal nodes.
    #---------------------------------------------------------------------------
    def __down_pass(self):
        for node in self.__tree.traverse("postorder"):
            #Check for internal nodes that have been visted - marked as "Ancestor"
            if node.name is "Ancestor":
                if not node.is_root():
                    #If the parent node of the current ancestor node is unvisited,
                    #attach the character state of this node to its ancestor
                    if node.up.name is "":
                        node.up.add_feature("anadromy", node.anadromy)
                        node.up.add_feature("aqp3", node.aqp3)
                        node.up.name = "Ancestor"
                    #If the node has an intersection with its ancestor, set it
                    if node.aqp3.issubset(node.up.aqp3) or node.aqp3.issuperset(node.up.aqp3):
                        node.up.add_feature("aqp3", node.up.aqp3.intersection(node.aqp3))
                    else: #Otherwise, it's a union of two states
                        node.up.add_feature("aqp3", node.up.aqp3.union(node.aqp3))
                    #If the node has an intersection with its ancestor, set it
                    if node.anadromy.issubset(node.up.anadromy) or node.anadromy.issuperset(node.up.anadromy):
                        node.up.add_feature("anadromy", node.up.anadromy.intersection(node.anadromy))
                    else: #Otherwise, it's a union of two states
                        node.up.add_feature("anadromy", node.up.anadromy.union(node.anadromy))
            else: #Otherwise, it could be an unnamed internal node, or a terminal node
                #If it's a terminal node, grab its states from the lookup
                if node.name in self.__anadromy_lookup:
                    isAnadromous = set([self.__anadromy_lookup[node.name][self.ANAD_INDEX]])
                    isAqp3 = set([self.__anadromy_lookup[node.name][self.AQP3_INDEX]])
                    node.add_feature("anadromy", isAnadromous)
                    node.add_feature("aqp3", isAqp3)

                    if node.up.name is "": #If the internal node is not yet named, it is unvisited
                        node.up.add_feature("anadromy", isAnadromous)
                        node.up.add_feature("aqp3", isAqp3)
                        node.up.name = "Ancestor" #Tag internal nodes as Ancestor to easily identify visited nodes

                    if self.__anadromy_lookup[node.name][self.AQP3_INDEX] in node.up.aqp3:
                        node.up.add_feature("aqp3", node.aqp3.intersection(node.up.aqp3))
                    else:
                        node.up.add_feature("aqp3", node.up.aqp3.union(node.aqp3))

                    if self.__anadromy_lookup[node.name][self.ANAD_INDEX] in node.up.anadromy:
                        node.up.add_feature("anadromy", node.anadromy.intersection(node.up.anadromy))
                    else:
                        node.up.add_feature("anadromy", node.up.anadromy.union(node.anadromy))
                node.name = self.__anadromy_lookup[node.name][self.COMMON_INDEX]
    #end __down_pass

    #----------------------------__up_pass--------------------------------------
    # Description: Private method to perform up-pass to clear any union in
    #              ancestor nodes by sinding the intersection of the
    #              ancestor and its parent node.
    #---------------------------------------------------------------------------
    def __up_pass(self): #Up-pass to clear any union in ancestor nodes
        for node in self.__tree.traverse("preorder"):
            if node.name is "Ancestor":
                if not node.is_root():
                    if len(node.anadromy) > 1:
                        node.add_feature("anadromy", node.anadromy.intersection(node.up.anadromy))
                    if len(node.aqp3) > 1:
                        node.add_feature("aqp3", node.aqp3.intersection(node.up.aqp3))
    #end __up_pass

    #--------------------------__clean_tree-------------------------------------
    # Description: Private function to clear the sets in the attributes for
    #              anadromy and AQP3 in each node and turn them into integers.
    #---------------------------------------------------------------------------
    def __clean_tree(self):
        for node in self.__tree.traverse("preorder"):
            character_state_anad = next(iter(node.anadromy))
            character_state_aqp3 = next(iter(node.aqp3))
            node.add_feature("anadromy", character_state_anad)
            node.add_feature("aqp3", character_state_aqp3)
    #end __clean_tree

    #-------------------------__find_char_states---------------------------------
    # Description: Private function to find the number of branches, as well as
    #              find the number of character states - both individual and
    #              branches with both andromy and AQP3.
    #---------------------------------------------------------------------------
    def __find_char_states(self):
        for node in self.__tree.traverse("preorder"):
            self.__num_of_branches += 1
            if node.anadromy == 1 and node.aqp3 == 1:
                self.__num_anad_and_aqp3 += 1
            if node.anadromy == 1:
                self.__num_anad += 1
            if node.aqp3 == 1:
                self.__num_aqp3 += 1
        self.__num_of_branches -= 1 #Not counting the root as a separate branch
    #end __find_char_states

#end ASRTree

 #--------------------------------main------------------------------------------
 # Description: The main/driver file to support the ASRTree. Interacts with user
 #              to import files, create the ASR, and display the result.
 #------------------------------------------------------------------------------
newASR = ASRTree()
user_input = 69
print("\n\nWelcome to Anadromy Determinator 1000")
input("\nPress Enter/Return to begin")
while user_input != -1:
    print("\n\n\tMain Menu")
    user_input = int(input("\nChoose one of the following options:\n[1] Build Tree\
    \n[2] Import Look-Up File\n[3] Run Maximum Parsimony\n[4] Tree Information\
    \n[5] Display Tree\n[6] Run Monte Carlo Simulations\n[7] Show Histogram\
    \n[8] Get P-Value\n[0] Exit Program\n\n"))
    if user_input == 1:
        newASR.build_tree("RAxML_bestTree.result")
    elif user_input == 2:
        path = input("\nPlease input the file path for the look-up file, fish_anadromy.xlsx: ")
        newASR.import_lookup(path)
    elif user_input == 3:
        newASR.run_max_parsimony()
    elif user_input == 4:
        print(newASR.to_string())
    elif user_input == 5:
        newASR.show_tree()
    elif user_input == 6:
        newASR.monte_carlo_sim(1000)
    elif user_input == 7:
        newASR.plot_histogram()
    elif user_input == 8:
        print("\nP-Value:", newASR.get_p_value())
    elif user_input == 0:
        break
    else:
        print("\nInvalid Entry. Please try again.")

print("\n\nThank you for using Anadromy Determinator 1000\n\n")
