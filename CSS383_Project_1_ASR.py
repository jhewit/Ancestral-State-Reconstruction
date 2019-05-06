from ete3 import Tree
from ete3 import EvolTree
import xlrd

class ASRTree:
    #Attributes
    __tree = None
    __charStateChanges = 0
    __numOfTaxa = 0
    __anadromyLookUp = dict() #Dictionary matching FASTA file names (key) to a list of taxa names and character states
    scientificIndex = 0
    commonIndex = 1
    stateIndex = 2

    #Constructor
    def __init__(self):
        self.__tree = None

    #Public Methods
    def buildTree(self, path): #Builds newick tree from an aligned and filtered FASTA file
        raxFile = open(path, "r")
        if raxFile.mode == "r":
            contents = raxFile.read()
            self.__tree = Tree(contents)
            print("\nRAxML tree imported successully.")
        else:
            print("\nRAxML tree failed to import successfully. Please check the file path and try again.")

    def runMaxParsimony(self): #Calls private functions for Fitch's algorithm of maximum parsimony
        if self.__tree is None:
            print("\n****************Error****************\nTree has not been imported. Please run buildTree method first.")
        else:
            self.__tree.resolve_polytomy() #Transform tree to bifurcating - does nothing if already bifurcating
            self.__downPass()
            self.__upPass()
            self.__findCharStateChanges()

    def getNumOfTaxa(self): #Returns the number of taxa
        return self.__numOfTaxa

    def getCharStateChanges(self): #Returns the number of character state changes
        return self.__charStateChanges

    def importLookUp(self, path): #Imports the look-up file for assigning character state changes and taxa names
        importFile = xlrd.open_workbook(path)
        file = importFile.sheet_by_index(0)
        values = list() #Local list for holding cell row information

        for row in range(1, file.nrows): #Nested loops to cover entire spreadsheet
            for col in range(file.ncols): #Creates a list of the scientific names, common names and character states for each fish in file
                if col == 0:
                    fileName = file.cell_value(row, col)
                    values.append(fileName)
                elif col == 1:
                    scientificName = file.cell_value(row, col)
                    values.append(scientificName)
                elif col == 2:
                    commonName = file.cell_value(row, col)
                    values.append(commonName)
                else:
                    anadromous = int(file.cell_value(row, col))
                    values.append(anadromous)
                    self.__anadromyLookUp[values[0]] = values[1:]
                    values.clear()

        __numOfTaxa = len(self.__anadromyLookUp)

    def showTree(self):
        print(self.__tree.get_ascii(attributes=["name", "anadromy"], show_internal=True))
        self.__tree.show()
        print("Char State Changes:", self.__charStateChanges)

    def toString(self):
        if self.__tree == None or self.__charStateChanges == 0:
            return "\n****************Error****************\nTree not constructed, or maximum parsimony not yet run. Please run methods and try again."

        count = 0
        asrInfo = "\n\t\tTaxa\n"
        for key in self.__anadromyLookUp:
            count += 1
            asrInfo += str(count) + ": " + self.__anadromyLookUp[key][self.scientificIndex]
            asrInfo += " (" + self.__anadromyLookUp[key][self.commonIndex] + ")\n"
        asrInfo += "\nCharacter States Changes: " + str(self.__charStateChanges)
        return asrInfo

    #Private Methods
    def __downPass(self): #Down-pass to assign character state to tips and internal nodes
        for node in self.__tree.traverse("postorder"):
            #Check for internal nodes that have been visted - marked as "Ancestor"
            if node.name is "Ancestor":
                if not node.is_root():
                    if node.up.name is "":
                        node.up.add_feature("anadromy", node.anadromy)
                        node.up.name = "Ancestor"

                    elif node.anadromy.issubset(node.up.anadromy) or node.anadromy.issuperset(node.up.anadromy):
                        node.up.add_feature("anadromy", node.up.anadromy.intersection(node.anadromy))

                    else:
                        node.up.add_feature("anadromy", node.up.anadromy.union(node.anadromy))
            else:
                if node.name in self.__anadromyLookUp:
                    isAnadromous = set([self.__anadromyLookUp[node.name][self.stateIndex]])
                    node.add_feature("anadromy", isAnadromous)

                    if node.up.name is "": #If the internal node is not yet named, it is unvisited
                        node.up.add_feature("anadromy", isAnadromous)
                        node.up.name = "Ancestor" #Tag internal nodes as Ancestor to easily identify visited nodes

                    elif self.__anadromyLookUp[node.name][self.stateIndex] in node.up.anadromy:
                        node.up.add_feature("anadromy", node.anadromy.intersection(node.up.anadromy))

                    else:
                        node.up.add_feature("anadromy", node.up.anadromy.union(node.anadromy))
                node.name = self.__anadromyLookUp[node.name][self.commonIndex]

    def __upPass(self): #Up-pass to clear any union in ancestor nodes
        for node in self.__tree.traverse("preorder"):
            if node.name is "Ancestor":
                if not node.is_root():
                    if len(node.anadromy) > 1:
                        node.add_feature("anadromy", node.anadromy.intersection(node.up.anadromy))

    def __findCharStateChanges(self):
        characterState = 0
        for node in self.__tree.traverse("preorder"):
            if node.is_root():
                characterState = next(iter(node.anadromy))
            else:
                if not (characterState in node.anadromy):
                    self.__charStateChanges += 1
#Main - Driver for ASRTree
newASR = ASRTree()
userInput = 69
print("\n\nWelcome to Anadromy Determinator 1000")
input("\nPress any key to begin")
while userInput != -1:
    print("\n\n\tMain Menu")
    userInput = int(input("\nChoose one of the following options:\n[1] Build Tree\n[2] Import Look-Up File\n[3] Run Maximum Parsimony\n[4] Tree Information\n[5] Display Tree\n[0] Exit Program\n\n"))
    if userInput == 1:
        newASR.buildTree("RAxML_bestTree.result")
    elif userInput == 2:
        path = input("\nPlease input the file path for the look-up file, fish_anadromy.xlsx: ")
        newASR.importLookUp(path)
    elif userInput == 3:
        newASR.runMaxParsimony()
    elif userInput == 4:
        print(newASR.toString())
    elif userInput == 5:
        newASR.showTree()
    elif userInput == 0:
        break
    else:
        print("\nInvalid Entry. Please try again.")

print("\n\nThank you for using Anadromy Determinator 1000\n\n")
