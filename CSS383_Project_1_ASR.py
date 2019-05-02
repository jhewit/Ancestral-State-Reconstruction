from ete3 import Tree
from ete3 import EvolTree
import xlrd

class ASRTree:
    #Attributes
    __charStateChanges = 0
    __numOfTaxa = 0
    __anadromyLookUp = dict()
    __tree = None

    #Constructor
    def __init__(self):
        self.__tree = None

    #Public Methods
    def buildTree(self, path):
        raxFile = open(path, "r")
        if raxFile.mode == "r":
            contents = raxFile.read()
            self.__tree = Tree(contents)
            print("\nRAxML tree imported successully.")
        else:
            print("\nRAxML tree failed to import successfully. Please check the file path and try again.")

    def runMaxParsimony(self):
        if self.__tree is None:
            print("\n****************Error****************\nTree has not been imported. Please run buildTree method first.")
        else:
            self.__tree.resolve_polytomy()
            self.__downPass()
            self.__upPass()

    def getNumOfTaxa(self):
        return self.__numOfTaxa

    def getCharStateChanges(self):
        return self.__charStateChanges

    def importLookUp(self, path):
        importFile = xlrd.open_workbook(path)
        file = importFile.sheet_by_index(0)
        values = list()

        for row in range(1, file.nrows):
            for col in range(file.ncols):
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

    #Private Methods
    def __downPass(self):
        for node in self.__tree.traverse("postorder"):
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
                    isAnadromous = set([self.__anadromyLookUp[node.name][2]])
                    node.add_feature("anadromy", isAnadromous)

                    if node.up.name is "":
                        node.up.add_feature("anadromy", isAnadromous)
                        node.up.name = "Ancestor"

                    elif self.__anadromyLookUp[node.name][2] in node.up.anadromy:
                        node.up.add_feature("anadromy", node.anadromy.intersection(node.up.anadromy))

                    else:
                        node.up.add_feature("anadromy", node.up.anadromy.union(node.anadromy))
                node.name = self.__anadromyLookUp[node.name][1]

    def __upPass(self):
        for node in self.__tree.traverse("preorder"):
            if node.name is "Ancestor":
                if not node.is_root():
                    if len(node.anadromy) > 1:
                        node.add_feature("anadromy", node.anadromy.intersection(node.up.anadromy))

newASR = ASRTree()
path = ("C:/Users/johna/Google Drive/School/UW Bothell/CSS 383/Bioinformatics Team/Project 1/Python Resources/fish_anadromy.xlsx")
newASR.importLookUp(path)
newASR.buildTree("RAxML_bestTree.result")
newASR.runMaxParsimony()
newASR.showTree()
