from ete3 import Tree
from ete3 import EvolTree

class ASRTree:
    charStateChanges = 0
    def __init__(self):
        self.tree = None

    def buildTree(self):
        tree = None
        return tree

#Map of Taxa to a 1 if anadromous and 0 if non-anadromous.
anadromy = {'sea_lion': 0, 'seal': 0, 'monkey': 0, 'cat': 1, 'weasel': 1, 'dog': 1, 'raccoon': 0, 'bear': 1}

tree = Tree("(tr|B2LJ88|B2LJ88_SALAL:0.00000100000050002909,(tr|B2LIV4|B2LIV4_ONCNE:0.00465520075489249049,(tr|B2LIV8|B2LIV8_ONCTS:0.00000100000050002909,(tr|B2LIU8|B2LIU8_ONCMY:0.00000100000050002909,tr|B2LIU1|B2LIU1_ONCKI:0.00000100000050002909):0.00000100000050002909):0.00000100000050002909):0.00467063957289102442,tr|B2LJ81|B2LJ81_SALSA:0.00000100000050002909):0.0;")

tree.resolve_polytomy()

# Downpass for ASR
for node in tree.traverse("postorder"):
    if node.name in anadromy:
        isAnadromous = set([anadromy[node.name]])
        node.add_feature("anadromy", isAnadromous)

        if node.up.name is "":
            node.up.add_feature("anadromy", isAnadromous)
            node.up.name = "visited"

        elif anadromy[node.name] in node.up.anadromy:
            node.up.add_feature("anadromy", node.anadromy.intersection(node.up.anadromy))

        else:
            node.up.add_feature("anadromy", node.up.anadromy.union(node.anadromy))

    if node.name is "visited":
        if not node.is_root():
            if node.up.name is "":
                node.up.add_feature("anadromy", node.anadromy)
                node.up.name = "visited"

            elif node.anadromy.issubset(node.up.anadromy) or node.anadromy.issuperset(node.up.anadromy):
                node.up.add_feature("anadromy", node.up.anadromy.intersection(node.anadromy))

            else:
                node.up.add_feature("anadromy", node.up.anadromy.union(node.anadromy))

# Up-pass for ASR
#for node in tree.traverse("preorder"):

print(tree.get_ascii(attributes=["name", "anadromy"], show_internal=True)) #internal will show internal nodes, removing name will only show attributes
tree.show()
