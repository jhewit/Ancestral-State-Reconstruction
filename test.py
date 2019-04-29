from ete3 import Tree
from ete3 import EvolTree

#Map of Taxa to a 1 if anadromous and 0 if non-anadromous.
anadromy = {'sea_lion': 0, 'seal': 0, 'monkey': 0, 'cat': 1, 'weasel': 1, 'dog': 1, 'raccoon': 0, 'bear': 1}
#Creates tree object; currently just an example string
tree = Tree("((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154);")
#Resolves any polytomy; converts tree to a bifurcating/binary tree
tree.resolve_polytomy()

#Downpass for ASR
for node in tree.traverse("postorder"):
    if node.name in anadromy: #Finds the tips
        isAnadromous = set([anadromy[node.name]])
        node.add_feature("anadromy", isAnadromous)

        if node.up.name is "": #Finds an unvisited parent/ancestor node
            node.up.add_feature("anadromy", isAnadromous)
            node.up.name = "visited"

        elif anadromy[node.name] in node.up.anadromy: #If the node has been visited, check if the character state is in ancestor's set; create intersection
            node.up.add_feature("anadromy", node.anadromy.intersection(node.up.anadromy))

        else: #Otherwise, we create the union of the two states
            node.up.add_feature("anadromy", node.up.anadromy.union(node.anadromy))

    if node.name is "visited": #Finds internal/parent nodes that have already been visited and makes sure it's not the root
        if not node.is_root():
            if node.up.name is "": #If the internal node above has not been visited, simply attach the state and mark it visited
                node.up.add_feature("anadromy", node.anadromy)
                node.up.name = "visited"

            elif node.anadromy.issubset(node.up.anadromy) or node.anadromy.issuperset(node.up.anadromy): #If the state exists in either, assign the intersection
                node.up.add_feature("anadromy", node.up.anadromy.intersection(node.anadromy))

            else: #Otherwise, join the sets in a union
                node.up.add_feature("anadromy", node.up.anadromy.union(node.anadromy))

#Up-pass for ASR
for node in tree.traverse("preorder"):
    if node.name is "visited":
        if node.is_root():
            node.anadromy #do something here to select root node value...
        else:
            if len(node.anadromy) > 1: #If there are two in the set (a union), then we attempt to resolve it by taking the intersection of both internal nodes' states
                node.add_feature("anadromy", node.anadromy.intersection(node.up.anadromy))


print(tree.get_ascii(attributes=["name", "anadromy"], show_internal=True)) #internal will show internal nodes, removing name will only show attributes
tree.show() #displays a separate visual for the tree outside of command line console
