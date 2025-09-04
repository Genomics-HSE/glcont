
from anytree.importer import JsonImporter
from anytree import Node
from anytree.search import find
from anytree import find_by_attr, PreOrderIter, Walker

def find_by_name(name, tree):
    return find(tree, lambda node : node.name == name)


def calc_genomic_distance(node1: str, node2: str, tree) -> int:
    """
    Calculate the genomic distance between two nodes in the phylogenetic tree.
    walk_paths = w.walk(node1, node2)
    for node in walk_paths[0]:
    node1 (Node): The first node.
    node2 (Node): The second node.

    Returns:
    returns number of mutations (snps and indels) betwee two haplogroup
    """
    
    
    node1 = find_by_name(node1, tree)
    node2 = find_by_name(node2, tree)
    w = Walker()
    mutations = list()
    # if node1.is_root or node2.is_root:
        # print('AAAA')
    walker = w.walk(node1, node2)
    for node in walker[0]:
        mutations += node.snps
        mutations += node.snps_back
        mutations += node.insertion
        mutations += node.deletion
    for node in walker[2]:
        mutations += node.snps
        mutations += node.snps_back
        mutations += node.insertion
        mutations += node.deletion
            
    return len(mutations)


def import_tree(json_fn):
    """
    Make tree from json file

    Parameters:
    json_fn : str
    path to json file

    Returns:
    root : Node
    root node of the tree
    """
    try:
        with open(json_fn) as f:
            json_text = f.read()
    except OSError:
        return "No such tree file"

    importer = JsonImporter()
    root = importer.import_(json_text)
    return root