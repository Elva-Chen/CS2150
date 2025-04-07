#!/usr/bin/env python3
from Bio import Phylo
from itertools import combinations
from treeswift import read_tree_newick, read_tree_nexml, read_tree_nexus
import dendropy
import ete3
import numpy
from dendropy import TaxonNamespace
import re
NA = "NA" # what to print when function is not implemented

# get memory usage
def memory():
    from os import getpid; from psutil import Process
    return Process(getpid()).memory_info().rss

# distance matrix in ETE Toolkit
# obtained from https://github.com/linsalrob/EdwardsLab/blob/master/trees/tree_to_cophenetic_matrix.py
def distance_matrix_ete(tree):
    leaves = tree.get_leaves()
    paths = {x:set() for x in leaves}
    for n in leaves:
        if n.is_root():
            continue
        movingnode = n
        while not movingnode.is_root():
            paths[n].add(movingnode); movingnode = movingnode.up
    leaf_distances = {x.name:{} for x in leaves}
    for (leaf1, leaf2) in combinations(leaves, 2):
        uniquenodes = paths[leaf1] ^ paths[leaf2]
        distance = sum(x.dist for x in uniquenodes)
        leaf_distances[leaf1.name][leaf2.name] = distance
        leaf_distances[leaf2.name][leaf1.name] = distance
    return leaf_distances

# main code
from sys import argv
from time import time
if len(argv) != 4 or argv[2] not in {'treeswift_newick', 'treeswift_nexml', 'treeswift_nexus'}:
    print("USAGE: %s <tree_file> <treeswift_tool> <task>" % argv[0]); exit(1)

# Load tree file
if argv[1].lower().endswith('.gz'):
    from gzip import open as gopen
    treestr = gopen(argv[1]).read().decode().strip()
else:
    treestr = open(argv[1]).read().strip()

# convert to format    
def newick_to_treeswift_nexml(newick_str):
    """Convert Newick to treeswift-compatible NeXML"""
    # Parse Newick with explicit taxa
    taxon_namespace = TaxonNamespace()
    tree = dendropy.Tree.get(
        data=newick_str,
        schema="newick",
        taxon_namespace=taxon_namespace)
    
    # Generate minimal NeXML
    nexml_lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<nexml xmlns="http://www.nexml.org/2009">',
        '  <otus id="taxa">'
    ]
    
    # Add OTUs (taxa)
    for taxon in taxon_namespace:
        nexml_lines.append(f'    <otu id="t{taxon.label}" label="{taxon.label}"/>')
    
    nexml_lines.extend([
        '  </otus>',
        '  <trees id="trees" otus="taxa">',
        '    <tree id="tree1">'
    ])
    
    # Add nodes with explicit root
    node_id_map = {}
    for i, node in enumerate(tree.preorder_node_iter()):
        node_id = f"n{i}"
        node_id_map[node] = node_id
        is_root = (node == tree.seed_node)
        label = node.taxon.label if node.taxon else None
        nexml_lines.append(
            f'      <node id="{node_id}"' + 
            (f' label="{label}"' if label else '') +
            (' root="true"' if is_root else '') + 
            '/>'
        )
    
    # Add root edge length
    if tree.seed_node.edge_length is not None:
        nexml_lines.append(
            f'      <rootedge target="{node_id_map[tree.seed_node]}" ' +
            f'length="{tree.seed_node.edge_length}"/>'
        )
    
    # Add all other edges
    for node in tree.preorder_node_iter():
        if node.parent_node is not None:
            nexml_lines.append(
                f'      <edge source="{node_id_map[node.parent_node]}" '
                f'target="{node_id_map[node]}" '
                f'length="{node.edge_length}"/>'
            )
    
    # Close tags
    nexml_lines.extend([
        '    </tree>',
        '  </trees>',
        '</nexml>'
    ])
    
    return '\n'.join(nexml_lines)

def count_taxa(newick_str):
    """Count the number of taxa in a Newick string"""
    # Remove comments and whitespace
    clean_str = re.sub(r'\[.*?\]', '', newick_str)  # Remove comments
    clean_str = re.sub(r'\s+', '', clean_str)  # Remove whitespace
    
    # Count unique taxa labels (anything before : or , or )
    taxa = set()
    current_label = []
    for char in clean_str:
        if char in [':', ',', ')', ';']:
            if current_label:
                taxa.add(''.join(current_label))
                current_label = []
        elif char not in ['(', ';']:
            current_label.append(char)
    return len(taxa)


def newick_to_nexus(newick_str, tree_name="TREE_1", translate_dict=None, taxlabels=None):
    """
    Convert a Newick string to Nexus format compatible with treeswift's read_tree_nexus()
    
    Args:
        newick_str (str): Newick format tree string
        tree_name (str): Name to give the tree in the Nexus file
        translate_dict (dict): Optional translation dictionary {id: label}
        taxlabels (list): Optional list of taxon labels
        
    Returns:
        str: Nexus format string
    """
    # Validate Newick string ends with semicolon
    if not newick_str.strip().endswith(';'):
        newick_str = newick_str + ';'
    
    # Build Nexus file content
    nexus_lines = [
        "#NEXUS",
        "BEGIN TAXA;",
        f"DIMENSIONS NTAX={count_taxa(newick_str) if taxlabels is None else len(taxlabels)};",
    ]
    
    # Add taxlabels if provided
    if taxlabels is not None:
        nexus_lines.append("TAXLABELS")
        nexus_lines.extend([f"    {label}" for label in taxlabels])
        nexus_lines.append(";")
    
    nexus_lines.append("END;")
    nexus_lines.append("BEGIN TREES;")
    
    # Add translate block if provided
    if translate_dict is not None:
        nexus_lines.append("TRANSLATE")
        translate_items = [f"    {k} {v}" for k, v in translate_dict.items()]
        # Add commas to all but last item
        translate_items = [item + "," for item in translate_items[:-1]] + translate_items[-1:]
        nexus_lines.extend(translate_items)
        nexus_lines.append(";")
    
    # Add the tree
    nexus_lines.append(f"TREE {tree_name} = {newick_str}")
    nexus_lines.append("END;")
    
    return "\n".join(nexus_lines)


nexml_tree = newick_to_treeswift_nexml(treestr)
nexus_tree = newick_to_nexus(treestr)


# ladderize
def ladderize(m):
    if m == 'treeswift_newick':
        tree = read_tree_newick(treestr)
        t_start = time()
        tree.ladderize()
        t_end = time()
    elif m == 'treeswift_nexml':
        tree = read_tree_nexml(nexml_tree)
        t_start = time()
        tree.ladderize()
        t_end = time()
    elif m == 'treeswift_nexus':
        tree = read_tree_nexus(nexus_tree)
        t_start = time()
        tree.ladderize()
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# distance matrix
def distance_matrix(m):
    if m == 'treeswift_newick':
        tree = read_tree_newick(treestr)
        t_start = time()
        tree.distance_matrix()
        t_end = time()
    elif m == 'treeswift_nexml':
        tree = read_tree_nexml(nexml_tree)
        t_start = time()
        tree.distance_matrix()
        t_end = time()
    elif m == 'treeswift_nexus':
        tree = read_tree_nexus(nexus_tree)
        t_start = time()
        tree.distance_matrix()
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# inorder traversal
def inorder(m):
    if m == 'treeswift_newick':
        tree = read_tree_newick(treestr)
        t_start = time()
        for node in tree.traverse_inorder():
            pass
        t_end = time()
    elif m == 'treeswift_nexml':
        tree = read_tree_nexml(nexml_tree)
        t_start = time()
        for node in tree.traverse_inorder():
            pass
        t_end = time()
    elif m == 'treeswift_nexus':
        tree = read_tree_nexus(nexus_tree)
        t_start = time()
        for node in tree.traverse_inorder():
            pass
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# level-order traversal
def levelorder(m):
    if m == 'treeswift_newick':
        tree = read_tree_newick(treestr)
        t_start = time()
        for node in tree.traverse_levelorder():
            pass
        t_end = time()
    elif m == 'treeswift_nexml':
        tree = read_tree_nexml(nexml_tree)
        t_start = time()
        for node in tree.traverse_levelorder():
            pass
        t_end = time()
    elif m == 'treeswift_nexus':
        tree = read_tree_nexus(nexus_tree)
        t_start = time()
        for node in tree.traverse_levelorder():
            pass
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# MRCA
def mrca(m):
    if m == 'treeswift_newick':
        tree = read_tree_newick(treestr)
        t_start = time()
        leaves = {str(l) for l in tree.traverse_leaves()}
        tree.mrca(leaves)
        t_end = time()
    elif m == 'treeswift_nexml':
        tree = read_tree_nexml(nexml_tree)
        t_start = time()
        leaves = {str(l) for l in tree.traverse_leaves()}
        tree.mrca(leaves)
        t_end = time()
    elif m == 'treeswift_nexus':
        tree = read_tree_nexus(nexus_tree)
        t_start = time()
        leaves = {str(l) for l in tree.traverse_leaves()}
        tree.mrca(leaves)
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# postorder traversal
def postorder(m):
    if m == 'treeswift_newick':
        tree = read_tree_newick(treestr)
        t_start = time()
        for node in tree.traverse_postorder():
            pass
        t_end = time()
    elif m == 'treeswift_nexml':
        tree = read_tree_nexml(nexml_tree)
        t_start = time()
        for node in tree.traverse_postorder():
            pass
        t_end = time()
    elif m == 'treeswift_nexus':
        tree = read_tree_nexus(nexus_tree)
        t_start = time()
        for node in tree.traverse_postorder():
            pass
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# preorder traversal
def preorder(m):
    if m == 'treeswift_newick':
        tree = read_tree_newick(treestr)
        t_start = time()
        for node in tree.traverse_preorder():
            pass
        t_end = time()
    elif m == 'treeswift_nexml':
        tree = read_tree_nexml(nexml_tree)
        t_start = time()
        for node in tree.traverse_preorder():
            pass
        t_end = time()
    elif m == 'treeswift_nexus':
        tree = read_tree_nexus(nexus_tree)
        t_start = time()
        for node in tree.traverse_preorder():
            pass
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# root distance order traversal
def rootdistorder(m):
    if m == 'treeswift_newick':
        tree = read_tree_newick(treestr)
        t_start = time()
        for node in tree.traverse_rootdistorder():
            pass
        t_end = time()
    elif m == 'treeswift_nexml':
        tree = read_tree_nexml(nexml_tree)
        t_start = time()
        for node in tree.traverse_rootdistorder():
            pass
        t_end = time()
    elif m == 'treeswift_nexus':
        tree = read_tree_nexus(nexus_tree)
        t_start = time()
        for node in tree.traverse_rootdistorder():
            pass
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# total branch length
def total_branch_length(m):
    if m == 'treeswift_newick':
        tree = read_tree_newick(treestr)
        t_start = time()
        tree.edge_length_sum()
        t_end = time()
    elif m == 'treeswift_nexml':
        tree = read_tree_nexml(nexml_tree)
        t_start = time()
        tree.edge_length_sum()
        t_end = time()
    elif m == 'treeswift_nexus':
        tree = read_tree_nexus(nexus_tree)
        t_start = time()
        tree.edge_length_sum()
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# load tree
def load_tree(m):
    if m == 'treeswift_newick':
        t_start = time()
        tree = read_tree_newick(treestr)
        t_end = time()
    elif m == 'treeswift_nexml':
        t_start = time()
        tree = read_tree_nexml(nexml_tree)
        t_end = time()
    elif m == 'treeswift_nexus':
        t_start = time()
        tree = read_tree_nexus(nexus_tree)
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# memory
def measure_memory(m):
    if m == 'treeswift_newick':
        m_start = memory()
        t = read_tree_newick(treestr)
        m_end = memory()
    elif m == 'treeswift_nexml':
        m_start = memory()
        t = read_tree_nexml(nexml_tree)
        m_end = memory()
    elif m == 'treeswift_nexus':
        m_start = memory()
        t = read_tree_nexus(nexus_tree)
        m_end = memory()
    else:
        assert False, "Invalid tool: %s"%m
    return m_end-m_start

TASKS = {
    'distance_matrix':distance_matrix,
    'inorder':inorder,
    'ladderize':ladderize,
    'levelorder':levelorder,
    'load_tree':load_tree,
    'memory':measure_memory,
    'mrca':mrca,
    'postorder':postorder,
    'preorder':preorder,
    'rootdistorder':rootdistorder,
    'total_branch_length':total_branch_length,
}

# run
if argv[3] not in TASKS:
    print("Invalid task. Valid options: %s"%', '.join(sorted(TASKS.keys()))); exit(1)
print(TASKS[argv[3]](argv[2]))
