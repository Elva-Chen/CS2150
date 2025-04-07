#!/usr/bin/env python3
from Bio import Phylo
from itertools import combinations
from treeswift import read_tree_newick
import dendropy
import ete3
import numpy
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
from io import StringIO
from sys import argv
from time import time
if len(argv) != 4 or argv[2] not in {'treeswift','dendropy','biophylo','ete3'}:
    print("USAGE: %s <tree_file> <treeswift_or_dendropy_or_biophylo_or_ete3> <task>"%argv[0]); exit(1)
if argv[1].lower().endswith('.gz'):
    from gzip import open as gopen
    treestr = gopen(argv[1]).read().decode().strip()
else:
    treestr = open(argv[1]).read().strip()
treeio = StringIO(treestr) # for Bio.Phylo

# ladderize
def ladderize(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        tree.ladderize()
        t_end = time()
    elif m == 'biophylo':
        tree = Phylo.read(treeio, 'newick')
        t_start = time()
        tree.ladderize()
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        tree.ladderize()
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr,format=1)
        t_start = time()
        tree.ladderize()
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# distance matrix
def distance_matrix(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        tree.phylogenetic_distance_matrix()
        t_end = time()
    elif m == 'biophylo':
        # from http://biopython.org/wiki/Phylo_cookbook
        tree = Phylo.read(treeio, 'newick')
        t_start = time()
        allclades = list(tree.find_clades(order='level'))
        lookup = {}
        for i, elem in enumerate(allclades):
            lookup[elem] = i
        distmat = numpy.repeat(numpy.inf, len(allclades)**2)
        distmat.shape = (len(allclades), len(allclades))
        for parent in tree.find_clades(terminal=False, order='level'):
            for child in parent.clades:
                if child.branch_length:
                    distmat[lookup[parent], lookup[child]] = child.branch_length
        if not tree.rooted:
            distmat += distmat.transpose()
        numpy.matrix(distmat)
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        tree.distance_matrix()
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr,format=1)
        t_start = time()
        distance_matrix_ete(tree)
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# inorder traversal
def inorder(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        for node in tree.inorder_node_iter():
            pass
        t_end = time()
    elif m == 'biophylo':
        return NA
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        for node in tree.traverse_inorder():
            pass
        t_end = time()
    elif m == 'ete3':
        return NA
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# level-order traversal
def levelorder(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        for node in tree.levelorder_node_iter():
            pass
        t_end = time()
    elif m == 'biophylo':
        tree = Phylo.read(treeio, 'newick')
        t_start = time()
        for node in tree.find_clades(order='level'):
            pass
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        for node in tree.traverse_levelorder():
            pass
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr,format=1)
        t_start = time()
        for node in tree.traverse(strategy='levelorder'):
            pass
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# MRCA
def mrca(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        leaves = {l.taxon for l in tree.leaf_node_iter()}
        tree.mrca(taxa=leaves)
        t_end = time()
    elif m == 'biophylo':
        tree = Phylo.read(treeio, 'newick')
        t_start = time()
        leaves = tree.get_terminals()
        tree.common_ancestor(leaves)
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        leaves = {str(l) for l in tree.traverse_leaves()}
        tree.mrca(leaves)
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr,format=1)
        t_start = time()
        leaves = tree.get_leaf_names()
        tree.get_common_ancestor(leaves)
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# postorder traversal
def postorder(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        for node in tree.postorder_node_iter():
            pass
        t_end = time()
    elif m == 'biophylo':
        tree = Phylo.read(treeio, 'newick')
        t_start = time()
        for node in tree.find_clades(order='postorder'):
            pass
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        for node in tree.traverse_postorder():
            pass
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr,format=1)
        t_start = time()
        for node in tree.traverse(strategy='postorder'):
            pass
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# preorder traversal
def preorder(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        for node in tree.preorder_node_iter():
            pass
        t_end = time()
    elif m == 'biophylo':
        tree = Phylo.read(treeio, 'newick')
        t_start = time()
        for node in tree.find_clades(order='preorder'):
            pass
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        for node in tree.traverse_preorder():
            pass
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr,format=1)
        t_start = time()
        for node in tree.traverse(strategy='preorder'):
            pass
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# root distance order traversal
def rootdistorder(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        tree.calc_node_ages(is_force_max_age=True)
        for node in tree.ageorder_node_iter(descending=True):
            pass
        t_end = time()
    elif m == 'biophylo':
        return NA
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        for node in tree.traverse_rootdistorder():
            pass
        t_end = time()
    elif m == 'ete3':
        return NA
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# total branch length
def total_branch_length(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        tree.length()
        t_end = time()
    elif m == 'biophylo':
        tree = Phylo.read(treeio, 'newick')
        t_start = time()
        tree.total_branch_length()
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        tree.edge_length_sum()
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr,format=1)
        t_start = time()
        sum(node.dist for node in tree.traverse(strategy='preorder'))
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# load tree
def load_tree(m):
    if m == 'dendropy':
        t_start = time()
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_end = time()
    elif m == 'biophylo':
        t_start = time()
        tree = Phylo.read(treeio, 'newick')
        t_end = time()
    elif m == 'treeswift':
        t_start = time()
        tree = read_tree_newick(treestr)
        t_end = time()
    elif m == 'ete3':
        t_start = time()
        tree = ete3.Tree(treestr,format=1)
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# memory
def measure_memory(m):
    if m == 'dendropy':
        m_start = memory()
        t = dendropy.Tree.get(data=treestr, schema='newick')
        t.encode_bipartitions()
        m_end = memory()
    elif m == 'biophylo':
        m_start = memory()
        t = Phylo.read(treeio, 'newick')
        m_end = memory()
    elif m == 'treeswift':
        m_start = memory()
        t = read_tree_newick(treestr)
        m_end = memory()
    elif m == 'ete3':
        m_start = memory()
        t = ete3.Tree(treestr,format=1)
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
