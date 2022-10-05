'''
Comprehensive construction of empCpds via generic tree structures.
Each empCpd = annoTree.
1. construct all isotopologues into a list of trees
2. attach in-src modifications (adducts) to nodes in each tree
3. for peaks without isotopologues, add to new trees and attach adducts

This applies to regular LC-MS data, but also
enables easy analysis of isotope tracing and chemical labeling data.

Example data:

    >>> node2tree = construct_isotopic_trees(plist)
    Found 873 isotopic pairs, 703 trees and 170 in branches.
    >>> node2tree['F795'].root
    'F649'
    >>> node2tree['F795'].nodes
    {'F649': Node(tag=179.1123@172.9, identifier=F649, data=None), 'F755': Node(tag=180.1156@172.9, identifier=F755, data=None), 
    'F794': Node(tag=181.1198@171.7, identifier=F794, data=None), 'F795': Node(tag=181.1198@173.3, identifier=F795, data=None), 
    'F839': Node(tag=182.1233@171.7, identifier=F839, data=None)}
    >>> node2tree['F795'].show()
    179.1123@172.9
    └── 180.1156@172.9
        ├── 181.1198@171.7
        │   └── 182.1233@171.7
        └── 181.1198@173.3


Trees use peak IDs, but full peak data are in json peak lists.

Shuzhao 2022-09-22

#
# to-dos: Above only includes isotopic peaks; will add peaks without C13 counterparts into empCpd trees
#
# 

The trees can be annotated by customized methods, e.g.
- in-house compound library
- targeted pathway search
- public databases
- linking tandem mass spec data


'''

import os
import json
# from itertools import combinations
import treelib

from .search import build_centurion_tree, find_all_matches_centurion_indexed_list


def make_peak_tag(peak):
    '''
    peak format: {'id': 'F1', 'mz': 60.0808, 'rtime': 117.7, 'intensities': [250346.0], 'representative_intensity': 250346.0}
    '''
    return str(round(peak['mz'], 4)) + '@' + str(round(peak['rtime'], 1))

def make_peak_dict(peak_list):
    peak_dict = {}
    for p in peak_list:
        peak_dict[p['id']] = p
    return peak_dict

def construct_isotopic_trees(peak_list, 
                    search_patterns = [(1.003355, '13C/12C', (0, 0.8))], 
                    mz_tolerance_ppm=5, 
                    isotope_rt_tolerance=2, 
                    check_isotope_ratio = False,
                    tree_depth_limit=10):
    '''
    Make isotopic trees from peak_list. 
    tree_depth_limit: limit of tree depth, i.e. steps on a branch.
    For other parameters, see get_isotopic_pairs.
    Return node2tree, a dict of peak to tree mapping.
    No peak is assigned to more than one trees. But some close peaks can't be distinguished here.

    Example
    =======
    >>> node2tree['F515'].show()
    173.0922@172.9
    └── 174.0956@171.5
        ├── 175.0991@172.4
        │   ├── 176.1031@171.7
        │   │   ├── 177.1065@171.7
        │   │   └── 177.1065@174.1
        │   └── 176.1031@174.1
        └── 175.0991@173.1

    >>> node2tree['F796'].show()
    180.1244@171.2
    └── 181.1279@170.1
        └── 182.1313@170.1
    '''
    peak_dict = make_peak_dict(peak_list)
    mztree = build_centurion_tree(peak_list)
    isotopologues = get_isotopic_pairs(peak_list, mztree, search_patterns, mz_tolerance_ppm, 
                    isotope_rt_tolerance, check_isotope_ratio)
    # isotopologues format: [ ((195, 'anchor'), (206, '13C/12C')), ...]
    # build trees
    annoTrees, branches = [], []
    all_target_nodes = set([x[1][0] for x in isotopologues])
    for pair in isotopologues:
        if pair[0][0] not in all_target_nodes:   # if source_node appears in any of the target_nodes, it's not a root
            tree = treelib.Tree()
            tree.create_node(make_peak_tag(peak_dict[pair[0][0]]), pair[0][0], data=pair[0][1])
            tree.create_node(make_peak_tag(peak_dict[pair[1][0]]), pair[1][0], parent=pair[0][0], data=pair[1][1])
            annoTrees.append(tree)
        else:
            branches.append(pair)
    
    print("Found %d isotopic pairs, %d trees and %d in branches." %(len(isotopologues), len(annoTrees), len(branches)))
    
    node2tree = {}
    for tree in annoTrees:
        for node in tree.nodes:
            node2tree[node] = tree

    # do branches now
    remaining = []
    for pair in branches:          
        if pair[0][0] in node2tree:                  # pair[0] already in a tree
            try:
                this_tree = node2tree[pair[0][0]]
                this_tree.create_node( make_peak_tag(peak_dict[pair[1][0]]), pair[1][0], parent=pair[0][0], data=pair[1][1] )
                node2tree[pair[1][0]] = this_tree
            except treelib.exceptions.DuplicatedNodeIdError:
                # print("already included ", pair)
                pass                 # pair already in a tree
        else:
            remaining.append(pair)

    steps = 0
    while remaining and steps < tree_depth_limit:
        tmp = []
        for pair in remaining:          
            if pair[0][0] in node2tree: 
                try:
                    this_tree = node2tree[pair[0][0]]
                    this_tree.create_node( make_peak_tag(peak_dict[pair[1][0]]), pair[1][0], parent=pair[0][0], data=pair[1][1] )
                    node2tree[pair[1][0]] = this_tree
                except treelib.exceptions.DuplicatedNodeIdError:
                    # print("already included ", pair)
                    pass               # pair already in a tree
            else:
                tmp.append(pair)
        remaining = tmp
        steps += 1

    return node2tree

def get_isotopic_pairs(list_peaks, 
                    mztree, 
                    search_patterns = [(1.003355, '13C/12C', (0, 0.8))],
                    mz_tolerance_ppm=5, 
                    isotope_rt_tolerance=2, 
                    check_isotope_ratio = False,
                    ):
    '''
    To find all isotope pairs. Similar to search.get_seed_empCpd_signatures, but return unidirectional pairs only.
    If input peaks have overlaps/duplicates, result will contain the redundant overlap peaks.

    Input
    =====
    list_peaks: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'rtime': 654, 'height': 14388.0, 'id': 555}, ...]
    mztree: indexed list_peaks
    mz_tolerance_ppm: ppm tolerance in examining m/z patterns.
    seed_search_patterns: a list in the format of [(mz difference, notion, (ratio low limit, ratio high limit)), ..]
            This can be obtained through search.isotopic_patterns. The ratios are optional, because 
            1) naturally occuring constrains are based on chemical formula;
            2) rules are different when isotope tracers are introduced to the experiments.
    isotope_rt_tolerance: tolerance threshold for deviation in retetion time, arbitrary unit depending on input data.

    Return
    ======
    list of lists of peak pairs that match search_patterns patterns, e.g.
    [ ((195, 'anchor'), (206, '13C/12C')), ...]. 
    '''
    signatures = []
    for P1 in list_peaks:
        matched = [  ] 
        for _pair in search_patterns:
            (mass_difference, relation) = _pair[:2]
            tmp = find_all_matches_centurion_indexed_list(P1['mz'] + mass_difference, mztree, mz_tolerance_ppm)
            for P2 in tmp:
                if abs(P1['rtime']-P2['rtime']) <= isotope_rt_tolerance:
                    if check_isotope_ratio and len(_pair) > 2:  # checking abundance ratio
                        (abundance_ratio_min, abundance_ratio_max) = _pair[2]
                        if abundance_ratio_min*P1['height'] < P2['height'] < abundance_ratio_max*P1['height']:
                            matched.append( ((P1['id'], 'anchor'), (P2['id'], relation)) )
                    else:
                        matched.append( ((P1['id'], 'anchor'), (P2['id'], relation)) )
 
        signatures += matched

    return signatures


def rt_matched_by_tolerance(P1, P2, rt_tolerance):
    return abs(P2['rtime']-P1['rtime']) < rt_tolerance

def rt_compared_by_values(P1, P2, rt_tolerance=None):
    return P2['rtime'] > P1['rtime']

def merge_trees_by_modifications(trees, list_peaks, search_patterns, rt_verify_function,
                    mz_tolerance_ppm=5, rt_tolerance=10):
    '''
    Merge trees that are from same compound but with  modifications. 
    rt_verify_function applies rules on retention time.
    '''
    peak_dict = make_peak_dict(list_peaks)
    tree_dict = {}
    for tree in trees: tree_dict[tree.root] = tree
    print("Merging adducts on %d trees..." %len(trees))
    matched_pairs = []
    relevant_peak_centTree = build_centurion_tree( [peak_dict[p] for p in tree_dict.keys()] )
    for p in tree_dict:
        P1 = peak_dict[p]
        matched = [  ] 
        for _pair in search_patterns:
            (mass_difference, relation) = _pair[:2]
            tmp = find_all_matches_centurion_indexed_list(P1['mz'] + mass_difference, 
                                                    relevant_peak_centTree, mz_tolerance_ppm)
            for P2 in tmp:
                if rt_verify_function(P1, P2, rt_tolerance):
                    matched.append( (P1['id'], P2['id'], relation) )
        matched_pairs += matched

    list_primary = set([x[0] for x in matched_pairs])
    list_adducts = set([x[1] for x in matched_pairs])
    overlap = list_primary.intersection(list_adducts)
    union = list_primary.union(list_adducts)
    if overlap: print("Unresolved multiple relationships: ", overlap)

    good_trees = [tree_dict[x] for x in tree_dict if x not in union]
    for P in matched_pairs:
        p1, p2 = P[:2]
        tree_dict[p2].get_node(p2).data = P[2]      # add relation note
        try:
            tree_dict[p1].paste(p1, tree_dict[p2])
        except ValueError:                  # likely Duplicated nodes
            # print(tree_dict[p1].show(), tree_dict[p2].show())
            for node in tree_dict[p2].nodes:
                if node not in tree_dict[p1].nodes:
                    tree_dict[p1].add_node(tree_dict[p2].nodes[node], parent=tree_dict[p1].root)
        good_trees.append(tree_dict[p1])

    print("Got %d merged trees." %len(good_trees))
    return good_trees


def merge_trees_by_insrc_modifications(trees, list_peaks,
                    search_patterns = [(1.0078, 'H'), (21.9820, 'Na/H'), (41.026549, 'Acetonitrile')], 
                    mz_tolerance_ppm=5, rt_tolerance=10):
    '''
    Merge trees that are from same compound but with diff adducts or neutral loss,
    with user supplied search_patterns (can be from search.common_adducts).
    These in-source modifications must have same retention time.
    '''
    return merge_trees_by_modifications(trees, list_peaks,
                    search_patterns, rt_verify_function=rt_matched_by_tolerance,
                    mz_tolerance_ppm=mz_tolerance_ppm, rt_tolerance=rt_tolerance)

def merge_trees_by_derivatization(trees, list_peaks,
                    search_patterns = [(161.08407, 'DmPA'), (233.05105, 'Dens'), (247.07793, 'DnsHz')], 
                    mz_tolerance_ppm=5, ):
    '''
    Merge trees that are from same compound but with chemical labeling/derivatization,
    which leads to greater retention time.
    '''
    return merge_trees_by_modifications(trees, list_peaks,
                    search_patterns, rt_verify_function=rt_compared_by_values,
                    mz_tolerance_ppm=mz_tolerance_ppm, rt_tolerance=None)


def export_json_trees(trees, outfile="export_annoTree.tsv"):
    pass


def export_tsv_trees(trees, outfile="export_annoTree.tsv"):
    s = 'Feature_ID\tFeature_tag\troot\troot_tag\trelation\n'
    for tree in trees:
        for node in tree.nodes:
            ND, ROOT = tree.get_node(node), tree.get_node(tree.root)
            s += '\t'.join([ND.identifier, ND.tag, ROOT.identifier, ROOT.tag, ND.data or '']) + '\n'
    with open(outfile, 'w') as O:
        O.write(s)




# ------------------------------------------------------------

class empTree:
    '''
    Not used now. 
    Container for empirical tree, as ambiguous branches don't fit a strict tree structure.  

    Trees break down here when leaves are shared:
        118.0652@109.9
        ├── 119.0685@109.2
        └── 159.0917@109.9
            └── 160.0952@110.4

        159.0917@110.4
        └── 160.0952@110.4

    '''
    def __init__(self) -> None:
        self.ambiguous = False
        self.uniqe_nodes = []
        self.trees = []


