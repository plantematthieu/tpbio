#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file, 'r') as f:
        lines = f.readlines()
        
    i = 1
    while i < len(lines):
        yield  lines[i][:-1]
        i += 4


def cut_kmer(read, kmer_size):
    i = 0
    while i <= (len(read) - kmer_size):
        yield read[i: kmer_size + i]
        i += 1



def build_kmer_dict(fastq_file, kmer_size):
    dic = dict()
    seq_generator = read_fastq(fastq_file)
    
    for seq in seq_generator:
        kmer_generator = cut_kmer(seq, kmer_size)
        for kmer in kmer_generator:
            try:
                dic[kmer] += 1
            except:
                dic[kmer] = 1
    return dic


def build_graph(kmer_dict):
    G = nx.DiGraph()
    
    for key in kmer_dict:
        G.add_node(key[:-1])
        G.add_node(key[1:])
        G.add_edge(key[:-1], key[1:], weight=kmer_dict[key])

    return G
    


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        
        if delete_entry_node:
            if delete_sink_node:
                graph.remove_nodes_from(path)
            else:
                graph.remove_nodes_from(path[:-1])
                
        else:
            if delete_sink_node:
                graph.remove_nodes_from(path[1:])
            else:
                graph.remove_nodes_from(path[1:-1])
            
    return graph


def std(data):
    return statistics.stdev(data)



def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    
    best_weight = max(weight_avg_list)
    best = [i for i, weight in enumerate(weight_avg_list) if weight == best_weight]
    
    if len(best) > 1:
        best2 = []
        max_length = max(path_length)
        for i in best:
            if path_length[i] == max_length:
                best2.append(i)

        path = random.choice(best2)
    
    else:
        path = best[0]
        
    sub = path_list[:path] + path_list[path +1:]
    return remove_paths(graph, sub, delete_entry_node, delete_sink_node)



def path_average_weight(graph, path):
    mean_list = []
    for i in range(len(path) -1):
        mean_list.append(graph[path[i]][path[i+1]]["weight"])
    return statistics.mean(mean_list)


def solve_bubble(graph, ancestor_node, descendant_node):
    paths = list(nx.shortest_simple_paths(graph, ancestor_node, descendant_node))
    path_length = [len(path) for path in paths]
    weight_avg_list = [path_average_weight(graph, path) for path in paths]
    
    return select_best_path(graph, paths, path_length, weight_avg_list)

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    l = []
    found = False
    for edge in graph.edges:
        found = False
        for edge1 in graph.edges:
            if edge[0] == edge1[1]:
                found = True
                break
        if found == False:
            l.append(edge[0])
    return l
    

def get_sink_nodes(graph):
    l = []
    found = False
    for edge in graph.edges:
        found = False
        for edge1 in graph.edges:
            if edge[1] == edge1[0]:
                found = True
                break
        if found == False:
            l.append(edge[1])
    return l

def get_contigs(graph, starting_nodes, ending_nodes):
    l = []
    for start in starting_nodes:
        for end in ending_nodes:
            for path in nx.all_simple_paths(graph, source=start, target=end):
                char = path[0]
                for c in path[1:]:
                    char += c[-1]
                l.append((char, len(char)))
                
    return l


def save_contigs(contigs_list, output_file):
    with open(output_file, 'a') as f:
        i = 0
        for contig in contigs_list:
            output = f'>contig_{i} len={contig[1]}\n{contig[0]}\n'
            f.write(fill(output))
            i += 1
    
def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    file = args.fastq_file
    size = args.kmer_size
    output = args.output_file
    
    graph = build_graph(build_kmer_dict(file, size))
    
    simp = simplify_bubbles(graph)
    
    starts = get_starting_nodes(simp)
    sinks = get_sink_nodes(simp)
    #graph = solve_entry_tips(res, starts)
    #graph = solve_out_tips(res, sinks)

    contigs = get_contigs(simp, get_starting_nodes(simp), get_sink_nodes(simp))
    save_contigs(contigs, output)

    
    
if __name__ == '__main__':
    main()


