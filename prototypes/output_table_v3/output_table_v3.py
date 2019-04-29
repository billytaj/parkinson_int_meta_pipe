#This is output_table_v3.  with rationales and possibly cleaner ways to do things.
import sys
import matplotlib
from matplotlib import cm
import os

if __name__ == "__main__":
    gene2read = sys.argv[1]
    gene2read_dict = dict()
    infile = open(gene2read, "r")
    for line in infile:
        cols = line.split("\t")
        gene = cols[0]
        gene_len = cols[1]
        reads = []
        for read in cols[3:]:
            reads.append(read.strip("\n"))
        #mapped_reads += len(reads)
        if gene in gene2read_dict:
            gene2read_dict[gene][1].extend(reads)
        else:
            gene2read_dict[gene] = (gene_len, reads)

    i = 0
    for item in gene2read_dict:
        print("KEY:", item)
        print(gene2read_dict[item])
        i += 1
        if i > 10:
            break
    """    
    # Two methods to define taxa in order of increasing priority: Cutoff or Tax ID list
    cutoff = sys.argv[1]                    #IN: Proportion of annotated reads -> 0.01 from the example command list.
    ID_list = sys.argv[2]                   #IN: literally a list of taxid we can optionally use.  If nothing, it's empty and it's fine.
    nodes = sys.argv[3]                     #IN: nodes.dmp
    names = sys.argv[4]                     #IN: names.dmp
    gene2read = sys.argv[5]                 #IN: genes -> reads (from GA)
    read2taxid = sys.argv[6]                #IN: reads -> taxid (from TA)
    gene2EC = sys.argv[7]                   #IN: genes -> EC mapping (from EC, EC.All)
    show_unclassified_flag = sys.argv[8]    #IN: either true or false.  Pull from config.  default: true
    raw_count = sys.argv[9]                 #OUT
    RPKM = sys.argv[10]                     #OUT
    cytoscape = sys.argv[11]                #OUT
    """
    