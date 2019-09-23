#!/usr/bin/env python
#what's changed?  DETECT-2 no longer needs to be filtered.  

import os.path
import sys
from datetime import datetime as dt

def create_swissprot_map(SWISS_PROT_MAP):
    mapping_dict = {}
    with open(SWISS_PROT_MAP, "r") as mapping:
        for line in mapping.readlines():
            line_as_list = line.split("\t")
            mapping_dict[line_as_list[0]] = set(line_as_list[2:])
    return mapping_dict
    
#not used anymore. We don't filter DETECT-2
def filter_and_export_detect_predictions(detect_file, detect_ECs):
    with open(detect_file, "r") as topred:
        with open(detect_ECs, "w") as cutoff:
            for line in topred.readlines():
                line_as_list = line.split("\t")
                if line_as_list[2] == "probability":
                    continue
                if float(line_as_list[2]) >= 0.2 and int(line_as_list[3]) > 5:
                    cutoff.write(line)
                    
                    
def filter_and_export_priam_predictions(priam_file, priam_ECs):
    with open(priam_file, "r") as ECs:
        with open(priam_ECs, "w") as processedECs:
            ID = None
            EC = None
            for line in ECs.readlines():
                if line.startswith(">"):
                    ID = line.split(" ")[0][1:].strip("\n")
                    continue
                elif line == "\n":
                    continue
                elif ID and not line.startswith("#"):
                    EC = line.split(" ")[0]
                    processedECs.write("\t".join([ID, EC]))
                    processedECs.write("\n")
                    ID = None
                    EC = None
                    continue 


def filter_and_export_diamond_predictions(diamond_file, diamond_ECs):
    with open(diamond_file, "r") as blastout:
        with open(diamond_ECs, "w") as ecout:
            for line in blastout.readlines():
                line_as_list = line.strip().split("\t")
                for EC in mapping_dict:
                    if line_as_list[1] in mapping_dict[EC]:
                        ecout.write("\t".join([line_as_list[0], EC + "\n"]))


def import_detect_ec(detect_fbeta_file):
    #DETECT-2's fbeta file gets left as-is.
    gene_ec_dict = dict()
    with open(detect_fbeta_file, "r") as detect_fbeta:
        for line in detect_fbeta:
            list_line = line.split("\t")
            if(list_line[2] == "probability"):
                continue
            else:
                key = list_line[0]
                EC_val = list_line[1]
                if(key in gene_ec_dict):
                    gene_ec_dict[key].append(EC_val)
                else:
                    gene_ec_dict[key] = [EC_val]
    return gene_ec_dict
    
def import_priam_ec(priam_sequence_ec):
    gene_ec_dict = dict():
    query_name = "None"
    with open(priam_sequence_ec, "r") as priam_ec:
        for line in priam_ec:
            if(line.startswith(">")):
                query_name = line
            else:
                if not(line.startswith("#")):
                    list_line = line.split("\t")
                    ec = list_line[0]
                    probability = list_line[1]
            if(query_name is "None"):
                print(dt.today(), "This shouldn't be happening.  a line was skipped")
            
def import_diamond_ec(diamond_proteins_blastout, swissprot_map_dict):
    gene_ec_dict = dict()
    with open(diamond_proteins_blastout, "r") as diamond_ec:
        for item in diamond_ec:
            #note that this DIAMOND call has special parameters (for some reason)
            list_line = item.split("\t")
            query_name = list_line[0]
            swissprot_name = list_line[1]
            e_value = list_line[6]
            bitscore = list_line[7]
            identity_score = list_line[10]
            
        

if __name__ == "__main__":
    Input_File = sys.argv[1]
    Input_Name = os.path.splitext(os.path.basename(Input_File))[0]
    detect_file = sys.argv[2]
    detect_dir = os.path.dirname(detect_file)
    priam_file = sys.argv[3]
    priam_dir = os.path.dirname(priam_file)
    diamond_file = sys.argv[4]
    diamond_dir = os.path.dirname(diamond_file)
    SWISS_PROT = sys.argv[5]
    SWISS_PROT_MAP = sys.argv[6]
    Output_Dir = sys.argv[7]
    
    detect_ECs = os.path.join(detect_dir, Input_Name + ".toppred.cutoff")
    priam_ECs = os.path.join(priam_dir, Input_Name + ".ECs")
    diamond_ECs = os.path.join(diamond_dir, Input_Name + ".ECs")
    
    
    swissprot_map_dict = create_swissprot_map(SWISS_PROT_MAP)
    
    