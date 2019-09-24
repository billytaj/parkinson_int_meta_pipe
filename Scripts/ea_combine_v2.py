#!/usr/bin/env python
#what's changed?  DETECT-2 no longer needs to be filtered.  

import os.path
import sys
from datetime import datetime as dt
import multiprocessing as mp

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


def import_detect_ec(detect_fbeta_file, gene_ec_dict):
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
    #return gene_ec_dict
    
def import_priam_ec(priam_sequence_ec, gene_ec_dict):
    gene_ec_dict = dict():
    query_name = "None"
    ec_list = []
    with open(priam_sequence_ec, "r") as priam_ec:
        for line in priam_ec:
            if(line.startswith(">")):
                #new line
                if(len(ec_list) > 0):
                    if(query_name is "None"):
                        print(dt.today(), "This shouldn't happen.  full EC list.  no query")
                    else:
                        gene_ec_dict[query_name] = ec_list
                query_name = line.strip(">")
            else:
                if not(line.startswith("#")):
                    list_line = line.split("\t")
                    ec = list_line[0]
                    probability = float(list_line[1])
                    if(probability >= 0.5):
                        ec_list.append(ec)
                        
            if(query_name is "None"):
                print(dt.today(), "This shouldn't be happening.  a line was skipped")
    #return gene_ec_dict
    
def import_diamond_ec(diamond_proteins_blastout, swissprot_map_dict, gene_ec_dict):
    gene_ec_dict = dict()
    with open(diamond_proteins_blastout, "r") as diamond_ec:
        for item in diamond_ec:
            #note that this DIAMOND call has special parameters (for some reason)
            list_line = item.split("\t")
            query_name = list_line[0]
            swissprot_name = list_line[1]
            query_start = int(list_line[2])
            query_end = int(list_line[3])
            query_length = query_end - query_start
            bitscore = int(list_line[7])
            coverage_percent = float(list_line[8])
            identity_percent = float(list_line[10])
            
            if(query_length >= 100):
                if(bitscore >= 60):
                    #take this hit
                    EC_value = swissprot_map_dict[swissprot_name]
                    if(query_name in gene_ec_dict):
                        gene_ec_dict[query_name].append(EC_val)
                    else:
                        gene_ec_dict[query_name] = EC_val
            else:
                if((identity_percent >= 85) and (coverage_percent >= 65)):
                    #take this hit
                    EC_value = swissprot_map_dict[swissprot_name]
                    if(query_name in gene_ec_dict):
                        gene_ec_dict[query_name].append(EC_val)
                    else:
                        gene_ec_dict[query_name] = EC_val
        #return gene_ec_dict
                    
            
            
            
            
        

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
    
    ec_process_list = []
    #-----------------------------------------
    # import the EC annotations
    manager = mp.Manager()
    diamond_ec_manager_dict = manager.dict()
    swissprot_map_dict = create_swissprot_map(SWISS_PROT_MAP)
    diamond_ec_process = mp.Process(
        target = import_diamond_ec, 
        args = (diamond_file, swissprot_map_dict, diamond_ec_manager_dict)
    )
    diamond_ec_process.start()
    ec_process_list.append(diamond_ec_process)
    
    priam_ec_manager_dict = manager.dict()
    priam_ec_process = mp.Process(
        target = import_priam_ec, 
        args = (priam_file, priam_ec_manager_dict)
    )
    priam_ec_process.start()
    ec_process_list.append(priam_ec_process)
    
    detect_ec_manager_dict = manager.dict()
    detect_ec_process = mp.Process(
        target = import_detect_ec,
        args = (detect_file, detect_ec_manager_dict)
    )
    detect_ec_process.start()
    ec_process_list.append(detect_ec_process)
    
    for item in ec_process_list:
        item.join()
    ec_process_list[:] = []
    
    
    diamond_ec_dict = dict(diamond_ec_manager_dict)
    priam_ec_dict = dict(priam_ec_manager_dict)
    detect_ec_dict = dict(detect_ec_manager_dict)
    
    #--------------------------------------------------
    # merge the results
    diamond_keys = set(diamond_ec_dict.keys())
    priam_keys = set(priam_ec_dict.keys())
    common_keys = diamond_keys & priam_keys
    
    common_dict = dict()
    for item in common_keys:
        inner_ec_list = []
        for priam_ec in priam_ec_dict[item]:
            inner_ec_list.append(priam_ec)
        for diamond_ec in diamond_ec_dict[item]:
            inner_ec_list.append(diamond_ec)
        common_dict[item]= inner_ec_list
    
    for key in detect_ec_dict.keys():
        if(key in common_dict):
            for detect_ec in detect_ec_dict[key]:
                common_dict[key].append(detect_ec)
        else:
            common_dict[key] = detect_ec_dict[key]
            
    common_dict = sorted(common_dict)        
    #----------------------------------------------
    #export the final ec list
    with open(os.path.join(Output_Dir, Input_Name + ".ECs_All"), "w") as ec_out:
    ec_out.writelines(sorted(All_preds))
    #export the final ec list
        
    
    
    
    
    