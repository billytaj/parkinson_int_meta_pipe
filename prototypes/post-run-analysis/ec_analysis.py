import os
import sys

import pandas as pd

#this script pulls out the ECs, matches against each metabolic pathway, and gives us a breakdown of how much each pathway is covered
#it also lays out which tool did what, and how


#read in RPKM, and thepathway map
def make_result(ec_list, ec_path_df):
    #select all paths with the EC present
    frames_list = []
    for item in ec_list:
        selected_df = ec_path_df[ec_path_df["Ecs"].apply(lambda x: item in x)]
        selected_df["tally"] = 1
        frames_list.append(selected_df)    

    #then glue it together
    result = pd.concat(frames_list)
    
    #do a groupby so that the paths' tallies collect into a sum
    result = result.groupby("Pathway", as_index = False).sum()
    result["num_Ecs"] = result["num_Ecs"].mask(result["tally"] > 1, result["num_Ecs"] / result["tally"])
    #result["coverage"] = 1
    #result["coverage"] = result["coverage"].mask(result["tally"] > 0, result["tally"] *100 / result["num_Ecs"])
    #print(result)
    return result


def turn_to_list(x):
    new_list = x.split("|")
    return new_list
    
#These 2 translate the gene abundance files into ECs, using their uniref90 code.    
    
def translate_uniref(ec_dict, key):
    try:
        return ec_dict[key]
    except:
        return 0
        
def humann2_gene_to_ec(humann2_gene_file, ec_map_file):
    #ec_map_file = sys.argv[1]
    #humann2_gene_file = sys.argv[2]
    #out_path = sys.argv[3]
    
    ec_dict = dict()
    
    ec_file = open(ec_map_file, "r")
    count = 0
    for item in ec_file:
        count += 1
        cleaned_item = item.strip("\n")
        item_list = cleaned_item.split("\t")
        ec_number = item_list[0]
        item_list.pop(0)
        for obj in item_list:
            ec_dict[str(obj)] = ec_number
    #ec_dict["UNMAPPED"] = 0
        
    #for item in ec_dict:
    #    print(item, ":", ec_dict[item])
    
    gene_map = pd.read_csv(humann2_gene_file, sep = "\t", error_bad_lines = False)
    #print(gene_map)
    
    ID_df = gene_map["# Gene Family"].str.split("|", expand = True)[0].drop_duplicates()
    
    ID_df = pd.DataFrame(ID_df)
    ID_df.columns = (["uniref"])
    ID_df["ec"] = ID_df["uniref"].apply(lambda x: translate_uniref(ec_dict, x))
    #out_name = out_path + "_uniref_to_ec.csv"
    #ID_df = ID_df[ID_df.ec != 0]
    #ID_df.to_csv(out_name, mode = "w", index = False)
    #print(ID_df)
    return ID_df    
    
    
def divide_num_ec(x):
    if(x["humann2"] == 0 and x["both"] > 1 and x["mpro"] > 1):
        return True
    elif(x["humann2"] > 1 and x["both"] > 1 and x["mpro"] == 0):
        return True
    elif(x["humann2"] > 1 and x["both"] == 0 and x["mpro"] > 1):
        return True
    else:
        return False
        

if __name__ == "__main__":
    ec_path_file = sys.argv[1]
    ec_map_file = sys.argv[2]
    rpkm_file = sys.argv[3]
    humann2_gene_file = sys.argv[4]
    write_dir = sys.argv[5]
    name_prefix = sys.argv[6]
    

    #get our ec path map searchable
    ec_path_df = pd.read_csv(ec_path_file)
    ec_path_df["Ecs"] = ec_path_df["Ecs"].apply(lambda x: turn_to_list(x))
    #print(ec_path_df)
     
    num_ecs = ec_path_df["num_Ecs"]
    
    
    
    #get our annotated ECs
    rpkm_df = pd.read_csv(rpkm_file, sep = "\t")
    rpkm_df = rpkm_df[["GeneID", "EC#"]]
    rpkm_df = rpkm_df.drop_duplicates("EC#")
    rpkm_df = rpkm_df[rpkm_df["EC#"] != "0.0.0.0"]    
    
    rpkm_ec = rpkm_df["EC#"].tolist()
    rpkm_ec_df = rpkm_df["EC#"].to_frame()#pd.DataFrame(rpkm_df["EC#"])
    rpkm_ec_df["tally"] = 1
    rpkm_ec_df["EC#"] = rpkm_ec_df["EC#"].apply(lambda x: str(x))
    rpkm_ec_df.columns = ["ec", "tally"]
    rpkm_ec_df.reset_index(inplace = True)
    #print("RPKM DF")
    #print(rpkm_ec_df)
    
    #make the results
    mpro_ec_coverage_df = make_result(rpkm_ec, ec_path_df)
    #mpro_write_path = os.path.join(write_dir, "mpro_ec.csv")
    #mpro_ec_coverage_df.to_csv(mpro_write_path, index = False)
    
    #get humann2's list
    humann2_df = humann2_gene_to_ec(humann2_gene_file, ec_map_file)#pd.read_csv(humann2_file)
    humann2_df = humann2_df[humann2_df["ec"] != "0"]
    humann2_df = humann2_df.drop_duplicates("ec")
    humann2_ec = humann2_df["ec"].tolist()
    
    humann2_ec_coverage_df = make_result(humann2_ec, ec_path_df)
    #humann2_write_path = os.path.join(write_dir, "humann2_ec.csv")
    #humann2_ec_coverage_df.to_csv(humann2_write_path, index = False)
    
    humann2_ec_df = humann2_df["ec"].to_frame()#pd.DataFrame(humann2_df["ec"])
    humann2_ec_df["tally"] = 2
    humann2_ec_df["ec"] = humann2_ec_df["ec"].apply(lambda x: str(x)) #type correction so we can do merges and groupbys
    humann2_ec_df.reset_index(inplace = True)
    #print("HUMANN2 df")
    #print(humann2_ec_df)
    
    #get the intersection
    common_ec = list(set(rpkm_ec) & set(humann2_ec))
    
    #print("RPKM")
    #print(sorted(set(rpkm_ec)))
    #print("HUMANN2")
    #print(sorted(set(humann2_ec)))
    #print("COMMON")
    #print(sorted(set(common_ec)))
    #print("common ECs")
    #for item in common_ec:
    #    print(item)
        
    common_ec_coverage_df = make_result(common_ec, ec_path_df)
    #common_ec_coverage_df["who"] = "both"
    common_ec_coverage_df["both"] = common_ec_coverage_df["tally"]
    common_ec_coverage_df.drop("tally", axis = 1, inplace = True)
    #common_write_path = os.path.join(write_dir, "common_ec.csv")
    #common_ec_coverage_df.to_csv(common_write_path, index = False)
    
    
    #rpkm_ec_df and humann2_ec_df are going to contain ECs and a tally number for groupby
    
    mpro_only_ec = list(set(rpkm_ec).difference(set(humann2_ec)))
    #print("MPRO ONLY")
    #print(sorted(set(mpro_only_ec)))
    
    mpro_only_ec_coverage_df = make_result(mpro_only_ec, ec_path_df)
    #mpro_only_ec_coverage_df["who"] = "mpro"
    mpro_only_ec_coverage_df["mpro"] = mpro_only_ec_coverage_df["tally"]
    #print(mpro_ec_coverage_df)
    mpro_only_ec_coverage_df.drop("tally", axis = 1, inplace = True)
    #mpro_only_write_path = os.path.join(write_dir, "mpro_only.csv")
    #mpro_only_ec_coverage_df.to_csv(mpro_only_write_path, index = False)
    
    humann2_only_ec = list(set(humann2_ec).difference(set(rpkm_ec)))
    #print("HUMANN2 only")
    #print(sorted(set(humann2_only_ec)))
    humann2_only_coverage_df = make_result(humann2_only_ec, ec_path_df)
    #humann2_only_coverage_df["who"] = "humann2"
    humann2_only_coverage_df["humann2"] = humann2_only_coverage_df["tally"]
    humann2_only_coverage_df.drop(["tally"], axis = 1, inplace = True)
    print(humann2_only_coverage_df)
    #humann2_only_write_path = os.path.join(write_dir, "humann2_only.csv")
    #humann2_only_coverage_df.to_csv(humann2_only_write_path, index = False)
    #next, we need a EC -> pathway 1:1 list. aka: a dict
    #no, this can't happen.  some ecs exist in many paths
    
    #results_df = common_ec_coverage_df
    #results_df["mpro"] = mpro_only_ec_coverage_df["tally"]
    #results_df["humann2"] = humann2_ec_coverage_df["tally"]
    
    results_df = pd.concat([common_ec_coverage_df, mpro_only_ec_coverage_df, humann2_only_coverage_df])
    #num_Ecs = pd.DataFrame(results_df[["num_Ecs","Pathway"]])
    #num_Ecs = num_Ecs[~num_Ecs.index.duplicated()]
    print("--------------------------------------")
    #print(num_Ecs)
    
    
    results_df.sort_values(by = ["Pathway"], inplace = True)
    avg_results_df = results_df.groupby("Pathway", as_index = False).mean()
    print("AVERAGE")
    print(avg_results_df)
    results_df = results_df.groupby("Pathway", as_index = False).sum()
    
    #results_df["num_Ecs"] = results_df["num_Ecs"].mask(results_df["num_Ecs"] > 0, results_df["num_Ecs"] / ( (results_df["mpro"]/results_df["mpro"]) + (results_df["humann2"] / results_df["humann2"]) + (results_df["both"] / results_df["both"])))
    
    #results_df["division"] = (results_df["mpro"]/results_df["mpro"]) + (results_df["humann2"] / results_df["humann2"]) + (results_df["both"] / results_df["both"])
    results_df["num_Ecs"] = avg_results_df["num_Ecs"]
    results_df["sample"] = name_prefix
    new_name = os.path.join(write_dir, name_prefix + "_merged_results.csv")
    results_df.to_csv(new_name, index = False)
    print("--------------------------------------")
    print("final")
    print(results_df)
    