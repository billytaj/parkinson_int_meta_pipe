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
    result["coverage"] = 1
    result["coverage"] = result["coverage"].mask(result["tally"] > 0, result["tally"] *100 / result["num_Ecs"])
    #print(result)
    return result


def turn_to_list(x):
    new_list = x.split("|")
    return new_list
    

if __name__ == "__main__":
    ec_path_file = sys.argv[1]
    rpkm_file = sys.argv[2]
    humann2_file = sys.argv[3]
    write_dir = sys.argv[4]
    name_prefix = sys.argv[5]
    

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
    humann2_df = pd.read_csv(humann2_file)
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
    common_ec_coverage_df["who"] = "both"
    #common_write_path = os.path.join(write_dir, "common_ec.csv")
    #common_ec_coverage_df.to_csv(common_write_path, index = False)
    
    
    #rpkm_ec_df and humann2_ec_df are going to contain ECs and a tally number for groupby
    
    mpro_only_ec = list(set(rpkm_ec).difference(set(humann2_ec)))
    #print("MPRO ONLY")
    #print(sorted(set(mpro_only_ec)))
    
    mpro_only_ec_coverage_df = make_result(mpro_only_ec, ec_path_df)
    mpro_only_ec_coverage_df["who"] = "mpro"
    #mpro_only_write_path = os.path.join(write_dir, "mpro_only.csv")
    #mpro_only_ec_coverage_df.to_csv(mpro_only_write_path, index = False)
    
    humann2_only_ec = list(set(humann2_ec).difference(set(rpkm_ec)))
    #print("HUMANN2 only")
    #print(sorted(set(humann2_only_ec)))
    humann2_only_coverage_df = make_result(humann2_only_ec, ec_path_df)
    humann2_only_coverage_df["who"] = "humann2"
    #humann2_only_write_path = os.path.join(write_dir, "humann2_only.csv")
    #humann2_only_coverage_df.to_csv(humann2_only_write_path, index = False)
    #next, we need a EC -> pathway 1:1 list. aka: a dict
    #no, this can't happen.  some ecs exist in many paths
    
    
    results_df = pd.concat([common_ec_coverage_df, mpro_only_ec_coverage_df, humann2_only_coverage_df])
    results_df.sort_values(by = ["Pathway"], inplace = True)
    new_name = os.path.join(write_dir, name_prefix + "_merged_results.csv")
    results_df.to_csv(new_name, index = False)