import sys
import os
import pandas as pd

def patch_name(name):
    if(name == "UNMAPPED"):
        return "UNMAPPED|UNMAPPED"
    else:
        return name

def name_check(name):
    name_sublist = name.split("|")
    return len(name_sublist)


if __name__ == "__main__":
    gene_file = sys.argv[1]
    new_name = sys.argv[2]
    gene_df = pd.read_csv(gene_file, sep = "\t",  error_bad_lines = False)
    gene_df.columns = (['GeneID', 'RPK'])
    
    selected_df = gene_df#[gene_df["GeneID"].str.contains("g__")]
    selected_df["GeneID"] = selected_df["GeneID"].apply(lambda x: patch_name(x))
    selected_df["name_length"] = selected_df["GeneID"].apply(lambda x: name_check(x))
    
    selected_df = selected_df[selected_df["name_length"] == 2]
    selected_df.drop("name_length", axis = 1, inplace = True)
    print(selected_df)
    selected_df["taxa"] = selected_df["GeneID"].apply(lambda x:x.split("|")[1])
    new_taxa_name = os.path.join(new_name, new_name + "_humann2_taxa_found_0.csv")
    selected_df.to_csv(new_taxa_name, index = False)
    selected_df = selected_df.groupby(["taxa"], as_index = False).sum()
    total_rpk = selected_df["RPK"].sum()
    selected_df["percentage_rpk"] = 1
    selected_df["percentage_rpk"] = selected_df["percentage_rpk"].mask(selected_df["percentage_rpk"] > 0, 100 * selected_df["RPK"] / total_rpk)
    print(selected_df)
    new_taxa_name = os.path.join(new_name, new_name + "_humann2_taxa_found_1.csv")
    selected_df.to_csv(new_taxa_name, index = False)