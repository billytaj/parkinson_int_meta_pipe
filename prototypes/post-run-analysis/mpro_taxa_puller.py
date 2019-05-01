import sys
import os
import pandas as pd


def make_genus_label (taxa):
    taxa_tree = taxa.split(";")
    if(taxa == "unclassified"):
        return "unclassified"
    elif(len(taxa_tree) < 7):
        return "g_unknown"
    else:
        if(taxa_tree[6] == "g_"):
            return "g_unknown"
        else:
            return taxa_tree[6]



if __name__ == "__main__":
    gene_file = sys.argv[1]
    new_name = sys.argv[2]
    gene_df = pd.read_csv(gene_file, sep = "\t",  error_bad_lines = False)
    
    selected_df = gene_df#gene_df[~gene_df["Taxonomy"].str.contains("unclassified")]
    selected_df["RPK"] = 1
    selected_df["RPK"] = selected_df["RPK"].mask(selected_df["RPK"] > 0, selected_df["Reads"] / (selected_df["Length"] / 1000))
    selected_df.sort_values("RPK", inplace = True, ascending = False)
    new_taxa_name = os.path.join(new_name, new_name + "_mpro_taxa_found_0.csv")
    selected_df.to_csv(new_taxa_name, index = False)
    selected_df = selected_df.groupby(["Taxonomy"], as_index = False).sum()
    selected_df.sort_values("RPK", inplace = True, ascending = False)
    print(selected_df)
    new_taxa_name = os.path.join(new_name, new_name + "_mpro_taxa_found_1.csv")
    selected_df.to_csv(new_taxa_name, index = False)
    
    
    #selected_df = selected_df.groupby(["Taxonomy"], as_index = False).sum()
    
    selected_df["genus"] = selected_df["Taxonomy"].apply(lambda x: make_genus_label(x))
    print(selected_df)
    
    new_taxa_name = os.path.join(new_name, new_name + "_mpro_taxa_found_2.csv")
    selected_df = selected_df.groupby(["genus"], as_index = False).sum()
    selected_df["percentage_rpk"] = 1
    total_rpk = selected_df["RPK"].sum()
    selected_df["percentage_rpk"] = selected_df["percentage_rpk"].mask(selected_df["percentage_rpk"] > 0, 100 * selected_df["RPK"] / total_rpk)
    selected_df.sort_values("percentage_rpk", ascending = False, inplace = True)
    
    selected_df.to_csv(new_taxa_name, index = False)