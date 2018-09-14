import pandas as pd
import sys
import os

if __name__ == "__main__":
    ec_pathway_file = sys.argv[1]   #path->EC file  -> EC_pathway.txt
    rpkm_table_file = sys.argv[2]   #EC-> gene file -> RPKM_table.tsv
    pathway_superpathway_file = sys.argv[3] #path->superpath file -> pathway_to_superpathway.csv
    
    #imports the path->ec file
    path_df = pd.read_csv(ec_pathway_file, header=None, names=None, sep = '\t', skip_blank_lines = False)
    path_df.columns = (["path", "EC"])
    path_df["path"] = path_df["path"].apply(lambda x: x.split(":")[1])  
    path_df.index = path_df.groupby(["path"]).cumcount()
    temp_columns = path_df.index
    path_df = path_df.pivot(values='EC', columns = 'path')
    path_df = path_df.T
    path_df["EC_list"] = path_df[path_df.columns[:]].apply(lambda x: ",".join(x.dropna()), axis=1)
    path_df.drop(temp_columns, axis=1, inplace=True)
    
    rpkm_df = pd.read_csv(rpkm_table_file, sep = '\t', skip_blank_lines = False)
    
    #imports the path->superpath file
    super_df = pd.read_csv(pathway_superpathway_file, sep = ',', skip_blank_lines = False)
    super_df.drop(columns = ["Pathway"], inplace=True)
    
    
    super_df = super_df.join(path_df, how = 'left', on = 'Pathway ID')
    super_df.to_csv("this_sucks.csv", mode="w")
    super_df.drop(columns = ["Pathway ID"], inplace = True)
    
    super_df.index = super_df.groupby("Superpathway").cumcount()
    super_df = super_df.pivot(values = "EC_list", columns = 'Superpathway')
    
    super_df = super_df.T
    
    super_temp_columns = super_df.columns
    super_df["combined"] = super_df[super_df.columns[:]].apply(lambda x: ",".join(x.dropna()), axis=1)
    super_df.drop(super_temp_columns, axis=1, inplace=True)
    
    super_enzyme_df = pd.DataFrame(super_df.combined.str.split(',').tolist(), index = super_df.index)
    super_enzyme_df["count"] = super_enzyme_df.count(axis=1)
    #super_enzyme_df["Superpathway"] = super_enzyme_df.index
    #super_enzyme_df.index = range(super_enzyme_df.shape[0])
    super_enzyme_df.to_csv("superpath_enzyme_reg.csv", mode="w")
    superpath_list = list(super_enzyme_df.index)
    print(superpath_list)
    super_enzyme_df = super_enzyme_df.T
    for item in superpath_list:
        super_enzyme_df[item].fillna(9999)
        super_enzyme_df[item] = super_enzyme_df[item].unique()
        #print(item, super_enzyme_df[item].iloc[0])
    #print(list(super_enzyme_df.columns.values))
    super_enzyme_df.to_csv("superpath_enzyme.csv", mode="w")
    print("done superpath-to-enzyme")
    #print(super_enzyme_df["count"])
    
    #we have superpathway -> enzyme now
    
    path_df = pd.read_csv(ec_pathway_file, header=None, names=None, sep = '\t', skip_blank_lines = False)
    path_df.columns = ["path", "EC"]
    path_df["path"] = path_df["path"].apply(lambda x: x.split(":")[1])  
    path_df = path_df.loc[path_df["path"].str.contains("map")]
    #path_df.to_csv("path_df.csv", mode="w") 
    
    super_df = pd.read_csv(pathway_superpathway_file, sep = ',', skip_blank_lines = False)
    super_df.drop(columns = ["Pathway"], inplace=True)
    super_df.index = super_df["Pathway ID"]
    #print(super_df)
    
    path_df = path_df.join(super_df, on = 'path')
    #path_df.index = range(path_df.shape[0])
    path_df.drop(columns = ['Pathway ID', 'path'], axis=1, inplace=True)
    #print(path_df)
    
    #path_df.to_csv("new_path.csv", mode="w")
    enzyme_super_df = path_df
    enzyme_super_df.index = enzyme_super_df["EC"]
    print(enzyme_super_df)
    
    #print(rpkm_df)
    
    #This df now contains all genes -> enzymes, listed with their superpathway
    actual_read_df = rpkm_df[["GeneID", "EC#"]]
    #actual_read_df["EC#"] = actual_read_df["EC#"].apply(lambda x: "ec:"+x)
    actual_read_df["EC#"] = "ec:" + actual_read_df["EC#"]
    actual_read_df = actual_read_df.join(enzyme_super_df, on = "EC#")
    
    #print(actual_read_df)
    actual_read_df.to_csv("rpkm_match.csv", mode="w")
    actual_read_df.dropna(inplace=True)
    #now we need to figure out what the coverage is for the superpathway
    #superpath_coverage_df = actual_read_df.loc[actual_read_df.Superpathway is not None]
    #print(superpath_coverage_df)
    
    actual_read_df.drop(columns = ["GeneID", "EC#"], axis=1, inplace=True)
    actual_read_df.sort_values(by = ["Superpathway", "EC"], inplace = True)
    actual_read_df = actual_read_df[~actual_read_df.index.duplicated(keep='first')] #removes the weird copy effect we're seeing with the same index multiple times
    #each line is a new, different gene.  many genes can be a part of an enzyme.  
    #actual_read_df.index = actual_read_df.groupby("Superpathway").cumcount()
    print(actual_read_df.iloc[0:10])
    actual_read_df.to_csv("pathway_summary.csv", mode="w")
    