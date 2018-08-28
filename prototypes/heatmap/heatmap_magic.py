import pandas as pd
import sys


if __name__ == "__main__":
    ec_pathway_file = sys.argv[1]   #path->EC file
    rpkm_table_file = sys.argv[2]   #EC-> gene file
    pathway_superpathway_file = sys.argv[3] #path->superpath file
    path_df = pd.read_csv(ec_pathway_file, header=None, names=None, sep = '\t', skip_blank_lines = False)
    path_df.columns = ["path", "EC"]
    
    path_df["path"] = path_df["path"].apply(lambda x: x.split(":")[1])  
    
    path_df.index = path_df.groupby(["path"]).cumcount()
    temp_columns = path_df.index
    path_df = path_df.pivot(values='EC', columns = 'path')
    path_df = path_df.T
    path_df["EC_list"] = path_df[path_df.columns[:]].apply(lambda x: ",".join(x.dropna()), axis=1)
    path_df.drop(temp_columns, axis=1, inplace=True)
    #print(path_df)
    
    rpkm_df = pd.read_csv(rpkm_table_file, sep = '\t', skip_blank_lines = False)
    super_df = pd.read_csv(pathway_superpathway_file, sep = ',', skip_blank_lines = False)
    super_df.drop(columns = ["Pathway"], inplace=True)
    
    
    #super_df.to_csv("super_path.csv", mode=  "w")
    #path_df.to_csv("path_enzymes.csv", mode = "w")
    #print(super_df | path_df)
    
    super_df = super_df.join(path_df, how = 'left', on = 'Pathway ID')
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
    print(super_enzyme_df["count"])
    
    #we have superpathway -> enzyme now
    
    path_df = pd.read_csv(ec_pathway_file, header=None, names=None, sep = '\t', skip_blank_lines = False)
    path_df.columns = ["path", "EC"]
    path_df["path"] = path_df["path"].apply(lambda x: x.split(":")[1])  
    path_df = path_df.loc[path_df["path"].str.contains("map")]
     
    
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
    
    actual_read_df = rpkm_df[["GeneID", "EC#"]]
    actual_read_df["EC#"] = actual_read_df["EC#"].apply(lambda x: "ec:"+x)
    actual_read_df = actual_read_df.join(enzyme_super_df, on = "EC#")
    print(actual_read_df)
    actual_read_df.to_csv("rpkm_match.csv", mode="w")
    #print(enzyme_super_df[0:50])
    #print(super_df[0:10])
    