import pandas as pd
import sys


if __name__ == "__main__":
    ec_pathway_file = sys.argv[1]   #path->EC file
    rpkm_table_file = sys.argv[2]   #EC-> gene file
    pathway_superpathway_file = sys.argv[3] #path->superpath file
    path_df = pd.read_csv(ec_pathway_file, header=None, names=None, sep = '\t', skip_blank_lines = False)
    path_df.columns = ["path", "EC"]
    
    path_df["path"] = path_df["path"].apply(lambda x: x.split(":")[1])  
    
    path_df = path_df.loc[path_df["path"].str.contains("map")]
    path_df.index = path_df.groupby(["path"]).cumcount()
    temp_columns = path_df.index
    path_df = path_df.pivot(values='EC', columns = 'path')
    path_df = path_df.T
    #print(path_df)
    
    rpkm_df = pd.read_csv(rpkm_table_file, sep = '\t', skip_blank_lines = False)
    super_df = pd.read_csv(pathway_superpathway_file, sep = ',', skip_blank_lines = False)
    super_df.drop(columns = ["Pathway"], inplace=True)
    super_df.index = super_df.groupby("Superpathway").cumcount()
    super_df = super_df.pivot(values = 'Pathway ID', columns = 'Superpathway')
    super_df = super_df.T
    #print(super_df | path_df)
    print(super_df.iloc[0:10])
    print("--------------------------------------------")
    print("path df")
    print(path_df.iloc[0:10])
    #print(super_df)
    
    #pathway_ec_df = path_df.loc[path_df["path"] == super_df["Pathway ID"]]
    #print(pathway_ec_df)
    """
    super_pathways = list(super_df["Superpathway"].unique())
    path_list = set(super_df.loc[super_df["Superpathway"] == super_pathways[0]].values.flatten())
    path_list.discard(super_pathways[0])
    path_list = list(path_list)
    superpath_df = pd.DataFrame.from_dict({super_pathways[0]: path_list})
    for item in super_pathways:
        if(item == super_pathways[0]):
            continue
        temp_path_list = set(super_df.loc[super_df["Superpathway"] == item].values.flatten())
        temp_path_list.discard(item)
        path_list = list(temp_path_list)
        temp_df = pd.DataFrame.from_dict({item: path_list})
        #print(temp_df)
        superpath_df = superpath_df.merge(temp_df, left_index = True, right_index = True, how = 'outer')
    
    #superpath_df.dropna(inplace=True)
    superpath_df.fillna('0', inplace=True)
    superpath_df = superpath_df.T
    superpath_df.to_csv("superpath_redrawn.csv", header=False) 
    print(superpath_df)
    """
    
    """
    print("===============================")
    print("ec paths:")
    print(path_df.iloc[0:10])
    print("==============================")
    print("rpkm:")
    print(rpkm_df.iloc[0:10])
    """
    
    #print("===================================")
    #print("superpathway file")
    #print(super_df.iloc[0:10])
    
    