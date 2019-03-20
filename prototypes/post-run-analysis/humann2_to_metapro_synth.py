import pandas as pd
import os 
import sys


def parse_hum2_files(name):
    hum2_df = pd.read_csv(name, sep='\t', error_bad_lines=False)
    hum2_df = hum2_df["# Gene Family"].str.split("|", 0, expand=True)
    #hum2_df[1] = hum2_df[1].str.split(".", 0, expand=True)
    hum2_df.dropna(inplace=True)
    hum2_df.reset_index(drop=True, inplace=True)
    hum2_df.columns = ["uniref", "taxa"]
    
    hum2_df["GeneID"] = hum2_df["uniref"] + "|" + hum2_df["taxa"]
    hum2_df.drop(columns = ["uniref", "taxa"], inplace = True)
    return hum2_df


    
def parse_metapro(name):
    mpro_df = pd.read_csv(name, sep='\t', error_bad_lines = False)
    mpro_df = mpro_df["GeneID"].str.split("|", 0, expand=True)
    cols = [0, 1, 2, 3, 4, 5, "taxa", "uniref90", 8]
    mpro_df.columns = cols
    
    mpro_df["GeneID"] = mpro_df["uniref90"] + "|" + mpro_df["taxa"]
    mpro_df.dropna(inplace=True)
    mpro_df.drop(columns  = cols, inplace = True)
    
    return mpro_df
    
if __name__ == "__main__":
    
    hum2_df = parse_hum2_files(sys.argv[1])
    mpro_df = parse_metapro(sys.argv[2])
    
    #taxa_df = mpro_df["taxa"].str.split(".", 0, expand=True)
    #mpro_df.drop(columns = ["taxa"])
    print("HUMANN2")
    print(hum2_df)
    
    print("===============================")
    print("MetaPro")
    print(mpro_df)
    
    
    print("=======================================")
    print("humann2 caught | metapro missed")
    
    hum2_hit_df = hum2_df[~hum2_df.GeneID.isin(mpro_df.GeneID)]
    print(hum2_hit_df)
    
    print("====================================")
    print("humann2 missed | metapro caught")
    mpro_hit_df = mpro_df[~mpro_df.GeneID.isin(hum2_df.GeneID)]
    print(mpro_hit_df)