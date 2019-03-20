import pandas as pd
import os 
import sys


def parse_hum2_files(name):
    hum2_df = pd.read_csv(name, sep='\t', error_bad_lines=False)
    hum2_df = hum2_df.iloc[::2, :]
    hum2_df = hum2_df["# Gene Family"].str.split("|", 0, expand=True)
    hum2_df.columns = ["GeneID", "junk"]
    hum2_df.drop(columns=["junk"], inplace = True)
    hum2_df.reset_index(drop = True, inplace = True)
    hum2_df.drop_duplicates(inplace = True)
    return hum2_df


if __name__ == "__main__":
    mpro_rpkm               = sys.argv[1]
    pair_1_hum2_rpkm        = sys.argv[2]
    singletons_hum2_rpkm    = sys.argv[3]
    raw_hum2_rpkm           = sys.argv[4]
    
    
    mpro_df = pd.read_csv(mpro_rpkm, sep='\t', error_bad_lines=False)
    gene_id_df = mpro_df["GeneID"].str.split("|", 0, expand=True)
    cols = [0, 1, 2, 3, 4, 5, 6, "GeneID", 8]
    #gene_id_df.drop(gene_id_df.columns[cols], axis=1, inplace=True)
    gene_id_df.columns = cols#["GeneID"]
    #gene_id_df.drop_duplicates(inplace = True)
    gene_id_df.reset_index(drop = True, inplace = True)
    
    unique_gene_id_df = gene_id_df.drop_duplicates()
    print("UNIQUE gene ID DF")
    print(unique_gene_id_df)
    
    print("REGULAR")
    
    print(gene_id_df)
    
    pair_1_df       = parse_hum2_files(pair_1_hum2_rpkm)
    singletons_df   = parse_hum2_files(singletons_hum2_rpkm)
    raw_df          = parse_hum2_files(raw_hum2_rpkm)
    #print("pair 1")
    #print(pair_1_df)
    #print("singletons")
    #print(singletons_df)
    #print(raw_df)
    
    
    cleaned_hum2_df = pd.merge(pair_1_df, singletons_df, how='outer', on = "GeneID")
    cleaned_hum2_df.reset_index(drop = True, inplace = True)
    
    print("cleaned")
    print(cleaned_hum2_df)
    
    unique_hum2_df = cleaned_hum2_df.drop_duplicates()
    print("humann2 unique")
    print(unique_hum2_df)
    
    #print("RAW")
    #print(raw_df)
    
    
    like_findings_df = gene_id_df[gene_id_df.GeneID.isin(cleaned_hum2_df.GeneID)]
    print("HUMANN2 and MetaPro caught")
    print(like_findings_df)
    
    mpro_missed_df = cleaned_hum2_df[~cleaned_hum2_df.GeneID.isin(gene_id_df.GeneID)]
    print("HUMANN2 caught | MetaPro missed")
    print(mpro_missed_df)
    
    
    hum2_missed_df = gene_id_df[~gene_id_df.GeneID.isin(cleaned_hum2_df.GeneID)]
    print("MetaPro caught | Humann2 missed")
    print(hum2_missed_df)