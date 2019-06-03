#this pulls all of the EC charts into 1 large table. and also averages it

import os
import sys
import pandas as pd

if __name__ == "__main__":
    list_of_paths = sys.argv[1]
    paths_wanted_file = open(list_of_paths, "r")
    list_of_items_to_search_for = []
    for line in paths_wanted_file:
        list_of_items_to_search_for.append(line)
    
    final_results = []
    list_dir = os.listdir(os.getcwd())
    for item in list_dir:
        if(item.endswith(".csv")):
            print("---------------------------------")
            print("WORKING ON:", item)
            #file_path = os.path.join(os.getcwd(), item, str(os.listdir(item)[0]))
            file_df = pd.read_csv(item)
            frames_list = []
            for line in list_of_items_to_search_for:
                print("searching for:", line)
                path_to_find = line
                if(path_to_find.endswith("\n")):
                    path_to_find = path_to_find[:-1]
                selected_df = file_df[file_df["Pathway"] == path_to_find]
                #print(line)
                print(selected_df)
                frames_list.append(selected_df)
            if(frames_list):
                result_df = pd.concat(frames_list)
                result_df["source"] = item
                final_results.append(result_df)
                print(result_df)
            
            #print(os.listdir(item))
    final_collection_df = pd.concat(final_results)
    final_collection_df = final_collection_df.groupby("Pathway", as_index = False).mean()
    
    final_collection_df = final_collection_df.round()
    final_collection_df["remaining"] = 1
    final_collection_df["remaining"] = final_collection_df["remaining"].mask(final_collection_df["remaining"] > 0, final_collection_df["num_Ecs"] - (final_collection_df["both"] + final_collection_df["humann2"] + final_collection_df["mpro"]))
    print("==============================")
    print("FINAL COLLECTION")
    print(final_collection_df)
    
    final_collection_df.to_csv("final_collection.csv", mode = "w", index = False)