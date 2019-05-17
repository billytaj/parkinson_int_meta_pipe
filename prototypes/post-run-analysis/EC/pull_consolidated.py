import os
import sys
import pandas as pd

if __name__ == "__main__":
    list_of_paths = sys.argv[1]
    paths_wanted_file = open(list_of_paths, "r")
    
    
    list_dir = os.listdir(os.getcwd())
    for item in list_dir:
        if(item.startswith("5")):
        
            print("WORKING ON:", item)
            file_path = os.path.join(os.getcwd(), item, os.listdir(item)[0])
            file_df = pd.read_csv(file_path)
            frames_list = []
            for line in paths_wanted_file:
                path_to_find = line
                if(path_to_find.endswith("\n")):
                    path_to_find = path_to_find[:-1]
                selected_df = file_df[file_df["Pathway"] == path_to_find]
                #print(line)
                #print(selected_df)
                frames_list.append(selected_df)
            result_df = pd.concat(frames_list)
            result_df["source"] = item
            
            print(result_df)
            
            #print(os.listdir(item))
    