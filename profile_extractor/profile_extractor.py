import csv
import os
import datetime
import sys
#This program takes in files ending in extension "_OUT" from qsub exports, extracts the walltime, 
#and displays it as a useful metric

#this tool needs to be upgraded to consider overlappings.
#some stages do try to run in parallel

class cloud_time():
    min_start = 0
    max_end = 0
    avg = 0 
    sample = 0    

def package_min_max(message, min, max):
    return [message, str(datetime.timedelta(seconds = (max - min))), (max - min) ]

def get_min_max(item, min_start_time, max_end_time, sum_time):
    # of a collection of common profiling files, find the min, and max of all of them
    process_duration = 0
    sample_counter = 0
    with open(item, "r") as time_file:
        time_reader = csv.reader(time_file, delimiter = ' ')
        
        for row in time_reader:
            
            if(len(row) > 0):
                
                if(row[0] == "Begin" and row[2] == "Prologue"):
                    
                    #print(row)
                    process_start_time = int(row[-1])
                    process_duration = process_start_time
                    if(min_start_time == 0):
                        min_start_time = process_start_time
                    elif(min_start_time > process_start_time):
                        min_start_time = process_start_time
                elif(row[0] == "End" and row[2] == "Epilogue"):
                    #print(row)
                    process_end_time = int(row[-1])
                    if(process_duration == 0 or process_duration >= process_end_time):
                        print("Something wrong here")
                        break
                    else:
                        process_duration = (process_end_time - process_duration)
                    if(max_end_time == 0):
                        max_end_time = process_end_time
                    elif(max_end_time < process_end_time):
                        max_end_time = process_end_time
                        
                
                    sum_time += process_duration        
    return min_start_time, max_end_time, sum_time
    
    
def main():
    full_profile_list = list()
    total_time = 0
    BLAT_start_time = 0
    BLAT_end_time = 0
    BLAT_avg_time = 0
    BLAT_sample_count = 0
    
    DIAMOND_start_time = 0
    DIAMOND_end_time = 0
    DIAMOND_avg_time = 0
    DIAMOND_sample_count = 0
    #unpaired_infernal_start = 0
    #unpaired_infernal_end = 0
    #paired_infernal_start = 0
    #paired_infernal_end = 0
    infernal_start_time = 0
    infernal_end_time = 0
    infernal_avg_time = 0
    infernal_sample_count = 0
    for item in os.listdir():
        # BLAT and DIAMOND jobs are launched independently.  Scinet allows multiple jobs at once. 
        # so we care about the min start, and the max end, which indicates how long the whole stage takes
        if("Annotate_BLAT" in item and "_OUT" in item and "Annotate_BLAT_Post" not in item):
            BLAT_start_time, BLAT_end_time, BLAT_avg_time = get_min_max(item, BLAT_start_time, BLAT_end_time, BLAT_avg_time)
            BLAT_sample_count += 1
            
        elif("Annotate_DMD" in item and "_OUT" in item and "Annotate_DMD_Post" not in item):
            DIAMOND_start_time, DIAMOND_end_time, DIAMOND_avg_time = get_min_max(item, DIAMOND_start_time, DIAMOND_end_time, DIAMOND_avg_time)
            DIAMOND_sample_count += 1
            
        elif("rRNA_Filter" in item and "_OUT" in item):
            infernal_start_time, infernal_end_time, infernal_avg_time = get_min_max(item, infernal_start_time, infernal_end_time, infernal_avg_time)
            infernal_sample_count += 1
            
        elif("_OUT" in item):
            with open(item, "r") as time_file:
                time_reader = csv.reader(time_file, delimiter = '\t')
                per_stage_profile = list()
                for time_row in time_reader:
                # the format is:  [name of qsub call, clock time, % of overall time]
                    if(len(time_row) > 0):
                        if(time_row[0] == "Job Name:"):
                            per_stage_profile.append(str(time_row[1]))
                        if(time_row[0] == "Resources:"):
                            extracted_time = time_row[1].split(",")[-1].split('=')[-1]
                            split_time = extracted_time.split(":")
                            split_hour = int(split_time[0]) * 3600
                            split_min = int(split_time[1]) * 60 
                            split_sec = int(split_time[2])
                            sum_sec = split_hour + split_min + split_sec
                            total_time += sum_sec
                            per_stage_profile.append(extracted_time)
                            per_stage_profile.append(sum_sec)
                    #time_row = time_row.strip('\t')
                    #if(len(time_row) > 0 and time_row[0] == "Resources:"):
                    #    print(time_row[1].split(",")[-1])
                full_profile_list.append(per_stage_profile)
                
    # needs total sum before we can apply the % breakdown
    total_time += (BLAT_end_time - BLAT_start_time) + (DIAMOND_end_time - DIAMOND_start_time) + (infernal_end_time - infernal_start_time)#(unpaired_infernal_end - unpaired_infernal_start) + (paired_infernal_end - paired_infernal_start)
    #paired_infernal_item = package_min_max("Paired infernal", paired_infernal_start, paired_infernal_end)
    #unpaired_infernal_item = package_min_max("Unpaired Infernal", unpaired_infernal_start, unpaired_infernal_end)
    infernal_item = package_min_max("Infernal", infernal_start_time, infernal_end_time)
    
    BLAT_item = package_min_max("BLAT annoation", BLAT_start_time, BLAT_end_time)
    DIAMOND_item = package_min_max("DIAMOND annotation", DIAMOND_start_time, DIAMOND_end_time)
    full_profile_list.extend((BLAT_item, DIAMOND_item, infernal_item))
    
    if(BLAT_sample_count == 0):
        BLAT_sample_count += 1
        
    if(DIAMOND_sample_count == 0):
        DIAMOND_sample_count += 1
        
    if(infernal_sample_count == 0):
        infernal_sample_count += 1
    print("infernal sample count:", infernal_sample_count)
    for item in full_profile_list:
        item[2] = str(round(((item[2] / total_time) * 100), 2)) + "%"
    full_profile_list.append(["total time:", str(datetime.timedelta(seconds=total_time))])    
    
    infernal_avg_pack = package_min_max("infernal average run:", 0, int(infernal_avg_time/infernal_sample_count))
    BLAT_avg_pack = package_min_max("BLAT avg run:", 0, int(BLAT_avg_time/BLAT_sample_count))
    DIAMOND_avg_pack = package_min_max("DIAMOND avg run:", 0, int(DIAMOND_avg_time/DIAMOND_sample_count))
    
    infernal_avg_pack[2] = "number of samples: ", infernal_sample_count
    
    full_profile_list.extend((infernal_avg_pack, BLAT_avg_pack, DIAMOND_avg_pack))
    return full_profile_list

if __name__ == "__main__":
    profile_dir = sys.argv[1]
    os.chdir("C:/Users/Billy/Desktop/Profile_info/" + profile_dir)
    
    summarized_timing = main()
    for item in summarized_timing:
        print(item)
    #other_thing()    