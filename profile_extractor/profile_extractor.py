import csv
import os
import datetime
#This program takes in files ending in extension "_OUT" from qsub exports, extracts the walltime, 
#and displays it as a useful metric


def main():
    full_profile_list = list()
    total_time = 0
    for item in os.listdir():
        if("_OUT" in item):
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
    for item in full_profile_list:
        item[2] = str(round(((item[2] / total_time) * 100), 2)) + "%"
    full_profile_list.append(["total time:", str(datetime.timedelta(seconds=total_time))])    
    return full_profile_list
if __name__ == "__main__":
    os.chdir("C:/Users/Billy/Desktop/Profile_info")
    summarized_timing = main()
    for item in summarized_timing:
        print(item)