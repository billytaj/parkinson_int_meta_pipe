import os
import sys
from datetime import datetime as dt

#this code imports the gene_report and gets the split proportion

def import_gene_report(gene_report_path):

    gene_segment_dict = dict()
    with open(gene_report_path, "r") as gene_report:
        contig_name = ""
        contig_segment_name_list = []
        start_of_loop = True
        gene_segment_sum = 0
        
        
        
        for line in gene_report:
            
            cleaned_line = line.strip("\n")
            cleaned_line = cleaned_line.split(" ")
            
            cleaned_line = [i for i in cleaned_line if i]
            if(len(cleaned_line) > 1):
                #print(cleaned_line)
                if(cleaned_line[0] == "FASTA"):
                    if(start_of_loop == False):
                        for item in contig_segment_name_list:
                            #print("gene segment length:", gene_segment_dict[item], "gene segment sum:", gene_segment_sum)
                            
                            gene_segment_dict[item] = gene_segment_dict[item] / gene_segment_sum
                        contig_segment_name_list[:] = []
                    contig_name = cleaned_line[3]
                    #print("contig:", contig_name)
                    start_of_loop = False
                    gene_segment_sum = 0
                if(cleaned_line[0] == "Model" or cleaned_line[0] == "Gene" or cleaned_line[0] == "#" or cleaned_line[0] == "Predicted"):
                    continue
                if(cleaned_line[0].isdigit()):
                    #print(cleaned_line)
                    contig_segment_name = "gene_" + cleaned_line[0] + "|" + contig_name
                    #print("contig segment name:", contig_segment_name)
                    gene_segment_sum += int(cleaned_line[4])
                    #print("gene segment sum:", gene_segment_sum)
                    contig_segment_name_list.append(contig_segment_name)
                    gene_segment_dict[contig_segment_name] = int(cleaned_line[4])
        for item in contig_segment_name_list:
        #    print("gene segment length:", gene_segment_dict[item], "gene segment sum:", gene_segment_sum)
            gene_segment_dict[item] = gene_segment_dict[item] / gene_segment_sum
        contig_segment_name_list[:] = []
        
    return gene_segment_dict

def import_contig_map(contig_map_path):
    contig_map_dict = dict()
    with open(contig_map_path, "r") as contig_map:
        
        for line in contig_map:
            line_split = line.split("\t")
            contig_name = line_split[0]
            read_count = line_split[1]
            contig_map_dict[contig_name] = read_count
    return contig_map_dict


def convert_gene_segment_dict(gene_segment_dict, contig_map_dict):
    for item in gene_segment_dict:
        contig_name = item.split("|")[1]
        gene_segment_percent = gene_segment_dict[item]
        contig_reads = contig_map_dict[contig_name]
        gene_segment_dict[item] = int(round(float(contig_reads) * gene_segment_percent))


if __name__ == "__main__":
    
    gene_report_path = sys.argv[1]
    contig_map_path = sys.argv[2]
    gene_segment_dict = import_gene_report(gene_report_path)
    contig_map_dict = import_contig_map(contig_map_path)
    
    convert_gene_segment_dict(gene_segment_dict, contig_map_dict)
            
    print(dt.today(), "converting")        
    for item in gene_segment_dict:
        print(item, gene_segment_dict[item])
        