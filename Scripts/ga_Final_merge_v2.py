#This code performs all of the merging needed for GA to be completed.

import os
import sys
import multiprocessing as mp
from datetime import datetime as dt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import time

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
    
def translate_gene_segement_map(gene_segment_dict, contig_map_dict):
    for item in gene_segment_dict:
        contig_name = item.split("|")[1]
        gene_segment_percent = gene_segment_dict[item]
        contig_reads = contig_map_dict[contig_name]
        gene_segment_dict[item] = int(round(float(contig_reads) * gene_segment_percent))
    
def import_contig_map(contig_map_path):
    contig_map_dict = dict()
    with open(contig_map_path, "r") as contig_map:
        
        for line in contig_map:
            line_split = line.split("\t")
            contig_name = line_split[0]
            read_count = line_split[1]
            contig_map_dict[contig_name] = read_count
    return contig_map_dict

def export_proteins(diamond_proteins_file, gene_trans_dict, final_proteins):
    with open(diamond_proteins_file, "a") as diamond_proteins:
        for item in gene_trans_dict:
            SeqIO.write(item, diamond_proteins, "fasta")
            
def scrub_duplicates(fasta_in, fasta_out):
    imported_seqs = dict()
    header = "None"
    seq_body = "None"
    skip_flag = False
    #inspect_key = ">gi|483984714|ref|NZ_KB892660.1|:c111601-110108|33035|g__Blautia.s__Blautia_producta|UniRef90_unknown|UniRef50_R5TWD5"
    with open(fasta_in, "r") as in_seq:
        for line in in_seq:
            cleaned_line = line.strip("\n")
            
            if(cleaned_line.startswith(">")):
                if(header == "None"):
                    header = cleaned_line
                if(header in imported_seqs):
                    #header already imported.  skip it
                    header = "None"
                    skip_flag = True
                else:
                    skip_flag = False
                    if(seq_body != "None"):
                        #gene data gathered
                        imported_seqs[header] = seq_body
                        #if(header == inspect_key):
                        #    print(dt.today(), "inspect key found:", header)#, seq_body)
                            
                        
                        seq_body = "None"
                        header = cleaned_line
                        
            else:
                if(skip_flag):
                    seq_body = "None"
                    
                else:
                    if(seq_body == "None"):
                        seq_body = cleaned_line
                    else:
                        seq_body += cleaned_line
    
    with open(fasta_out, "w") as out_seq:
        for item in imported_seqs:
            out_line = item + "\n" + imported_seqs[item] + "\n"
            out_seq.write(out_line)    
def convert_genes_to_proteins(mapped_gene_file):#, section, gene_trans_dict):
    # WRITE OUTPUT: BWA&BLAT&DMD-aligned gene/protIDs and aa seqs.  It's for a downstream tool.     
    # (.faa; fasta-format):
    # convert previously mapped genes to the AA format.
    #apparently seqIO can't handle duplicates in its file... it's kinda bullshit.
    
    
    gene_seqs = SeqIO.index(mapped_gene_file,"fasta")           # key=geneID, value=SeqRecord
    #inspect_key = "gi|483984714|ref|NZ_KB892660.1|:c111601-110108|33035|g__Blautia.s__Blautia_producta|UniRef90_unknown|UniRef50_R5TWD5"
    
    #gene_trans= []
    gene_trans_dict = dict()
    for gene in gene_seqs:                                  # Take each BWA&BLAT-aligned genes

        translated_seqRecord = SeqRecord(seq= gene_seqs[gene].seq.translate(stop_symbol=""), id= gene_seqs[gene].id, description= gene_seqs[gene].description)
        gene_id = gene
        gene_trans_dict[gene_id] = translated_seqRecord.seq                                                        #  and translate its SeqRecord sequence to aa.
        if(len(gene_trans_dict[gene_id]) * 3 > len(gene_seqs[gene].seq)):
            print(dt.today(), gene_id, "protein longer than bp seq.  protein;", len(gene_trans_dict[gene_id]), "bp:", len(gene_seqs[gene].seq))
            sys.exit("this is real bad")

    return gene_trans_dict
        

def concatenate_gene_maps_v2(path, gene_map_dict, context):
    #merge all gene maps given a directory
    #the gene length stays constant.  The number of reads do not.
    path_contents = os.listdir(path)
    duplicates_found = False
    #gene_map_dict = dict()
    for content in path_contents:
        item = os.path.join(path, content)
        base_name = os.path.basename(item)
        #print("base name:", base_name)
        if(base_name.endswith("_gene_map.tsv")):
            #print("context:", context, "checking:", item)
            #time.sleep(1)
            if(base_name.startswith(context)):
                #print("context:", context, "looking at:", item)
                with open(item, "r") as split_gene_map:
                    for line in split_gene_map:
                        line_split = line.strip("\n").split("\t")
                        number_of_reads = int(line_split[2])
                        gene_length = int(line_split[1])
                        gene_name = line_split[0]
                        reads = line_split[3:]
                  
                        if(gene_name in gene_map_dict):
                            old_pack = gene_map_dict[gene_name]
                            old_reads = old_pack[2:]
                            old_gene_length = int(old_pack[0])
                            
                            new_reads = list(set(old_reads + reads))
                            if(len(new_reads) < (len(old_reads) + len(reads))):
                                duplicates_found = True
                                print("---------------------------------")
                                print(dt.today(), gene_name, "duplicates removed", len(new_reads), "<-", len(old_reads) + len(reads))
                                print(gene_name, "new:", new_reads)
                                print(gene_name, "old:", old_reads)
                                print(gene_name, "added:", reads)
                                #time.sleep(3)
                            #print(gene_name, "old reads:", len(old_reads), "current reads:", len(reads), "combined:", len(new_reads))
                            #time.sleep(2)
                            
                            new_number_of_reads = len(new_reads)
                            final_gene_length = gene_length
                            if(int(old_gene_length) != int(gene_length)):
                                print("gene lengths for same gene in table do not match. This shouldn't happen")
                                print("taking larger length:", old_gene_length, gene_length)
                                print("for:", gene_name)
                                if(int(old_gene_length) > int(gene_length)):
                                    final_gene_length = old_gene_length
                                else:
                                    final_gene_length = gene_length
                            
                            gene_map_dict[gene_name] = [final_gene_length, new_number_of_reads] + new_reads
                            #print("merged entry:", gene_name, gene_map_dict[gene_name])
                        else:
                            
                            gene_map_dict[gene_name] = [gene_length, number_of_reads] + reads
    if(duplicates_found):
        print(dt.today(),"duplicates found")
    else:
        print(dt.today(), path, context, "OK! no duplicates found")
                

    
def merge_fastas(path_0, path_1, section, header, extension):
    #it's just a simple concatenation job. that also checks for duplicates
    IDs_used = list()
    skip_this_line = False
    path_contents = os.listdir(path_0)
    final_fasta_file = os.path.join(path_1, header + "_" + section + extension)
    with open(final_fasta_file, "w") as final_fasta:
        for content in path_contents:
            item = os.path.join(path_0, content)
            #print("fasta merge:", item)
            if(content.startswith(section)) and (content.endswith(extension)):
                with open(item, "r") as sample_fasta:
                    for line in sample_fasta:
                        if(line.startswith(">")):
                            ID = line
                            if(ID in IDs_used):
                                skip_this_line = True
                            else:
                                IDs_used.append(ID)
                                skip_this_line = False
                                
                        if not(skip_this_line):        
                            final_fasta.write(line)
                        #else:
                        #    print(dt.today(), ID, "already written.  skipping")
                            
    print(dt.today(), "find it at:", final_fasta_file)
    return final_fasta_file
    
    
def export_gene_map(gene_map, export_path, header = None):
    final_gene_map_file = ""
    if(header == None):
    
        final_gene_map_file = os.path.join(export_path, "gene_map.tsv")
    else:
        final_gene_map_file = os.path.join(export_path, header + "_gene_map.tsv")
    with open(final_gene_map_file, "w") as gene_map_out:
        #out_line = "geneID" + "\t" + "gene length" + "\t" + "Reads"
        for gene_name_key in gene_map:
            gene_entry = gene_map[gene_name_key]
            gene_ID = gene_name_key
            gene_length = gene_entry[0]
            number_of_reads = gene_entry[1]
            
            out_line = gene_ID + "\t" + str(gene_length) + "\t" + str(number_of_reads)
            for item in gene_entry[2:]:
                out_line += "\t" + item
            out_line += "\n"
            gene_map_out.write(out_line)
            
def make_merge_fasta_process(process_store, path_0, path_1, header, tail):
    section = ["pair_1", "pair_2", "contigs", "singletons"]
    for item in section:
        merge_process = mp.Process(target = merge_fastas, args = (path_0, path_1, item, header, tail))
        merge_process.start()
        process_store.append(merge_process)

def make_convert_proteins_process(process_store, path, gene_trans_dict):
    path_contents = os.listdir(path)
    
    for item in path_contents:
        if(item.endswith(".fna")):
            full_path = os.path.join(os.path.abspath(path), item)
            basename = os.path.splitext(item)[0]
            #print(full_path)
            #section = ["BWA_annotated", "BLAT_annotated"]
            gene_to_protein_process = mp.Process(target = convert_genes_to_proteins, args = (basename, full_path, gene_trans_dict))
            gene_to_protein_process.start()
            process_store.append(gene_to_protein_process)
        
def make_second_merge_process(process_store, path_0, path_1, header, tail):
    section = ["BWA_annotated", "BLAT_annotated"]
    for item in section:
        merge_process = mp.Process(target = merge_fastas, args = (path_0, path_1, item, header, tail))
        merge_process.start()
        process_store.append(merge_process)
        
        
def merge_all_proteins(path, gene_transcripts_dict, export_path):
    #this just exports.
    all_proteins_path = os.path.join(export_path, "all_proteins.faa")
    skip_this_line = False
    IDs_used = list()
    with open(all_proteins_path, "w") as proteins_out:
        for item in os.listdir(path):
            
            if(item.endswith(".faa")):
                full_path = os.path.join(os.path.abspath(path), item)
                ID = "None"
                with open(full_path, "r") as dmd_proteins:
                    for line in dmd_proteins:
                    
                        if(line.startswith(">")):
                            ID = line
                            if(ID in IDs_used):
                                #print(dt.today(), "protein already written. skipping")
                                skip_this_line = True
                            else:
                                IDs_used.append(ID)
                                skip_this_line = False
                                #print("IDs used:", IDs_used)
                                #time.sleep(0.5)
                            
                        if not skip_this_line:
                            proteins_out.write(line)
                        else:
                            print("not writing:", ID)
                            #print(line)
                            #time.sleep(0.5)
        for gene in gene_transcripts_dict:
            #print("EXPORTING gene transcript:", gene, gene_transcripts_dict[gene])
            out_line = ">" + gene + "\n" 
            out_line += str(gene_transcripts_dict[gene]) + "\n"
            
            proteins_out.write(out_line)
        #for gene in gene_transcripts_list:
        #    out_line = str(gene.id) + "\n" + str(gene.seq) + "\n"
        #    proteins_out.write(out_line)
        #SeqIO.write(gene_transcripts_list, proteins_out, "fasta")
        
def handle_final_proteins(final_path, export_path):
    #just a function dump so we can do this step in parallel with other things

    #scrub for duplicates (because seqIO index doesn't like duplicates)
    bwa_blat_proteins_file = merge_fastas(final_path, final_path, "all", "BWA_BLAT_proteins", ".fna")
    unique_fna_file = "unique_" + os.path.basename(bwa_blat_proteins_file)
    real_unique_fna_file = os.path.join(final_path, unique_fna_file)
    scrub_duplicates(bwa_blat_proteins_file, real_unique_fna_file)
    print(dt.today(), real_unique_fna_file, "should have no duplicates")
    #sys.exit("Duplicates removed")
    gene_transcripts_dict = convert_genes_to_proteins(real_unique_fna_file)
    #print("GENE TRANSCRIPT:", len(gene_transcripts_list))
    #for item in gene_transcripts_list:
    #    print(item)
    
    merge_all_proteins(final_path, gene_transcripts_dict, export_path)
    
def translate_contig_segments_to_reads(gene_map_dict, gene_segment_dict):
    #update the read count with the gene segments' proper read counts.   
    for item in gene_segment_dict:
        print(type(item), item)
        
    for gene in gene_map_dict:
        
        reads = gene_map_dict[gene][2:]
        number_of_reads = gene_map_dict[gene][1]
        gene_length = gene_map_dict[gene][0]
        for read_ID in reads:
            if(read_ID.startswith("gene")):
                read_ID_list = read_ID.split("|")
                tail = read_ID_list[-1]
                tail = tail.split(">")[-1]
                gene_segment_key = read_ID_list[0] + "|" + tail
                #print(dt.today(), "working with:", gene_segment_key)
                if gene_segment_key in gene_segment_dict:
                    #print("Read ID:", gene_segment_key, "length:", gene_segment_dict[gene_segment_key])
                    if(gene_segment_dict[gene_segment_key] > 1):
                        #print("read count added to count:", gene_segment_dict[gene_segment_key], "+", number_of_reads, "-1")
                        number_of_reads += gene_segment_dict[gene_segment_key] - 1
                        
                    gene_map_dict[gene] = [gene_length, number_of_reads] + reads
                    
def merge_fastas_v2(path, header, extension, name):
    list_of_file = os.listdir(path)
    combined_file = os.path.join(path, "all_" + header + extension)
    with open(combined_file, "w") as out_file:            
        for item in list_of_file:
            if(item.startswith(header) and item.endswith(extension)):
                
                working_file = os.path.join(path, item)
                print("file found matching pattern:", working_file)
                with open(working_file, "r") as in_file:
                    for line in in_file:
                        out_file.write(line)
    print(dt.today(), "combined file:", combined_file)
    return combined_file
                
def merge_dicts(list_of_dicts):
    #gene dict is: key: gene_name, val: [gene length, number of reads, reads]
    final_dict = {}
    
    for item in list_of_dicts:
        #if final dict is empty
        if(not bool(final_dict)):
            print("final dict empty")
            final_dict = item
            continue
        else:
            print("final dict not empty")
            for gene in item:
                if(gene in final_dict):
                    gene_length = item[gene][0]
                    
                    reads = item[gene][2:]
                    old_reads = final_dict[gene][2:]
                    new_reads = list(set(reads + old_reads))
                    new_entry = [gene_length, len(new_reads)] + new_reads
                    
                    if(len(new_reads) < (len(old_reads) + len(reads))):
                        print(gene, "dupes removed | old:", len(old_reads), "added:", len(reads), "combined:", len(new_reads))
                        print("new:", new_reads)
                        print("old:", old_reads)
                        print("added:", reads)
                        print("=====================================")
                    #else:
                    #    print(gene, "no dupes found | old:", len(old_reads), "added:", len(reads), "combined:", len(new_reads))
                    final_dict[gene] = new_entry
                    time.sleep(1)
                    
            else:
                #add new gene to map
                final_dict[gene] = item[gene]
    return final_dict
 
if __name__ == "__main__":
    assemble_path           = sys.argv[1]
    bwa_path                = sys.argv[2]
    blat_path               = sys.argv[3]
    diamond_path            = sys.argv[4]
    final_path              = sys.argv[5]
    export_path             = sys.argv[6]
    #diamond_proteins_file   = sys.argv[5]
    
    
    
    #gene_report_path = os.path.join(assemble_path, "data", "1_mgm", "gene_report.txt")
    gene_report_path = os.path.join(assemble_path, "gene_report.txt")
    
    contig_map_path = os.path.join(bwa_path, "contig_map.tsv")
    
    gene_segment_dict = import_gene_report(gene_report_path)
    contig_map_dict = import_contig_map(contig_map_path)
    translate_gene_segement_map(gene_segment_dict, contig_map_dict)
    
    
    
    manager = mp.Manager()
    mgr_bwa_pair_1_gene_map     = manager.dict()
    mgr_bwa_pair_2_gene_map     = manager.dict()
    mgr_bwa_singletons_gene_map = manager.dict()
    mgr_bwa_contig_gene_map     = manager.dict()

    mgr_blat_pair_1_gene_map    = manager.dict()
    mgr_blat_pair_2_gene_map    = manager.dict()
    mgr_blat_singletons_gene_map= manager.dict()
    mgr_blat_contig_gene_map    = manager.dict()
    
    mgr_dia_pair_1_gene_map     = manager.dict()
    mgr_dia_pair_2_gene_map     = manager.dict()
    mgr_dia_singletons_gene_map = manager.dict()
    mgr_dia_contig_gene_map     = manager.dict()
    
    
    process_store = []
    
    
    #merge the leftover reads, and annotated genes (for protein translation)
    #make_merge_fasta_process(process_store, diamond_path, export_path, "GA_leftover", ".fasta")  #leftover reads
    #make_merge_fasta_process(process_store, blat_path, final_path, "BLAT_annotated", ".fna")    #BLAT genes
    #make_merge_fasta_process(process_store, bwa_path, final_path, "BWA_annotated", ".fna")      #BWA genes
    #make_merge_fasta_process(process_store, diamond_path, final_path, "dmd", ".faa")            #DIAMOND proteins
    
    
    #merge the gene maps, but group them by data context (pair1, pair2, singletons, contigs) for each tool (BWA, BLAT, DIAMOND)
    process = mp.Process(target = concatenate_gene_maps_v2, args = (bwa_path, mgr_bwa_pair_1_gene_map, "pair_1"))
    process.start()
    process_store.append(process)
    """
    process = mp.Process(target = concatenate_gene_maps_v2, args = (bwa_path, mgr_bwa_pair_2_gene_map, "pair_2"))
    process.start()
    process_store.append(process)
    
    process = mp.Process(target = concatenate_gene_maps_v2, args = (bwa_path, mgr_bwa_singletons_gene_map, "singletons"))
    process.start()
    process_store.append(process)
    
    process = mp.Process(target = concatenate_gene_maps_v2, args = (bwa_path, mgr_bwa_contig_gene_map, "contig"))
    process.start()
    process_store.append(process)
    """
    #----------------------------------------------------------------------
    process = mp.Process(target = concatenate_gene_maps_v2, args = (blat_path, mgr_blat_pair_1_gene_map, "pair_1"))
    process.start()
    process_store.append(process)
    """
    process = mp.Process(target = concatenate_gene_maps_v2, args = (blat_path, mgr_blat_pair_2_gene_map, "pair_2"))
    process.start()
    process_store.append(process)
    
    process = mp.Process(target = concatenate_gene_maps_v2, args = (blat_path, mgr_blat_singletons_gene_map, "singletons"))
    process.start()
    process_store.append(process)
    
    process = mp.Process(target = concatenate_gene_maps_v2, args = (blat_path, mgr_blat_contig_gene_map, "contig"))
    process.start()
    process_store.append(process)
    """
    #----------------------------------------------------------------------
    process = mp.Process(target = concatenate_gene_maps_v2, args = (diamond_path, mgr_dia_pair_1_gene_map, "pair_1"))
    process.start()
    process_store.append(process)
    """
    process = mp.Process(target = concatenate_gene_maps_v2, args = (diamond_path, mgr_dia_pair_2_gene_map, "pair_2"))
    process.start()
    process_store.append(process)
    
    process = mp.Process(target = concatenate_gene_maps_v2, args = (diamond_path, mgr_dia_singletons_gene_map, "singletons"))
    process.start()
    process_store.append(process)
    
    process = mp.Process(target = concatenate_gene_maps_v2, args = (diamond_path, mgr_dia_contig_gene_map, "contig"))
    process.start()
    process_store.append(process)
    """
    
    
    print(dt.today(), "concatenating gene maps")
    for item in process_store:
        item.join()
    process_store[:] = []
    
    #-----------------------------------------------------------------------
    #convert 
    bwa_pair_1_gene_map     = dict(mgr_bwa_pair_1_gene_map)
    #bwa_pair_2_gene_map     = dict(mgr_bwa_pair_2_gene_map)
    #bwa_singleton_gene_map  = dict(mgr_bwa_singletons_gene_map)
    bwa_contig_gene_map     = dict(mgr_bwa_contig_gene_map)
    
    blat_pair_1_gene_map    = dict(mgr_blat_pair_1_gene_map)
    #blat_pair_2_gene_map    = dict(mgr_blat_pair_2_gene_map)
    #blat_singleton_gene_map = dict(mgr_blat_singletons_gene_map)
    blat_contig_gene_map    = dict(mgr_blat_contig_gene_map)
    
    dia_pair_1_gene_map     = dict(mgr_dia_pair_1_gene_map)
    #dia_pair_2_gene_map     = dict(mgr_dia_pair_2_gene_map)
    #dia_singleton_gene_map  = dict(mgr_dia_singletons_gene_map)
    dia_contig_gene_map     = dict(mgr_dia_contig_gene_map)
    
    
    pair_1_gene_map_list    = [bwa_pair_1_gene_map, blat_pair_1_gene_map, dia_pair_1_gene_map]
    #pair_2_gene_map_list    = [bwa_pair_2_gene_map, blat_pair_2_gene_map, dia_pair_2_gene_map]
    #singleton_gene_map_list = [bwa_singleton_gene_map, blat_singleton_gene_map, dia_singleton_gene_map]
    contig_gene_map_list    = [bwa_contig_gene_map, blat_contig_gene_map, dia_contig_gene_map]
    
    
    #sys.exit("los_ageless")
    pair_1_gene_map = merge_dicts(pair_1_gene_map_list)
    contig_gene_map = merge_dicts(contig_gene_map_list)
    
    print("contigs")
    for item in contig_gene_map:
        print(item, contig_gene_map[item])
    
    
    print(dt.today(), "done converting")
    #print(bwa_pair_1_gene_map)
    
    #all_dmd_file = merge_fastas_v2(final_path, "dmd", ".faa", "something")
    
    sys.exit("dialogue")
    
    
    #secondary combine on the genes (BWA and BLAT)
    make_second_merge_process(process_store, final_path, final_path, "all", ".fna")
    #make_merge_fasta_process(process_store, final_path, final_path, "dmd", ".faa")
    
    for item in process_store:
        item.join()
    process_store[:] = []
    
    
    
    
    finish_proteins_process = mp.Process(target = handle_final_proteins, args = (final_path, export_path))
    finish_proteins_process.start()
    process_store.append(finish_proteins_process)
    
    
    bwa_gene_map = dict(mgr_bwa_gene_map)
    blat_gene_map = dict(mgr_blat_gene_map)
    diamond_gene_map = dict(mgr_diamond_gene_map)
    
    #the maps contain contig segments that under-represent the count.  We need to convert.
    translate_contig_segments_to_reads(bwa_gene_map, gene_segment_dict)
    translate_contig_segments_to_reads(blat_gene_map, gene_segment_dict)
    translate_contig_segments_to_reads(diamond_gene_map, gene_segment_dict)
    
    
    gene_map_list = [bwa_gene_map, blat_gene_map, diamond_gene_map]
    #merge all gene maps
    final_gene_map_process = mp.Process(target = final_gene_map_merge, args = (gene_map_list, mgr_final_gene_map))
    final_gene_map_process.start()
    process_store.append(final_gene_map_process)
    
    for item in process_store:
        item.join()
    process_store[:] = []
    
    final_gene_map = dict(mgr_final_gene_map)
    
    export_gene_map(final_gene_map, export_path)
    
    print(dt.today(), "We're at the end")
    

        