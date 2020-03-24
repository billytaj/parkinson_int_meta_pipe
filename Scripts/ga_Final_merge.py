#This code performs all of the merging needed for GA to be completed.

import os
import sys
import multiprocessing as mp
from datetime import datetime as dt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


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
    header = ""
    seq_body = ""
    with open(fasta_in, "r") as in_seq:
        for line in in_seq:
            if line.startswith(">"):
                if(header == ""):
                    header = line.strip("\n")
                else:
                    if(header in imported_seqs):
                        print("skipping duplicate header:", header)
                        header = ""
                        seq_body = ""
                    else:
                        imported_seqs[header] = seq_body
                        header = ""
                        seq_body = ""
            else:
                seq_body += line.strip("\n")
    
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
    
    
    gene_trans= []
    for gene in gene_seqs:                                  # Take each BWA&BLAT-aligned genes
        print("working on:", gene)
        try:
            gene_trans.append(SeqRecord(seq= gene_seqs[gene].seq.translate(stop_symbol=""), id= gene_seqs[gene].id, description= gene_seqs[gene].description))
                                                            #  and translate its SeqRecord sequence to aa.
        except Exception:
            pass

    return gene_trans
        
        
        
def concatenate_gene_maps(path, gene_map_dict):
    #merge all gene maps given a directory
    #the gene length stays constant.  The number of reads do not.
    path_contents = os.listdir(path)
    #gene_map_dict = dict()
    for content in path_contents:
        item = os.path.join(path, content)
        
        if(item.endswith("_gene_map.tsv")):
            #print("looking at:", item)
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
                        new_reads = old_reads + reads
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
                 
                
                        
                        
                        #print("new entry:", gene_name, gene_map_dict[gene_name])
    #return gene_map_dict
def final_gene_map_merge(gene_map_list, final_gene_map):
    #merge all the interim gene maps that we've got
    #final_gene_map = dict()
    for gene_map in gene_map_list:
        for gene_name_key in gene_map:
            if(gene_name_key in final_gene_map):
                old_entry = final_gene_map[gene_name_key]
                #print("old entry:", old_entry)
                old_gene_length = int(old_entry[0])
                old_reads = old_entry[2:]
                added_reads = gene_map[gene_name_key][2:]
                added_gene_length = int(gene_map[gene_name_key][0])
                
                final_reads = old_reads + added_reads
                new_number_of_reads = len(final_reads)
                final_gene_length = old_gene_length
                if(old_gene_length != added_gene_length):
                    print(dt.today(), "final merge warning: gene lengths of the same key don't match")
                    print("on:", gene_name_key, "this shouldn't happen")
                    if(old_gene_length > added_gene_length):
                        final_gene_length = old_gene_length
                    else:
                        final_gene_length = added_gene_length
                
                final_gene_map[gene_name_key] = [final_gene_length, new_number_of_reads] + final_reads
            else:
                final_gene_map[gene_name_key] = gene_map[gene_name_key]
    #return final_gene_map


            
            
    
def merge_fastas(path_0, path_1, section, header, extension):
    #it's just a simple concatenation job.
    path_contents = os.listdir(path_0)
    final_fasta_file = os.path.join(path_1, header + "_" + section + extension)
    with open(final_fasta_file, "w") as final_fasta:
        for content in path_contents:
            item = os.path.join(path_0, content)
            #print("fasta merge:", item)
            if(content.startswith(section)) and (content.endswith(extension)):
                with open(item, "r") as sample_fasta:
                    for line in sample_fasta:
                        final_fasta.write(line)
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
        
        
def merge_all_proteins(path, gene_transcripts_list, export_path):
    all_proteins_path = os.path.join(export_path, "all_proteins.faa")
    with open(all_proteins_path, "w") as proteins_out:
        for item in os.listdir(path):
            
            if(item.endswith(".faa")):
                full_path = os.path.join(os.path.abspath(path), item)
                with open(full_path, "r") as dmd_proteins:
                    for line in dmd_proteins:
                        proteins_out.write(line)
        #for gene in gene_transcripts_list:
        #    out_line = str(gene.id) + "\n" + str(gene.seq) + "\n"
        #    proteins_out.write(out_line)
        SeqIO.write(gene_transcripts_list, proteins_out, "fasta")
        
def handle_final_proteins(final_path, export_path):
    #just a function dump so we can do this step in parallel with other things

    #scrub for duplicates (because seqIO index doesn't like duplicates)
    bwa_blat_proteins_file = merge_fastas(final_path, final_path, "all", "BWA_BLAT_proteins", ".fna")
    unique_fna_file = "unique_" + os.path.basename(bwa_blat_proteins_file)
    scrub_duplicates(bwa_blat_proteins_file, unique_fna_file)
    
    
    gene_transcripts_list = convert_genes_to_proteins(unique_fna_file)
    #print("GENE TRANSCRIPT:", len(gene_transcripts_list))
    #for item in gene_transcripts_list:
    #    print(item)
    
    merge_all_proteins(final_path, gene_transcripts_list, export_path)
    
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
                    

 
if __name__ == "__main__":
    assemble_path           = sys.argv[1]
    bwa_path                = sys.argv[2]
    blat_path               = sys.argv[3]
    diamond_path            = sys.argv[4]
    final_path              = sys.argv[5]
    export_path             = sys.argv[6]
    #diamond_proteins_file   = sys.argv[5]
    
    gene_report_path = os.path.join(assemble_path, "data", "1_mgm", "gene_report.txt")
    contig_map_path = os.path.join(bwa_path, "contig_map.tsv")
    
    gene_segment_dict = import_gene_report(gene_report_path)
    contig_map_dict = import_contig_map(contig_map_path)
    translate_gene_segement_map(gene_segment_dict, contig_map_dict)
    
    for item in gene_segment_dict:
        print(item, gene_segment_dict[item])
    #sys.exit("Early")
    manager = mp.Manager()
    mgr_bwa_gene_map            = manager.dict()
    mgr_blat_gene_map           = manager.dict()
    mgr_diamond_gene_map        = manager.dict()
    mgr_gene_transcripts_dict   = manager.dict()
    mgr_final_gene_map          = manager.dict()
    
    
    process_store = []
    
    
    #merge the leftover reads, and annotated genes (for protein translation)
    make_merge_fasta_process(process_store, diamond_path, export_path, "GA_leftover", ".fasta")  #leftover reads
    make_merge_fasta_process(process_store, blat_path, final_path, "BLAT_annotated", ".fna")    #BLAT genes
    make_merge_fasta_process(process_store, bwa_path, final_path, "BWA_annotated", ".fna")      #BWA genes
    make_merge_fasta_process(process_store, diamond_path, final_path, "dmd", ".faa")            #DIAMOND proteins
    
    #merge the gene maps
    bwa_gene_map_merge_process = mp.Process(target = concatenate_gene_maps, args = (bwa_path, mgr_bwa_gene_map))
    bwa_gene_map_merge_process.start()
    process_store.append(bwa_gene_map_merge_process)
    print(dt.today(), "launched BWA gene map merge")
    blat_gene_map_merge_process = mp.Process(target = concatenate_gene_maps, args = (blat_path, mgr_blat_gene_map))
    blat_gene_map_merge_process.start()
    process_store.append(blat_gene_map_merge_process)
    print(dt.today(), "launched BLAT gene map merge")
    diamond_gene_map_merge_process = mp.Process(target = concatenate_gene_maps, args = (diamond_path, mgr_diamond_gene_map))
    diamond_gene_map_merge_process.start()
    process_store.append(diamond_gene_map_merge_process)
    print(dt.today(), "launched DIAMOND gene map merge")
    
    
    for item in process_store:
        item.join()
    process_store[:] = []
    
    
    
    
    #secondary combine on the genes (BWA and BLAT)
    make_second_merge_process(process_store, final_path, final_path, "all", ".fna")
    make_merge_fasta_process(process_store, final_path, final_path, "dmd", ".faa")
    
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
    

        