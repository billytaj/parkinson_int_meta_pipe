#This code performs all of the merging needed for GA to be completed.

import os
import sys
import multiprocessing as mp
from datetime import datetime as dt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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
    #Then merge the DIAMOND proteins to the same file.
    #print(dt.today(), "writing fasta")
    return gene_trans
    #gene_trans_dict[section] = gene_trans
    #with open(prot_file,"a") as out_prot:
    #    SeqIO.write(genes_trans, out_prot, "fasta")         # Write aligned gene aa seqs
        #SeqIO.write(proteins, out_prot, "fasta")            # and aligned proteins aa seqs.
        
def concatenate_gene_maps(path, gene_map_dict):
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
    
    
def export_gene_map(gene_map, path):
    final_gene_map_file = os.path.join(path, "gene_map.tsv")
    
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
if __name__ == "__main__":
    bwa_path                = sys.argv[1]
    blat_path               = sys.argv[2]
    diamond_path            = sys.argv[3]
    final_path              = sys.argv[4]
    #diamond_proteins_file   = sys.argv[5]
    
    manager = mp.Manager()
    mgr_bwa_gene_map            = manager.dict()
    mgr_blat_gene_map           = manager.dict()
    mgr_diamond_gene_map        = manager.dict()
    mgr_gene_transcripts_dict   = manager.dict()
    mgr_final_gene_map          = manager.dict()
    
    
    process_store = []
    
    
    
    make_merge_fasta_process(process_store, diamond_path, final_path, "GA_leftover", ".fasta")
    make_merge_fasta_process(process_store, blat_path, final_path, "BLAT_annotated", ".fna")
    make_merge_fasta_process(process_store, bwa_path, final_path, "BWA_annotated", ".fna")
    make_merge_fasta_process(process_store, diamond_path, final_path, "dmd", ".faa")
    for item in process_store:
        item.join()
    process_store[:] = []
    
    
    #secondary combine
    make_second_merge_process(process_store, final_path, final_path, "all", ".fna")
    #make_merge_fasta_process(process_store, final_path, final_path, "dmd", ".faa")
    for item in process_store:
        item.join()
    process_store[:] = []
    #make_convert_proteins_process(process_store, final_path, mgr_gene_transcripts_dict)
    #
    #gene_trans_dict = dict(mgr_gene_transcripts_dict)
    #for item in gene_trans_dict:
    #    print(item)
    
    
    #tertiary combine
    bwa_blat_proteins_file = merge_fastas(final_path, final_path, "all", "BWA_BLAT_proteins", ".fna")
    unique_fna_file = "unique_" + os.path.basename(bwa_blat_proteins_file)
    scrub_duplicates(bwa_blat_proteins_file, unique_fna_file)
    
    
    gene_transcripts_list = convert_genes_to_proteins(unique_fna_file)
    print("GENE TRANSCRIPT:", len(gene_transcripts_list))
    #for item in gene_transcripts_list:
    #    print(item)
    
    merge_all_proteins(final_path, gene_transcripts_list, final_path)
    
    sys.exit("premature death")
    
    
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
    #------------------------------------------------------------------------------------
    #convert genes to proteins
    #convert_genes_to_proteins(section, mapped_gene_file, gene_trans_dict)
    
    
    
    
    
    
    
    
    
    
    
    ##export the transcripts
    #gene_trans_dict = dict(mgr_gene_transcripts_dict)
    #write_transcripts_process = mp.Process(target = export_proteins, args = (diamond_proteins_file, gene_trans_dict))
    #write_transcripts_process.start()
    #process_store.append(write_transcripts_process)
    
    
    bwa_gene_map = dict(mgr_bwa_gene_map)
    blat_gene_map = dict(mgr_blat_gene_map)
    diamond_gene_map = dict(mgr_diamond_gene_map)
    
    gene_map_list = [bwa_gene_map, blat_gene_map, diamond_gene_map]
    #merge all gene maps
    final_gene_map_process = mp.Process(target = final_gene_map_merge, args = (gene_map_list, mgr_final_gene_map))
    final_gene_map_process.start()
    process_store.append(final_gene_map_process)
    
    for item in process_store:
        item.join()
    process_store[:] = []
    
    final_gene_map = dict(mgr_final_gene_map)
    
    export_gene_map(final_gene_map, diamond_path)
    
    #for item in final_gene_map:
    #    print(item, final_gene_map[item])
        