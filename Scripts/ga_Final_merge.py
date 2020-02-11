#This code performs all of the merging needed for GA to be completed.

import bio
import os
import sys


def convert_genes_to_proteins(mapped_gene_file, out_prot):
    # WRITE OUTPUT: BWA&BLAT&DMD-aligned gene/protIDs and aa seqs.  It's for a downstream tool.     
    # (.faa; fasta-format):
    # convert previously mapped genes to the AA format.
    gene_seqs = SeqIO.index(mapped_gene_file,"fasta")           # key=geneID, value=SeqRecord
    genes_trans= []
    for gene in gene_seqs:                                  # Take each BWA&BLAT-aligned genes
        try:
            genes_trans.append(SeqRecord(seq= gene_seqs[gene].seq.translate(stop_symbol=""), id= gene_seqs[gene].id, description= gene_seqs[gene].description))
                                                            #  and translate its SeqRecord sequence to aa.
        except:
            pass
    #Then merge the DIAMOND proteins to the same file.
    print(dt.today(), "writing fasta")
    with open(prot_file,"w") as out_prot:
        SeqIO.write(genes_trans, out_prot, "fasta")         # Write aligned gene aa seqs
        SeqIO.write(proteins, out_prot, "fasta")            # and aligned proteins aa seqs.