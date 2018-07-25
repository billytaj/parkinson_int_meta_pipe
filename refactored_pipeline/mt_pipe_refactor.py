#!/usr/bin/env python

import sys
import os
import os.path
from argparse import ArgumentParser
import multiprocessing as mp
import mt_pipe_commands as mpcom
import mt_pipe_paths as mpp
import time


def make_folder(folder_path):
    if not (os.path.exists(folder_path)):
        os.makedirs(folder_path)


# Used to determine quality encoding of fastq sequences.
# Assumes Phred+64 unless there is a character within the first 10000 reads with encoding in the Phred+33 range.
def determine_encoding(fastq):
    encoding = 64
    with open(fastq) as fq:
        line_count = 0
        encoding_count = 0
        for line in fq:
            line_count += 1
            if line_count % 4 == 0:
                encoding_count += 1
                for char in line:
                    if ord(char) < 64:
                        encoding = 33
                        break
                if encoding_count == 10000 or encoding == 33:
                    break

    return encoding


# handles where to kill the pipeline, due to the prev step behaving badly
def check_where_kill(dep_job_label=None, dep_path=None):
    if dep_job_label is None:
        if dep_path is None:
            return True
        else:
            dep_job_path = dep_path
    else:
        dep_job_path = os.path.join(dep_job_label, "data", "final_results")

    file_list = os.listdir(dep_job_path)
    if len(file_list) > 0:
        for item in file_list:
            file_check_path = os.path.join(dep_job_path, item)
            if (os.path.getsize(file_check_path)) == 0:
                print("empty file detected: rerunning stage")
                sys.exit("bad dep")
        # run the job, silently
        return True
    else:
        print("stopping the pipeline.  dependencies don't exist")
        sys.exit("no dep")


# handles where to auto-resume the pipeline on a subsequent run
def check_where_resume(job_label=None, full_path=None, dep_job_path=None):
    check_where_kill(dep_job_path)
    if job_label:
        job_path = os.path.join(job_label, "data", "final_results")
    else:
        job_path = full_path

    print("looking at:", job_path)

    if os.path.exists(job_path):
        file_list = os.listdir(job_path)
        if len(file_list) > 0:
            for item in file_list:
                file_check_path = os.path.join(job_path, item)
                if (os.path.getsize(file_check_path)) == 0:
                    print("empty file detected: rerunning stage")
                    return False
            print("bypassing!")
            return True
        else:
            print("running")
            return False
    else:
        print("doesn't exist: running")
        return False


def main(pair_1_path, pair_2_path, single_path, output_folder_path, threads):
    if single_path:
        operating_mode = "single"
        quality_encoding = determine_encoding(single_path)
        print("OPERATING IN SINGLE-ENDED MODE")
    else:
        operating_mode = "paired"
        quality_encoding = determine_encoding(pair_1_path)
        print("OPERATING IN PAIRED-MODE")
    if threads == 0:
        thread_count = mp.cpu_count()
    else:
        thread_count = threads
    mp_store = []  # stores the multiprocessing processes

    # --------------------------------------------------
    # profiling vars

    start_time = time.time()

    if operating_mode == "paired":
        preprocess_label = "preprocess"
        rRNA_filter_label = "rRNA_filter"
        repop_job_label = "duplicate_repopulation"
        assemble_contigs_label = "assemble_contigs"
        gene_annotation_BWA_label = "gene_annotation_BWA"
        gene_annotation_BLAT_label = "gene_annotation_BLAT"
        gene_annotation_DIAMOND_label = "gene_annotation_DIAMOND"
        taxon_annotation_label = "taxonomic_annotation"
        ec_annotation_label = "enzyme_annotation"
        network_label = "RPKM_network"
        visualization_label = "visualization"

        # Creates our command object, for creating shellscripts.
        comm = mpcom.mt_pipe_commands(Quality_score=quality_encoding, Thread_count=thread_count,
                                      sequence_path_1=pair_1_path, sequence_path_2=pair_2_path,
                                      sequence_signle=single_path)  # start obj
        paths = mpp.tool_path_obj()

        # This is the format we use to launch each stage of the pipeline.
        # We start a multiprocess that starts a subprocess.
        # The subprocess is created from the comm object

        # The preprocess stage
        preprocess_start = time.time()
        preprocess_path = os.path.join(output_folder_path, preprocess_label)
        if not check_where_resume():
            process = mp.Process(
                target=comm.create_pbs_and_launch,
                args=(
                    preprocess_label,
                    comm.create_pre_double_command(preprocess_label),
                    True
                )
            )
            process.start()  # start the multiprocess
            process.join()  # wait for it to end
        preprocess_end = time.time()

        # rRNA removal stage
        rRNA_filter_start = time.time()

        rRNA_filter_path = os.path.join(output_folder_path, rRNA_filter_label)
        rRNA_filter_orphans_fastq_folder = os.path.join(output_folder_path, "rRNA_filter", "data", "orphans",
                                                        "orphans_fastq")
        rRNA_filter_pair_1_fastq_folder = os.path.join(output_folder_path, "rRNA_filter", "data", "pair_1", "pair_1_fastq")
        rRNA_filter_pair_2_fastq_folder = os.path.join(output_folder_path, "rRNA_filter", "data", "pair_2", "pair_2_fastq")

        if not check_where_resume(None, preprocess_path):
            process = mp.Process(
                target=comm.create_pbs_and_launch,
                args=(
                    rRNA_filter_label,
                    comm.create_rRNA_filter_prep_command(
                        rRNA_filter_label, int(mp.cpu_count() / 2), preprocess_label),
                    True
                )
            )
            process.start()
            process.join()

            orphans_mRNA_path = os.path.join(rRNA_filter_path, "data", "orphans", "orphans_mRNA")
            if not check_where_resume(orphans_mRNA_path, None):
                for item in os.listdir(rRNA_filter_orphans_fastq_folder):
                    file_root_name = os.path.splitext(item)[0]
                    inner_name = file_root_name + "_infernal"
                    process = mp.Process(
                        target=comm.create_pbs_and_launch,
                        args=(
                            "rRNA_filter",
                            comm.create_rRNA_filter_command("rRNA_filter", "orphans", file_root_name),
                            True,
                            inner_name
                        )
                    )
                    process.start()
                    mp_store.append(process)  # pack all the processes into a list
            for item in mp_store:
                item.join()  # wait for things to finish
            mp_store[:] = []  # clear the list

            pair_1_mRNA_path = os.path.join(rRNA_filter_path, "data", "pair_1", "pair_1_mRNA")
            if not check_where_resume(pair_1_mRNA_path, None):
                for item in os.listdir(rRNA_filter_pair_1_fastq_folder):
                    file_root_name = os.path.splitext(item)[0]
                    inner_name = file_root_name + "_infernal"
                    process = mp.Process(
                        target=comm.create_pbs_and_launch,
                        args=(
                            "rRNA_filter",
                            comm.create_rRNA_filter_command("rRNA_filter", "pair_1", file_root_name),
                            True,
                            inner_name
                        )
                    )
                    process.start()
                    mp_store.append(process)
            for item in mp_store:
                item.join()  # wait for things to finish
            mp_store[:] = []  # clear the list

            pair_2_mRNA_path = os.path.join(rRNA_filter_path, "data", "pair_2", "pair_2_mRNA")
            if not check_where_resume(pair_2_mRNA_path, None):
                for item in os.listdir(rRNA_filter_pair_2_fastq_folder):
                    file_root_name = os.path.splitext(item)[0]
                    inner_name = file_root_name + "_infernal"
                    process = mp.Process(
                        target=comm.create_pbs_and_launch,
                        args=(
                            "rRNA_filter",
                            comm.create_rRNA_filter_command("rRNA_filter", "pair_2", file_root_name),
                            True,
                            inner_name
                        )
                    )
                    process.start()
                    mp_store.append(process)

            for item in mp_store:
                item.join()  # wait for things to finish
            mp_store[:] = []  # clear the list

            inner_name = "rRNA_filter_post"
            process = mp.Process(
                target=comm.create_pbs_and_launch,
                args=(
                    rRNA_filter_label,
                    comm.create_rRNA_filter_post_command(rRNA_filter_label),
                    True,
                    inner_name
                )
            )
            process.start()
            process.join()

        rRNA_filter_end = time.time()
        # -------------------------------------------------------------
        # Duplicate repopulation
        repop_start = time.time()
        repop_job_path = os.path.join(output_folder_path, repop_job_label)
        if not check_where_resume(None, rRNA_filter_path):
            process = mp.Process(
                target=comm.create_pbs_and_launch,
                args=(
                    repop_job_label,
                    comm.create_repop_command(repop_job_label, preprocess_label, rRNA_filter_label),
                    True
                )
            )
            process.start()
            process.join()
        repop_end = time.time()
        # ----------------------------------------
        # Assemble contigs
        assemble_contigs_start = time.time()
        assemble_contigs_path = os.path.join(output_folder_path, assemble_contigs_label)
        if not check_where_resume(None, repop_job_path):
            process = mp.Process(
                target=comm.create_pbs_and_launch,
                args=(
                    assemble_contigs_label,
                    comm.create_assemble_contigs_command(assemble_contigs_label, repop_job_label),
                    True
                )
            )
            process.start()
            process.join()
        assemble_contigs_end = time.time()

        # ----------------------------------------------
        # BWA gene annotation
        GA_BWA_start = time.time()
        gene_annotation_BWA_path = os.path.join(output_folder_path, gene_annotation_BWA_label)
        if not check_where_resume(None, assemble_contigs_path):

            names = ["contigs", "orphans", "pair_1", "pair_2"]
            mp_store[:] = []
            for item in names:
                inner_name = "BWA_" + item
                process = mp.Process(
                    target=comm.create_pbs_and_launch,
                    args=(
                        gene_annotation_BWA_label,
                        comm.create_BWA_annotate_command(gene_annotation_BWA_label, assemble_contigs_label, item),
                        True,
                        inner_name
                    )
                )
                process.start()
                mp_store.append(process)

            for item in mp_store:
                item.join()
            mp_store[:] = []  # clear the list

            inner_name = "BWA_pp"
            process = mp.Process(
                target=comm.create_pbs_and_launch,
                args=(
                    gene_annotation_BWA_label,
                    comm.create_BWA_pp_command(gene_annotation_BWA_label, assemble_contigs_label),
                    True,
                    inner_name
                )
            )
            process.start()
            process.join()
        GA_BWA_end = time.time()

        # ------------------------------------------------
        # BLAT gene annotation
        GA_BLAT_start = time.time()
        gene_annotation_BLAT_path = os.path.join(output_folder_path, gene_annotation_BLAT_label)
        if not check_where_resume(None, gene_annotation_BWA_path):

            BlatPool = mp.Pool(int(thread_count / 2))
            sections = ["contigs", "orphans", "pair_1", "pair_2"]
            for section in sections:
                for fasta_db in os.listdir(paths.DNA_DB_Split):
                    if fasta_db.endswith(".fasta") or fasta_db.endswith(".ffn") or fasta_db.endswith(".fsa") or fasta_db.endswith(".fas") or fasta_db.endswith(".fna") or fasta_db.endswith(".ffn"):
                        inner_name = "BLAT_" + section + "_" + fasta_db
                        BlatPool.apply_async(comm.create_pbs_and_launch,
                                             args=(
                                                 gene_annotation_BLAT_label,
                                                 comm.create_BLAT_annotate_command(gene_annotation_BLAT_label,
                                                                                   gene_annotation_BWA_label, section,
                                                                                   fasta_db),
                                                 True,
                                                 inner_name
                                             )
                                             )
            BlatPool.close()
            BlatPool.join()

            for section in sections:
                inner_name = section + "_cat"
                process = mp.Process(
                    target=comm.create_pbs_and_launch,
                    args=(
                        gene_annotation_BLAT_label,
                        comm.create_BLAT_cat_command(gene_annotation_BLAT_label, section),
                        True,
                        inner_name
                    )
                )
                process.start()
                mp_store.append(process)
            for item in mp_store:
                item.join()
            mp_store[:] = []

            inner_name = "BLAT_pp"
            process = mp.Process(
                target=comm.create_pbs_and_launch,
                args=(
                    gene_annotation_BLAT_label,
                    comm.create_BLAT_pp_command(gene_annotation_BLAT_label, gene_annotation_BWA_label),
                    True,
                    inner_name
                )
            )
            process.start()
            process.join()
        GA_BLAT_end = time.time()

        # ------------------------------------------------------
        # Diamond gene annotation
        GA_DIAMOND_start = time.time()
        gene_annotation_DIAMOND_path = os.path.join(output_folder_path, gene_annotation_DIAMOND_label)
        if not check_where_resume(None, gene_annotation_BLAT_path):

            names = ["contigs", "orphans", "pair_1", "pair_2"]
            for item in names:
                inner_name = item + "_run_diamond"
                process = mp.Process(
                    target=comm.create_pbs_and_launch,
                    args=(
                        gene_annotation_DIAMOND_label,
                        comm.create_DIAMOND_annotate_command(gene_annotation_DIAMOND_label, gene_annotation_BLAT_label,
                                                             item),
                        True,
                        inner_name
                    )
                )
                process.start()
                mp_store.append(process)
            for item in mp_store:
                item.join()
            mp_store[:] = []

            inner_name = "diamond_pp"
            process = mp.Process(
                target=comm.create_pbs_and_launch,
                args=(
                    gene_annotation_DIAMOND_label,
                    comm.create_DIAMOND_pp_command(gene_annotation_DIAMOND_label, gene_annotation_BLAT_label),
                    True,
                    inner_name
                )
            )
            process.start()
            process.join()
        GA_DIAMOND_end = time.time()
        # ------------------------------------------------------
        # Taxonomic annotation
        TA_start = time.time()
        if not check_where_resume(None, gene_annotation_DIAMOND_path):
            process = mp.Process(
                target=comm.create_pbs_and_launch,
                args=(
                    taxon_annotation_label,
                    comm.create_taxonomic_annotation_command(taxon_annotation_label, assemble_contigs_label,
                                                             gene_annotation_DIAMOND_label),
                    True
                )
            )
            process.start()
            process.join()
        TA_end = time.time()

        # ------------------------------------------------------
        # Detect EC annotation
        EC_start = time.time()
        EC_DETECT_start = time.time()
        ec_annotation_path = os.path.join(output_folder_path, ec_annotation_label)
        if not check_where_resume(None, gene_annotation_DIAMOND_path):
            # Preparing folders for DETECT
            process = mp.Process(
                target=comm.create_pbs_and_launch,
                args=(
                    ec_annotation_label,
                    comm.create_EC_DETECT_prep(ec_annotation_label, gene_annotation_DIAMOND_label,
                                               int(mp.cpu_count() / 2)),
                    True
                )
            )
            process.start()
            process.join()

            # Running DETECT on split protein files
            proteins_path = os.path.join(output_folder_path, ec_annotation_label, "data", "0_proteins")
            for item in os.listdir(proteins_path):
                file_root_name = os.path.splitext(item)[0]
                inner_name = file_root_name + "_detect"
                process = mp.Process(
                    target=comm.create_pbs_and_launch,
                    args=(
                        ec_annotation_label,
                        comm.create_EC_DETECT_command(ec_annotation_label, file_root_name),
                        True,
                        inner_name
                    )
                )
                process.start()
                mp_store.append(process)  # pack all the processes into a list

            for item in mp_store:
                item.join()  # wait for things to finish
            mp_store[:] = []  # clear the list
        EC_DETECT_end = time.time()

        # --------------------------------------------------------------
        # Priam and Diamond EC annotation
        EC_PRIAM_DIAMOND_start = time.time()
        if not check_where_resume(None, gene_annotation_DIAMOND_path):
            process = mp.Process(
                target=comm.create_pbs_and_launch,
                args=(
                    ec_annotation_label,
                    comm.create_EC_PRIAM_DIAMOND_command(ec_annotation_label, gene_annotation_DIAMOND_label),
                    True
                )
            )
            process.start()
            process.join()

            inner_name = "ea_post"
            process = mp.Process(
                target=comm.create_pbs_and_launch,
                args=(
                    ec_annotation_label,
                    comm.create_EC_postprocess_command(ec_annotation_label, gene_annotation_DIAMOND_label),
                    True,
                    inner_name
                )
            )
            process.start()
            process.join()
        EC_PRIAM_DIAMOND_end = time.time()
        EC_end = time.time()

        # ------------------------------------------------------
        # RPKM Table and Cytoscape Network
        Cytoscape_start = time.time()
        network_path = os.path.join(output_folder_path, network_label)
        if not check_where_resume(None, ec_annotation_path):
            process = mp.Process(
                target=comm.create_pbs_and_launch,
                args=(
                    network_label,
                    comm.create_Network_generation_command(network_label, gene_annotation_DIAMOND_label,
                                                           taxon_annotation_label, ec_annotation_label),
                    True
                )
            )
            process.start()
            process.join()
        Cytoscape_end = time.time()

        # ------------------------------------------------------
        # Final Pie Charts
        Chart_start = time.time()
        if not check_where_resume(None, network_path):
            process = mp.Process(
                target=comm.create_pbs_and_launch,
                args=(
                    visualization_label,
                    comm.create_visualization_command(visualization_label, network_label),
                    True
                )
            )
            process.start()
            process.join()
        Chart_end = time.time()

        end_time = time.time()
        print("Total runtime:", '%1.1f' % (end_time - start_time), "s", "start:", '%1.1f' % start_time, "end:",
              '%1.1f' % end_time)
        print("preprocess:", '%1.1f' % (preprocess_end - preprocess_start), "s", "start:", '%1.1f' % preprocess_start,
              "end:", '%1.1f' % preprocess_end)
        print("rRNA filter:", '%1.1f' % (rRNA_filter_end - rRNA_filter_start), "s", "start:",
              '%1.1f' % rRNA_filter_start, "end:", '%1.1f' % rRNA_filter_end)
        print("repop:", '%1.1f' % (repop_end - repop_start), "s", "start:", '%1.1f' % repop_start, "end:",
              '%1.1f' % repop_end)
        print("assemble contigs:", '%1.1f' % (assemble_contigs_end - assemble_contigs_start), "s", "start:",
              '%1.1f' % assemble_contigs_start, "end:", '%1.1f' % assemble_contigs_end)
        print("GA BWA:", '%1.1f' % (GA_BWA_end - GA_BWA_start), "s", "start:", '%1.1f' % GA_BWA_start, "end:",
              '%1.1f' % GA_BWA_end)
        print("GA BLAT:", '%1.1f' % (GA_BLAT_end - GA_BLAT_start), "s", "start:", '%1.1f' % GA_BLAT_start, "end:",
              '%1.1f' % GA_BLAT_end)
        print("GA DIAMOND:", '%1.1f' % (GA_DIAMOND_end - GA_DIAMOND_start), "s", "start:", '%1.1f' % GA_DIAMOND_start,
              "end:", '%1.1f' % GA_DIAMOND_end)
        print("TA:", '%1.1f' % (TA_end - TA_start), "s", "start:", '%1.1f' % TA_start, "end:", '%1.1f' % TA_end)
        print("EC:", '%1.1f' % (EC_end - EC_start), "s", "start:", '%1.1f' % EC_start, "end:", '%1.1f' % EC_end)
        print("EC DETECT:", '%1.1f' % (EC_DETECT_end - EC_DETECT_start), "s", "start:", '%1.1f' % EC_DETECT_start,
              "end:", '%1.1f' % EC_DETECT_end)
        print("EC PRIAM + DIAMOND:", '%1.1f' % (EC_PRIAM_DIAMOND_end - EC_PRIAM_DIAMOND_start), "s", "start:",
              '%1.1f' % EC_PRIAM_DIAMOND_start, "end:", '%1.1f' % EC_PRIAM_DIAMOND_end)
        print("Cytoscape:", '%1.1f' % (Cytoscape_end - Cytoscape_start), "s", "start:", '%1.1f' % Cytoscape_start,
              "end:", '%1.1f' % Cytoscape_end)
        print("Charts: ", '%1.1f' % (Chart_end - Chart_start), "s", "start:", '%1.1f' % Chart_start, "end:",
              '%1.1f' % Chart_end)

    elif operating_mode == "single":
        print("not ready")


if __name__ == "__main__":
    # This is where the code starts
    # There's a few operating modes, mainly "docker", and "singularity".  These modes edit the pipeline filepaths

    parser = ArgumentParser(description="MetaPro - Meta-omic sequence processing and analysis pipeline"
                                        "Version 1.0 © 2018")

    parser.add_argument("-c", "--config", type=str,
                        help="Path to the configureation file")
    parser.add_argument("-1", "--pair1", type=str,
                        help="Path to the file containing the forward paired-end reads in fastq format")
    parser.add_argument("-2", "--pair2", type=str,
                        help="Path to the file containing the reverse paired-end reads in fastq format")
    parser.add_argument("-s", "--single", type=str,
                        help="Path to the file containing the single-end reads in fastq format")
    parser.add_argument("-o", "--output_folder", type=str, required=True,
                        help="Path of the folder for the output of the pipeline")
    parser.add_argument("-t", "--num_threads", type=int,
                        help="Maximum number of threads used by the pipeline")

    args = parser.parse_args()

    if (args.pair1 and not args.pair2) or (args.pair2 and not args.pair1):
        print("You must specify both forward and reverse reads for a paired-end run")
        sys.exit()
    elif args.single and (args.pair1 or args.pair2):
        print("You cannot specify both paired-end and single-end reads in a single run.")
        sys.exit()

    pair_1 = args.pair1 if args.pair1 else ""
    pair_2 = args.pair2 if args.pair2 else ""
    single = args.single if args.single else ""
    output_folder = args.output_folder
    num_threads = args.num_threads if args.num_threads else 0

    if not (os.path.exists(output_folder)):
        print("output folder does not exist.  Now building directory.")
        os.makedirs(output_folder)
    os.chdir(output_folder)
    main(pair_1, pair_2, single, output_folder, num_threads)
