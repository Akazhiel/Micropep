#!/usr/bin/env python

import os
import argparse
import numpy as np
import pandas as pd
from Bio import Seq
from Bio import SeqIO
import pyranges as pr


def obtain_fasta_chomosomes(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    chromosomes = []
    for record in records:
        chromosomes.append(record.id)
    return chromosomes


def filter_biotype(fasta_file, annotation_gtf, input_file, biotype):
    chromosomes = obtain_fasta_chomosomes(fasta_file)

    # Intersect chromosomes in fasta file with those in annotation gtf

    gtf_df = pr.read_gtf(annotation_gtf)
    gtf_df = gtf_df.subset(lambda df: df.Chromosome.isin(chromosomes))
    gtf_df = gtf_df[
        [
            "Chromosome",
            "Strand",
            "Start",
            "End",
            "gene_biotype",
            "Feature",
            "gene_id",
        ]
    ]
    gtf_df = gtf_df.df
    gtf_df["gene_biotype"] = gtf_df["gene_biotype"].ffill()
    gtf_df = pr.PyRanges(gtf_df)
    gtf_df = gtf_df[gtf_df.Feature == "exon"]
    # gtf_df = gtf_df.subset(lambda df: df.gene_biotype.isin(biotype))

    # Intersect chromosomes in fasta file with those in pipeline derived gtf

    input_df = pr.read_gtf(input_file)
    input_df = input_df.subset(lambda df: df.Chromosome.isin(chromosomes))
    input_df = input_df.drop(["Source", "Score", "Frame", "transcript_id"])

    return gtf_df, input_df


def filter_by_overlap(
    fasta_file, annotation_gtf, input_file, overlap, biotype
):
    gtf_df, input_df = filter_biotype(
        fasta_file, annotation_gtf, input_file, biotype
    )
    # result_df = input_df.coverage(gtf_df)
    # result_df = result_df[result_df.FractionOverlaps < overlap]
    # result_df = result_df.drop(["NumberOverlaps", "FractionOverlaps"])
    return input_df


def main(
    fasta_file, annotation_gtf, input_file, overlap, biotype, output, aa_length
):
    file_name = "Micropeptide_data.csv"
    protein_df = filter_by_overlap(
        fasta_file, annotation_gtf, input_file, overlap, biotype
    )
    protein_df.Chromosome = (
        protein_df.Chromosome.cat.remove_unused_categories()
    )
    protein_df.Fasta_seq = pr.get_fasta(protein_df, fasta_file)
    peptide = []

    for seq in protein_df.Fasta_seq:
        nuc = str(seq)
        aa_seqs = [str(Seq.translate(nuc[i:], to_stop=True)) for i in range(3)]
        peptide.append([str(seqs) for seqs in aa_seqs])

    start_codons = ["L", "M", "V", "I"]

    final_df = protein_df.df
    tmp_ = final_df
    tmp_["Sequence"] = peptide
    tmp_ = tmp_.explode("Sequence")

    seqs = []
    possible_seqs = []

    for seq in tmp_["Sequence"]:
        for x in start_codons:
            seq = str(seq)
            seq = seq[seq.find(x) :]
            if len(seq) == 1 or len(seq) == 0:
                seq = ""
            else:
                seq = seq
            seqs.append(seq)

    possible_seqs = [seqs[n : n + 12] for n in range(0, len(seqs), 12)]
    final_df["Sequence"] = possible_seqs
    final_df = final_df.explode("Sequence")
    final_df = final_df[final_df["Sequence"].str.strip().astype(bool)]
    results = (
        final_df[final_df["Sequence"].str.len() <= aa_length]
        .copy()
        .reset_index(drop=True)
    )
    results["Protein_ID"] = "Protein_" + pd.Series(
        np.arange(1, len(results) + 1, 1)
    ).astype(str)
    if not output:
        results.to_csv(
            os.path.join(os.getcwd(), file_name),
            header=True,
            index=False,
            sep="\t",
        )
    elif os.path.isdir(output):
        results.to_csv(
            os.path.join(output, file_name),
            header=True,
            index=False,
            sep="\t",
        )
    else:
        results.to_csv(output, header=True, index=False, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fasta_file",
        "-ff",
        type=str,
        required=True,
        help="Provide the genome fasta file to extract the nucleotide sequence from.",
    )

    parser.add_argument(
        "--annotation_gtf",
        "-agf",
        type=str,
        required=False,
        help="Provide the annotation GTF file.",
    )

    parser.add_argument(
        "--input_file",
        "-i",
        type=str,
        required=False,
        help="Provide the GTF file obtained from running the filter_samples script.",
    )

    parser.add_argument(
        "--overlap",
        type=float,
        required=False,
        default=0,
        help="Provide the percentage of overlap to filter the regions with.",
    )

    parser.add_argument(
        "--biotype",
        "-bt",
        type=str,
        # default=["protein_coding"],
        required=False,
        nargs="+",
        help="List of Ensembl gene types to keep (E.x protein_coding lincRNA). Defaults to protein_coding.",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=False,
        help="Output path to save the filtered dataset file. Default format is csv with tabular separation.",
    )

    parser.add_argument(
        "--length",
        "-l",
        type=int,
        required=False,
        default=100,
        help="Peptides with more length than the specified value will be filtered out.",
    )
    args = parser.parse_args()

    main(
        args.fasta_file,
        args.annotation_gtf,
        args.input_file,
        args.overlap,
        args.biotype,
        args.output,
        args.length,
    )
