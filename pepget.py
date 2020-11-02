#!/usr/bin/env python

import os
import argparse
import subprocess
import numpy as np
import pandas as pd
from Bio import Seq
import pyranges as pr
from pathlib import Path
from collections import defaultdict


def main(fasta_file, codon_list, input_file, output, aa_length):
    results_name = (
        output
        if output and not os.path.isfile(output)
        else "Micropeptide_data.csv"
    )

    if Path(fasta_file).suffix == "*.gz":
        print(
            "Unable to process compressed fasta. Please uncompress before starting."
        )
    else:
        chromosomes = (
            subprocess.run(
                f"grep '^>' {fasta_file} | cut -c 2- | sort -V | uniq",
                shell=True,
                capture_output=True,
            )
            .stdout.decode("utf-8")
            .strip()
            .split("\n")
        )

    input_df = pr.read_gtf(input_file)
    input_df = input_df.subset(lambda df: df.Chromosome.isin(chromosomes))
    input_df = input_df.drop(["Source", "Score", "Frame", "transcript_id"])
    input_df.Chromosome = input_df.Chromosome.cat.remove_unused_categories()
    input_df.Fasta = pr.get_fasta(input_df, fasta_file)

    codons = codon_list

    peptides = list()

    for seq in input_df.Fasta:
        orfs = [seq[i:] for i in range(3)]
        local_peptides = list()
        for orf in orfs:
            kmers = defaultdict(list)
            for i in range(0, len(orf) - 2, 3):
                kmers[orf[i : i + 3]].append(i)
            for c in codons:
                local_peptides += [
                    str(Seq.translate(orf[i:], to_stop=True)) for i in kmers[c]
                ]
        peptides.append(list(set(local_peptides)))

    final_df = input_df.df
    final_df["Sequence"] = peptides
    final_df = final_df.explode("Sequence")
    final_df = final_df[final_df["Sequence"].str.strip().astype(bool)]
    final_df = final_df[
        (final_df["Sequence"].str.len() <= aa_length[1])
        & (final_df["Sequence"].str.len() > aa_length[0])
    ].reset_index(drop=True)
    final_df["Protein_ID"] = "Protein_" + pd.Series(
        np.arange(1, len(final_df) + 1, 1)
    ).astype(str)
    final_df.to_csv(results_name, header=True, index=False, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fasta_file",
        "-ff",
        type=str,
        required=True,
        help="Provide the genome fasta file to extract the nucleotide sequence from. Must be uncompressed.",
    )

    parser.add_argument(
        "--codons",
        "-c",
        type=str,
        default=["ATG"],
        nargs="+",
        required=False,
        help="Provide the codons by which the micropeptides start in DNA sequence.",
    )

    parser.add_argument(
        "--input_file",
        "-i",
        type=str,
        required=False,
        help="Provide the GTF file obtained from running the filter_samples script.",
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
        nargs="+",
        required=False,
        default=[0, 100],
        help="Peptides with more length than the specified value will be filtered out.",
    )
    args = parser.parse_args()

    main(
        args.fasta_file,
        args.codons,
        args.input_file,
        args.output,
        args.length,
    )