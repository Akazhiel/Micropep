import os
import pandas as pd
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyranges as pr
import subprocess
import numpy as np

records = list(SeqIO.parse("../renamed.fna", "fasta"))
chromosomes = []
for record in records:
    chromosomes.append(record.id)

gtf_df = pr.read_gtf("../GCF_000001405.25_GRCh37.p13_genomic_ucsc.sorted.gtf")
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

input_df = pr.read_gtf("./gffmix.combined_sample_filtered.gtf")
input_df = input_df.subset(lambda df: df.Chromosome.isin(chromosomes))
input_df = input_df.drop(["Source", "Score", "Frame", "transcript_id"])
input_df.Chromosome = input_df.Chromosome.cat.remove_unused_categories()

result_df = input_df.join(gtf_df)
result_df = result_df.drop(
    [
        "Start_b",
        "End_b",
        "Strand_b",
        "Feature",
        "Feature_b",
        "gene_id",
        "TPM",
        "num_samples",
        "Strand",
    ]
)
result_df = result_df.df
# result_df.to_csv("./overlaps", sep="\t", index=False)
print(result_df)
# result_df = input_df.coverage(gtf_df)
# print(result_df)
print(result_df.groupby(["Chromosome", "Start", "End", "gene_biotype"]))