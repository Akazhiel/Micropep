import os
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main(input_csv, database, exclude, missmatch, output):

    file_name = (
        output if output and os.path.isfile(output) else "pblast_results.csv"
    )

    prot = pd.read_csv(input_csv, sep="\t")

    records = list()
    for index, seq in enumerate(prot["Sequence"]):
        records.append(
            SeqRecord(
                Seq(seq),
                id=("Protein_" + str(index + 1)),
                description="",
                name="",
            )
        )

    with open("Protein.fa", "w") as handle:
        SeqIO.write(records, handle, "fasta")

    subprocess.Popen(
        'blastp -query Protein.fa -db {0} -max_target_seqs 1 -taxids 9606 -outfmt "6" -out blast_res.tab -num_threads 4 -qcov_hsp_perc 100'.format(
            database
        ),
        stderr=subprocess.PIPE,
        shell=True,
    ).communicate()

    subprocess.run("rm Protein.fa", shell=True)

    blast_results = pd.read_csv("blast_res.tab", sep="\t", header=None)
    blast_results = blast_results[blast_results[4] <= missmatch][0]

    if exclude:
        prot = prot[~prot["Protein_ID"].isin(blast_results)]
    else:
        prot = prot[prot["Protein_ID"].isin(blast_results)]

    prot.to_csv(file_name, header=True, index=False, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--input",
        "-i",
        metavar="\b",
        type=str,
        required=True,
        help="Provide the path to the input file in which all aminoacids sequences are present under a column named 'Sequence'",
    )

    parser.add_argument(
        "--database",
        "-db",
        metavar="\b",
        type=str,
        required=True,
        help="Provide the path to the database the blast process will be run against.",
    )

    parser.add_argument(
        "--exclude",
        action="store_true",
        help="Indicate if the matching results from the blast will be excluded from the output file.",
    )

    parser.add_argument(
        "--mismatch",
        "-mm",
        type=int,
        default=1,
        required=False,
        metavar="\b",
        help="Blast results with more than %(default)s mismatches won't be excluded.",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=False,
        metavar="\b",
        help="Provide the output to the produced file. If not provided, result csv will be saved directly to the working directory.",
    )

    args = parser.parse_args()

    main(args.input, args.database, args.exclude, args.missmatch, args.output)
