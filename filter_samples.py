import os
import pandas as pd
import argparse
from pathlib import Path
import pyranges as pr
import re


def main(tracking_file, gtf_file, samples, output):
    results_name = (
        output
        if output and not os.path.isfile(output)
        else (re.sub(".gtf", "", Path(gtf_file).name) + "_filtered" + ".gtf")
    )
    if os.path.exists(gtf_file):
        gtf = pr.read_gtf(gtf_file)
        gtf = gtf[
            [
                "Chromosome",
                "Strand",
                "Start",
                "End",
                "num_samples",
                "transcript_id",
                "Feature",
            ]
        ]
        gtf_transcripts = gtf.subset(lambda df: df.Feature == "transcript")
        gtf_transcripts = gtf_transcripts[
            gtf_transcripts.num_samples.astype(int) >= filter
        ]
        tracking_df = pd.read_csv(tracking_file, sep="\t", header=None)
        tpm_cols = tracking_df.iloc[:, 4:].apply(
            lambda s: s.str.split("|", expand=True).iloc[:, 4]
        )
        tracking_df["tpm"] = tpm_cols.iloc[:, 0].str.cat(
            tpm_cols.iloc[:, 1:], sep=";", na_rep="None"
        )
        gtf_transcripts = gtf_transcripts.assign(
            "TPM",
            lambda df: df.transcript_id.map(tracking_df.set_index(0)["tpm"]),
        )
        gtf_transcripts.to_gtf(path=results_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tracking_file",
        "-tf",
        type=str,
        required=False,
        help="Provide the path to the tracking file obtained from gffcompare.",
    )
    parser.add_argument(
        "--gtf_file",
        "-gf",
        type=str,
        required=True,
        help="Provide the combined.gtf file obtained from gffcompare.",
    )

    parser.add_argument(
        "--samples",
        "-s",
        type=int,
        required=True,
        default=0,
        help="Provide the minimum number of samples for which you want to filter. Defaults to 0",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=False,
        help="Output path to save the filtered dataset file. Default format is a gtf file.",
    )

    args = parser.parse_args()

    main(args.tracking_file, args.gtf_file, args.samples, args.output)
