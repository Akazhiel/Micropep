import os
from pathlib import Path
import argparse
import pyranges as pr
import pandas as pd
import re
import numpy as np


def main(file_name, filter, output):
    results_name = (
        output
        if output and not os.path.isfile(output)
        else (re.sub(".gtf", "", Path(file_name).name) + "_filtered" + ".gtf")
    )
    print(results_name)
    if not (0 < filter < 1):
        print("Filter value must be bigger than 0 and smaller than 1.")

    else:
        gtf_df = pr.read_gtf(file_name).df
        gtf_df["TPM"] = pd.to_numeric(gtf_df["TPM"])
        gtf_df = gtf_df[gtf_df["Feature"] == "transcript"]
        gtf_filtered = gtf_df[
            gtf_df["TPM"] > np.quantile(gtf_df["TPM"], filter)
        ]
        gtf_pr = pr.PyRanges(gtf_filtered)

        gtf_pr.to_gtf(path=results_name)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--file",
        "-f",
        type=str,
        required=True,
        help="Input file/path for the script. Must be pointing to a GTF file obtained from Stringtie.",
    )
    parser.add_argument(
        "--tpm",
        type=float,
        required=True,
        help="Provide the percentile of TPMs by which the filter will be applied. Must be between 0 and 1.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=False,
        help="Output directory where you want to save the resulting filtered GTF file. If nothing is provided, the file name will be appended with an '_filtered'.",
    )
    args = parser.parse_args()

    main(args.file, args.tpm, args.output)