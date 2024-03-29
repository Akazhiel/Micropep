![Version](https://img.shields.io/badge/Version-0.1.0-green) ![License](https://img.shields.io/badge/license-MIT-green)

# Micropep

A pipeline to produce micropeptide sequences from **BAM** files derived from an alignment step.

<p align="center">
<img src="https://user-images.githubusercontent.com/34244350/97598261-0bccda80-1a07-11eb-9767-15e287ff5c85.png" height="600">
</p>

## Notes

- This pipeline must be run using the same genome and annotation file that were used in the generation of the BAM files in the alignment step.

## Requierements

- Python3 (>=3.6)
- pyranges
- GffCompare
- StringTie
- BioPython
- Pandas
- Numpy
- NCBI-Blast+

## Installation

As for the scripts present in this module, they can be installed running the following one-liner just once:

```bash
python setup.py install
```

Pyranges, StringTie, Biopython, Pandas, Numpy and NCBI-Blast+ can all be installed one by one through the anaconda package manager. Or creating a new environment with the `.yml` file that is provided.

GffCompare latest version is not yet available at the anaconda cloud, thus why it the source code must be downloaded from their [website](http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.12.1.tar.gz) and compiled in your machine. For the installation proceed with the following steps:

```bash
tar -xf gffcompare-*.tar.gz
cd gffcompare
make release
```
## Docker

A Dockerfile is provided with this repository. The image that is built using that file will create a docker container in which you can find an already activated conda environment with all the needed dependencies to run this pipeline.

On the other hand, an already built and functional image can be downloaded from the docker hub.

```bash
docker run -v /path/to/data_folder:/Micropeptide --name Micropeptide -it akazhiel/micro:latest
```

## Usage

1. Run StringTie to assemble the transcriptome guided with the annotation file

```bash
stringtie /path/to/bam -G /path/to/annotation_file -o /path/to/output
```

2. Run the filter_tpm.py to filter out the transcripts that have a value below the chosen percentile (default: 10%). This step must be run using the GTF derived from StringTie as input. Outputs a GTF file with the annotated and putative novel assembled transcripts.

```bash
filter_tpm.py -f /path/to/stringtie_gtf --tpm <percentile> -o /path/to/output
```

3. Run GffCompare using all the GTFs derived from the previous step as input, and the annotation file to guide the transcripts overlap. The input must be a file with the paths to the GTF files, one per line.

```bash
gffcompare -r /path/to/annotation_file -i <input_gtf_list>
```

4. Filter the `gffcmp.combined.gtf` by the minimum amount of samples in which you want the transcripts to be present in.

```bash
filter_samples.py -tf /path/to/gffcmp.tracking -gf /path/to/gffcmp.combined.gtf -s min_num_of_samples -o /path/to/output
```

5. Obtain the aminoacid sequence of the transcripts obtained in the previous step. In this step the filter `prot_size` can be applied. By default, all peptides that are above 100aa will be filtered out. In this step, we used the 4 most common start codons, as described in the paper by [Slavoff](https://europepmc.org/article/med/32209305). The output of this step is a tab-delimited csv file. As mentioned above, the fasta file to use for this step must be the same one that was used for aligning the reads.

```bash
pepget.py --fasta_file /path/to/genome_fasta --annotation_gtf /path/to/annotation_file -i /path/to/gtf_file -l <max_aas_length> -o /path/to/output
```

6. Blast the obtained protein sequences against the database of your choice. This script includes two flags to filter the results obtained. One flag is for the amount of mismatches allowed (Protein sequences with more than the specified amount will be filtered out), and a boolean flag to specify if we want to exclude these results from the final file or in the other hand, keep only these ones and filter out the rest.

```bash
pblast.py -i /path/to/input.csv -db /path/to/database --exclude --mismatch <num_of_mismatches> -o /path/to/output
```

In order to run this final step, the desired database must be downloaded in the machine. For more information on the databases available, please refer to the NCBI [website](https://www.ncbi.nlm.nih.gov/books/NBK62345/#blast_ftp_site.The_blastdb_subdirectory)

```bash
update_blast --decompress <db-name>
```

## License

MIT License, see LICENSE file.

## Author

Jonatan Gonzalez Rodriguez jonatan.gonzalez.r@outlook.com

