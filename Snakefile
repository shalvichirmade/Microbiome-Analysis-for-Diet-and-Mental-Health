## Pipeline to extract the relative abundance taxonomy from whole metagenome sequences.

## Author: Shalvi Chirmade
## Date created: July 5, 2022

# ---------------------------------------------------------------------

# User should have these software downloaded using conda:
    # fastqc
    # trimmomatic
    # fastq-join
    # metaphlan + required bowtie databases

# ---------------------------------------------------------------------

rule check_fastqc:
    input:
        "Snakemake_Trial/{sample}_R1_001.fastq.gz"
        "Snakemake_Trial/{sample}_R2_001.fastq.gz"

    output:
        "Snakemake_Trial/fastqc"
    
    shell:
        "fastqc {input} -o {output}"

