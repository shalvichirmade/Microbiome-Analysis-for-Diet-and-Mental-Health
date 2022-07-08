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

# Need to change path for hard drive

rule check_fastqc:
    input: 
        "Snakemake_Trial/{sample}_R1_001.fastq.gz"
        "Snakemake_Trial/{sample}_R2_001.fastq.gz"

    output:
        "Snakemake_Trial/fastqc"
    
    shell:
        "fastqc {input} -o {output}"

rule trim:
    input:
        "Snakemake_Trial/{sample}_R1_001.fastq.gz"
        "Snakemake_Trial/{sample}_R2_001.fastq.gz"

    output:
        "Snakemake_Trial/trim/{sample}.fastq.gz"

    shell:
        "java -jar "/Users/shalvichirmade/miniconda3/share/trimmomatic-0.39-2/trimmomatic.jar" PE -phred33 {input} -baseout {output} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"

rule join:
    input:
        "Snakemake_Trial/trim/*P.fastq.gz"

    output:
        "Snakemake_Trial/join/{sample}_%.fastq.gz"

    shell:
        "fastq-join {input} -o {output}"


