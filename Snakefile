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

# Define variables:

# SAMPLES = "N_S002_9_S57"

# rule all:
#     input:
#         expand("/Volumes/VanRaayLab/Weston Project/Pipeline_Evaluation/{sample}_R1_001.fastq.gz", sample=SAMPLES)
    # output:
    #     expand("/Volumes/VanRaayLab/Weston Project/Pipeline_Evaluation/join/{sample}_merged_fastq.gz", sample=SAMPLES)

# rule all:
#     params:
#         sample="N_S002_9_S57"

# configfile: "config.yaml"

rule check_fastqc:
    input: 
        "/Volumes/VanRaayLab/Weston_Project/Pipeline_Evaluation/N_S002_9_S57_R1_001.fastq.gz"
        "/Volumes/VanRaayLab/Weston_Project/Pipeline_Evaluation/N_S002_9_S57_R2_001.fastq.gz"

    # output:
    #     "/Volumes/VanRaayLab/Weston Project/Pipeline_Evaluation/fastqc"
    
    shell:
        "fastqc {input} -o /Volumes/VanRaayLab/Weston_Project/Pipeline_Evaluation/fastqc"

rule trim:
    input:
        "/Volumes/VanRaayLab/Weston_Project/Pipeline_Evaluation/N_S002_9_S57_R1_001.fastq.gz"
        "/Volumes/VanRaayLab/Weston_Project/Pipeline_Evaluation/N_S002_9_S57_R2_001.fastq.gz"

    output:
        "/Volumes/VanRaayLab/Weston_Project/Pipeline_Evaluation/trim/N_S002_9_S57.fastq.gz"

    shell:
        "trimmomatic PE -phred33 {input} -baseout {output} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"

rule join:
    input:
        "/Volumes/VanRaayLab/Weston_Project/Pipeline_Evaluation/trim/*P.fastq.gz"

    output:
        "/Volumes/VanRaayLab/Weston_Project/Pipeline_Evaluation/join/N_S002_9_S57_%.fastq.gz"

    shell:
        "fastq-join {input} -o {output}"

rule merge_join:
    input:
        "/Volumes/VanRaayLab/Weston_Project/Pipeline_Evaluation/join/N_S002_9_S57*.gz"

    output:
        "/Volumes/VanRaayLab/Weston_Project/Pipeline_Evaluation/join/N_S002_9_S57_merged_fastq.gz"

    shell:
        "cat {input} > {output}"

# rule metaphlan:
#     input: 
#         "/Volumes/VanRaayLab/Weston Project/Pipeline_Evaluation/join/{sample}_merged_fastq.gz"

#     output:
#         "/Volumes/VanRaayLab/Weston Project/Pipeline_Evaluation/metaphlan/{sample}_merged_profile.txt"

#     shell:
#         "conda activate mpa"
#         "metaphlan {input} --input_type fastq --bowtie2db "/Users/shalvichirmade/Documents/MBinf/BINF 6999/Weston Project/MetaPhlan 3.0" --nproc 2 > {output}"
