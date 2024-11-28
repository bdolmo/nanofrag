import os
import sys
import subprocess
import re
from collections import defaultdict
from bx.intervals.intersection import Intersecter, Interval
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from modules.utils import get_chromosome_sizes, create_windows_bed_wig
import gzip
import multiprocessing
import math
import subprocess
import rdata



def run_ichorcna_docker(input_bam, output_dir, docker_output, wig_file_path, sample_id="tumor_sample"):
    # Define paths within the container
    # wig_file_path = f"{output_dir}/{sample_id}.wig"
    
    # Step 1: Run readCounter to generate .wig file
    # readcounter_command = [
    #     "docker", "run", "-it", "-v", f"{input_bam}:/data/tumor.bam",
    #     "-v", f"{output_dir}:/output", "gavinhalab/ichorcna:1.0.0",
    #     "readCounter", "--window", "1000000", "--quality", "20",
    #     "--chromosome", "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y",
    #     "/data/tumor.bam", ">", "/output/tumor.wig"
    # ]
    print(input_bam)
    bam_name = os.path.basename(input_bam)
    fragment_folder = os.path.dirname(wig_file_path)
    # readcounter_command = [
    #     "docker", "run", "-it", "-v", f"{os.path.dirname(input_bam)}:/bam_dir",
    #     "-v", f"{output_dir}:/output", "gavinhalab/ichorcna:1.0.0",
    #     f' /bin/bash -c "readCounter --window 1000000 --quality 20 --chromosome chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY  /bam_dir/{bam_name}> /output/{sample_id}.wig"'
    # ]

    # if not os.path.isfile(wig_file_path):
    #     try:
    #         print("Running readCounter...")
    #         subprocess.run(" ".join(readcounter_command), shell=True, check=True)
    #         print(f"ReadCounter completed. Output saved to {wig_file_path}")
    #     except subprocess.CalledProcessError as e:
    #         print(f"Error in readCounter command: {e}")
    #         return
    wig_name = os.path.basename(wig_file_path)


    seg_file = os.path.join(output_dir, f"{sample_id}.seg.txt")
    if not os.path.isfile(seg_file):
        # Step 2: Run ichorCNA with the generated .wig file
        ichorcna_command = [
            "docker", "run", "-it", "-v", f"{docker_output}/CNA:/output", "-v", f"{docker_output}/FRAGMENTATION:/fragment_folder",
            "seqeralabs/ichorcna", "runIchorCNA.R",
            "--id", sample_id, "--WIG", f"/fragment_folder/{wig_name}", "--ploidy", "\"c(2,3)\"",
            "--normal", "\"c(0.5,0.6,0.7,0.8,0.9)\"", "--maxCN", "5",
            "--gcWig", "/opt/conda/share/r-ichorcna-0.1.0.20180710-0/extdata/gc_hg38_1000kb.wig",
            "--mapWig", "/opt/conda/share/r-ichorcna-0.1.0.20180710-0/extdata/map_hg38_1000kb.wig",
            "--centromere", "/opt/conda/share/r-ichorcna-0.1.0.20180710-0/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt",
            "--includeHOMD", "False", "--estimateNormal", "True",
            "--estimatePloidy", "True", "--estimateScPrevalence", "True",
            "--scStates", "\"c(1,3)\"", "--txnE", "0.9999", "--txnStrength", "10000",
            "--outDir", "/output"
        ]
        print(" ".join(ichorcna_command))

        try:
            print("Running ichorCNA analysis...")
            # subprocess.run(ichorcna_command, check=True)
            subprocess.run(" ".join(ichorcna_command), shell=True, check=True)
            print(f"ichorCNA analysis completed. Results saved in {output_dir}")
        except subprocess.CalledProcessError as e:
            print(f"Error in ichorCNA command: {e}")

    cna_plot = seg_file.replace(".seg.txt", ".cna.png")
    rdata_file = seg_file.replace(".seg.txt", ".RData")
    plot_cna_genomewide(sample_id, rdata_file, cna_plot)



def normalized_bed_to_dict(input_bed):
    """ """
    data_list = []
    header_dict = {}
    with open(input_bed) as f:
        for line in f:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            data_dict = {}

            if line.startswith("chr\t"):
                for idx,item in enumerate(tmp):
                    header_dict[idx] = item
                    data_dict[item] = ""
                continue
            for idx,item in enumerate(tmp):
                field_name = header_dict[idx]
                data_dict[field_name] = item
            data_list.append(data_dict)
    f.close()
    return data_list


def calculate_log2_ratios(sample_list, output_dir):
    """ """

    tumor_samples = []
    normal_samples = []
    for sample in sample_list:
        if sample.origin == "tumor":
            tumor_ratios_bed = os.path.join(output_dir, "CNA", f"{sample.name}.log2.ratios.bed")
            sample.add("log2_ratios", tumor_ratios_bed )
            tumor_samples.append(sample)
        else:
            normal_samples.append(sample)
    
    for tumor_sample in tumor_samples:

        # print(tumor_sample.name)

        normal_sample = normal_samples[0]

        tumor_data = normalized_bed_to_dict(tumor_sample.normalized_bed)
        normal_data = normalized_bed_to_dict(normal_sample.normalized_bed)

        header_tmp = []
        for field in tumor_data[0]:
            header_tmp.append(field)
        header_tmp.append("log2_ratio")
        header_str = "\t".join(header_tmp)
        # print(header_str)

        o = open(tumor_sample.log2_ratios, "w")
        o.write(header_str+"\n")
        for idx, row in enumerate(tumor_data):
            tumor_rd = float(tumor_data[idx]["short_fragments"])
            normal_rd = float(normal_data[idx]["short_fragments"])

            if tumor_rd == 0 or normal_rd == 0:
                log2_ratio = -3
            else:
                log2_ratio = math.log2(tumor_rd/normal_rd)
            tmp_data = []
            for field in tumor_data[idx]:
                value = tumor_data[idx][field]
                tmp_data.append(value)
            tmp_data.append(str(log2_ratio))
            tmp_str = "\t".join(tmp_data)
            o.write(tmp_str+"\n")
        o.close()

    return sample_list
            
 
def normalize_raw_depth(sample_name, input_bed, output_bed):
    """ """

    counts_by_length = []
    with open(input_bed) as f:
        for line in f:
            if line.startswith("chr\tpos"):
                # header_str = line
                continue
            tmp = line.split("\t")
            # raw_counts = int(tmp[3])
            raw_counts = int(tmp[5])

            bin_size = int(tmp[2])-int(tmp[1])
            counts_by_length.append(1000000*(raw_counts/bin_size))
    f.close()

    median_short_counts = np.median(counts_by_length)

    o = open(output_bed, "w")
    with open(input_bed) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("chr\tpos"):
                o.write(line+"\t"+"cn_status"+"\n")
                continue
            tmp = line.split("\t")
            # tmp[3] = str((int(tmp[3])/median_counts))
            tmp[5] = str((int(tmp[5])/median_short_counts))

            normalized_depth = float(tmp[5])
            cn_status = "Diploid"

            if normalized_depth < 0.9:
                cn_status = "Loss"

            if normalized_depth > 1.1:
                cn_status = "Gain"
            tmp.append(cn_status)
            line = "\t".join(tmp)
            o.write(line+"\n")
    f.close()
    o.close()


def plot_cn_profile_vs_baseline(sample_name, input_bed, output_png):

    df = pd.read_csv(input_bed, sep="\t", header=0, names=["chr", "pos", "end", "read_count", 
        "ultra_short_fragments", "short_fragments", "long_fragments","fragment_size_ratio", "cn_status", "log2_ratio"])

    chromosomes = df['chr'].tolist()
    chr_colors = {}
    chr_limits = {}
    
    idx = 0
    chr_count = 0
    unique_chromosomes = []
    ticks = []
    for chrom in chromosomes:
        chr_count += 1

        color = "#686868"
        idx = 0
        if not chrom in chr_colors:
            chr_colors[chrom] =color
            unique_chromosomes.append(chrom)
        if not chrom in chr_limits:
            chr_limits[chrom] = chr_count
            ticks.append(chr_count)

    cn_status_colors = {
        "Loss": "red",
        "Diploid": "#686868",
        "Gain": "green"
    }
    print(sample_name, "variance inter:", df["log2_ratio"].var())

    plt.figure(figsize=(20, 5))
    ax = sns.scatterplot( x=df.index, y=df["log2_ratio"], s=4, hue=df["cn_status"], palette=cn_status_colors)

    ax.set_xticks(ticks, unique_chromosomes, rotation=45)


    # Set titles and labels
    plt.title(f"CNA profile for sample {sample_name}", fontsize=16, weight='bold')
    plt.ylabel("Log2 Ratio", fontsize=14)
    plt.ylim(-1,1.5)
    # plt.figure(figsize=(20, 4))

    for chrom in chr_limits:
        plt.axvline(x=chr_limits[chrom], ymin=0, ymax=3, color="grey", linestyle="--")

    plt.legend([],[], frameon=False)

    # Save the plot
    plt.savefig(output_png)
    plt.close()


def plot_cna_genomewide(sample_name, rdata_file, output_png):
    """
    Plots the copy number profile for a given sample
    """

    # Convert RData to CSV if not already converted
    seg_csv = rdata_file.replace(".RData", ".csv")

    # if not os.path.isfile(seg_csv):
    parsed = rdata.parser.parse_file(rdata_file)
    converted = rdata.conversion.convert(parsed)

    for item in converted["results"]:
        data = (item["cna"][sample_name])
        break

    data.to_csv(seg_csv, sep='\t')
    max_idx = (converted['loglik']['loglik'].idxmax())-1
    selected_row = converted["loglik"].iloc[max_idx]

    # init                     n0.8-p2
    # n_est                       0.63
    # phi_est                    1.728
    # BIC                         <NA>
    # Frac_genome_subclonal       0.17
    # Frac_CNA_subclonal          0.26
    # loglik                    1363.0

    # ploidy = "."
    # tmp_ploidy = selected_row["init"].split("-")
    # if len(tmp_ploidy) > 1:
    #     ploidy = tmp_ploidy[1].replace("p", "")
    ploidy = selected_row["phi_est"]


    n_est = selected_row["n_est"]
    t_purity = round(1-float(n_est), 3)
    genome_subclonal_est = selected_row["Frac_genome_subclonal"]
    cna_subclonal_est = selected_row["Frac_CNA_subclonal"]
    loglik = selected_row["loglik"]


    # Load data
    df = pd.read_csv(seg_csv, sep="\t", header=0)

    # Ensure data types and add "chr" prefix
    df['chr'] = df['chr'].astype(str).apply(lambda x: f"chr{x}")

    # Create a cumulative genomic position for plotting
    df['cumulative_pos'] = df['start'] + df.groupby('chr')['start'].transform('min').cumsum()

    df["event"] = df["event"].map({'GAIN': 'Gain', 'AMP': 'Amplification', 
        'HLAMP': 'High-level Amplification', 'NEUT': 'Neutral', 'HETD': 'Loss'})

    # Define colors for the `event` categories
    event_colors = {
        "Gain": "#33b366",
        "Amplification": "green",
        "High-level Amplification": "#003200",
        "Neutral": "blue",
        "Loss": "red",
        # Add more event types if necessary
    }

    # Prepare chromosome ticks and colors
    chromosomes = df['chr'].tolist()
    chr_colors = {}
    chr_limits = {}
    ticks = []
    unique_chromosomes = []
    chr_count = 0

    for chrom in chromosomes:
        chr_count += 1
        color = "#686868"
        if chrom not in chr_colors:
            chr_colors[chrom] = color
            unique_chromosomes.append(chrom)
        if chrom not in chr_limits:
            chr_limits[chrom] = chr_count
            ticks.append(chr_count)

    # Set up the plot
    plt.figure(figsize=(18, 5))


    # Plot individual log2 ratios with event-based coloring
    sns.scatterplot(
        x='Unnamed: 0',
        y='logR',
        hue='event',
        palette=event_colors,
        data=df,
        s=9,
        alpha=0.9
    )
    sns.despine(left=False, bottom=True)
    # plt.legend(labels=['Gain', 'Amplification', 'High-level Amplification', 'Neutral', 'Loss'])


    plt.axhline(color="lightgrey")
    # Add vertical lines for chromosome boundaries
    for chrom, limit in chr_limits.items():
        plt.axvline(x=limit, color="grey", linestyle="--", linewidth=0.5)

    # Add chromosome labels
    plt.xticks(ticks, labels=unique_chromosomes, rotation=45, fontsize=14)
    plt.yticks([-2, -1, 0, 1, 2], fontsize=14) 

    # Set titles and labels
    plt.title(f"{sample_name}\nTumor purity: {t_purity} Ploidy: {ploidy}", fontsize=18)
    plt.ylabel("Log2 Ratio", fontsize=14)
    plt.xlabel("", fontsize=1)

    # plt.xlabel("Chromosomes", fontsize=14)
    plt.ylim(-2, 2)  # Adjust y-axis limits if needed

    # Add legend
    plt.legend(loc="upper right")

    # Save the plot
    plt.tight_layout()
    plt.savefig(output_png,  dpi=300)
    plt.close()



def plot_cn_profile_intrasample(sample_name, input_bed, output_png):
    """ """

    df = pd.read_csv(input_bed, sep="\t", header=0, names=["chr", "pos", "end", "read_count", 
        "ultra_short_fragments", "short_fragments", "long_fragments","fragment_size_ratio", "cn_status"])

    chromosomes = df['chr'].tolist()
    chr_colors = {}
    chr_limits = {}
    
    idx = 0
    chr_count = 0
    unique_chromosomes = []
    ticks = []
    for chrom in chromosomes:
        chr_count += 1
        # if idx == 0:
        #     color = "#e9e9e9"
        #     idx =+ 1
        # else:
        color = "#686868"
        idx = 0
        if not chrom in chr_colors:
            chr_colors[chrom] =color
            unique_chromosomes.append(chrom)
        if not chrom in chr_limits:
            chr_limits[chrom] = chr_count
            ticks.append(chr_count)

    cn_status_colors = {
        "Loss": "red",
        "Diploid": "#686868",
        "Gain": "green"
    }

    df["short_fragments"] = np.log2(df["short_fragments"])
    sdata = df["short_fragments"].dropna()
    # print(sample_name, "variance intra:", sdata.var())

    plt.figure(figsize=(18, 5))
    ax = sns.scatterplot( x=df.index, y=df["short_fragments"], s=4, hue=df["cn_status"], palette=cn_status_colors)

    # ax.set_xticklabels(unique_chromosomes, rotation=45)
    ax.set_xticks(ticks, unique_chromosomes, rotation=45)

    # Set titles and labels
    plt.title(f"CNA profile for sample {sample_name}", fontsize=16, weight='bold')
    plt.ylabel("Log2 Ratio", fontsize=14)
    plt.ylim(-1, 1.5)

    for chrom in chr_limits:
        plt.axvline(x=chr_limits[chrom], ymin=0, ymax=3, color="grey", linestyle="--")

    plt.legend([],[], frameon=False)

    # Save the plot
    plt.savefig(output_png)
    plt.close()


def run_cn_workflow(sample_list, docker_output, ann_dict,  output_dir):
    """ """

    cna_folder = os.path.join(output_dir, "CNA")
    if not os.path.isdir(cna_folder):
        os.mkdir(cna_folder)

    windows_bed = os.path.join(output_dir, "windows.1000kb.bed")
    windows_wig = os.path.join(output_dir, "windows.1000kb.wig")

    # create_windows_bed_wig(chrom_sizes_file, window_size, output_wig_file, output_bed_file)
    # create_windows_bed_wig(ann_dict["chromosomes"], 1000000, windows_wig, windows_bed)

    # for sample in sample_list:
    #     fragment_bed = os.path.join(output_dir, "FRAGMENTATION",
    #         f"{sample.name}.fragmentation.data.bed")
    #     normalized_bed = os.path.join(cna_folder,  f"{sample.name}.normalized.bed")
    #     sample.add("normalized_bed", normalized_bed)
    #     # normalize raw data
    #     normalize_raw_depth(sample.name, fragment_bed, normalized_bed)

    # sample_list = calculate_log2_ratios(sample_list, output_dir)

    # for sample in sample_list:
    #     if sample.origin == "tumor":
    #         # copy number plot
    #         cn_png =  os.path.join(cna_folder, f"{sample.name}.cn.png")
    #         # plot_cn_profile_vs_baseline(sample.name, sample.log2_ratios, cn_png)

    #     cn_png =  os.path.join(cna_folder, f"{sample.name}.intrasample.cn.png")
    #     # plot_cn_profile_intrasample(sample.name, sample.normalized_bed, cn_png)


    #plot_cn_profile(sample.name, normalized_bed, cn_png)
    for sample in sample_list:
        # if sample.origin == "tumor":
        run_ichorcna_docker(sample.bam, cna_folder, docker_output, sample.fragment_wig, sample.name)


    return sample_list
