import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import bisect
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from modules.utils import create_windows_bed_wig, tissues_dict
import numpy as np
import re
from scipy.signal import savgol_filter
import multiprocessing


def calculate_actual_fragment_size(cigar, sequence_length):
    """
    Calculate the actual fragment size by subtracting soft-clipped portions from the ends.
    
    Parameters:
    cigar (str): The CIGAR string from the SAM/BAM file.
    sequence_length (int): The original sequence length.
    
    Returns:
    int: The adjusted fragment size excluding soft-clipped ends.
    """
    # Regex to find soft-clipping at the start (beginning of the string) and end (end of the string)
    softclip_start = re.match(r"^(\d+)S", cigar)  # Matches soft-clipping at the start (e.g., 5S...)
    softclip_end = re.search(r"(\d+)S$", cigar)   # Matches soft-clipping at the end (e.g., ...5S)

    # Subtract soft-clipping at the start
    if softclip_start:
        sequence_length -= int(softclip_start.group(1))

    # Subtract soft-clipping at the end
    if softclip_end:
        sequence_length -= int(softclip_end.group(1))

    return sequence_length


def get_read_size_histogram(bam_file, output_txt, limit=5000000):
    """
        Plot a histogram of read sizes, skipping reads with mapping quality < 20, using the first 50 million reads.
    """

    if not os.path.isfile(output_txt):
        with open(output_txt, 'w') as out_file:
            # Use samtools to stream the BAM file and limit to first 50 million reads
            cmd = f"samtools view {bam_file} chr1"
            with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True) as proc:
                count = 0
                for line in proc.stdout:
                    if count >= limit:
                        break

                    fields = line.split("\t")
                    if len(fields) < 10:
                        continue

                    sequence = fields[9]  # Sequence
                    mapqual = int(fields[4])  # Mapping quality
                    cigar = fields[5]  # CIGAR string
                    # Get the actual fragment size by adjusting for soft-clipping
                    fragment_size = calculate_actual_fragment_size(cigar, len(sequence))
                    # Skip low-quality reads (mapqual < 20)
                    if mapqual >= 20:
                        out_file.write(f"{fragment_size}\n")
                        count+=1
                    
        msg = f" INFO: Fragment sizes saved to {output_txt}"
        print(msg)
    
    return output_txt


def plot_fragment_histogram(input_file, output_png, analysis_type):
    """ 
    """
    fragment_sizes = pd.read_csv(input_file, header=None, names=['Fragment_Size'])

    plt.figure(figsize=(10, 6))
    sns.histplot(fragment_sizes['Fragment_Size'], bins=7000, kde=False, color="blue")

    # Set titles and labels
    plt.title("Distribution of fragment size)", fontsize=16, weight='bold')
    plt.xlabel("Fragment Size (bp)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(0, 800)
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_histogram)

    msg = f" INFO: Histogram of read sizes saved to {output_histogram}"
    print(msg)


def calculate_fragment_counts(sample, fragment_folder, windows_bed, cxx_binary_path):

    fragment_bed = os.path.join(fragment_folder, f"{sample.name}.fragmentation.data.bed")
    fragment_wig = os.path.join(fragment_folder, f"{sample.name}.fragmentation.wig")
    sample.add("fragment_data", fragment_bed)
    sample.add("fragment_wig", fragment_wig)

    bam_file = sample.bam  # Assuming you have a method to get the BAM file for each sample
    command = [cxx_binary_path, bam_file, windows_bed, fragment_bed, fragment_wig]
    print(' '.join(command))
    # Execute the C++ binary
    if not os.path.isfile(fragment_bed) and not os.path.isfile(fragment_wig):
        subprocess.run(command, check=True)


def process_sample(sample, fragment_folder, windows_bed, cfdna_counter_path):
    # Define the output file for the current sample
    fragment_bed = os.path.join(fragment_folder, f"{sample.name}.fragmentation.data.bed")
    fragment_wig = os.path.join(fragment_folder, f"{sample.name}.fragmentation.wig")
    sample.add("fragment_data", fragment_bed)
    sample.add("fragment_wig", fragment_wig)

    # Call the function to calculate fragment counts (this runs the C++ binary)
    calculate_fragment_counts(sample, fragment_folder, windows_bed, cfdna_counter_path)

    # Generate the output PNG for the fragment size ratio plot
    fragment_size_ratio_png = os.path.join(fragment_folder, f"{sample.name}.fsr.png")
    plot_fragmentation_ratio(sample.name, fragment_bed, fragment_size_ratio_png)



def run_fragmentomic_analysis(sample_list, ann_dict, bin_dict, genome, output_dir, num_cpus, window_size=5000000):
    """ """

    fragment_folder = os.path.join(output_dir, "FRAGMENTATION")
    if not os.path.isdir(fragment_folder):
        os.mkdir(fragment_folder)

    windows_bed = os.path.join(output_dir, "windows.1000kb.bed")
    windows_wig = os.path.join(output_dir, "windows.1000kb.wig")

    create_windows_bed_wig(ann_dict["chromosomes"], 1000000, windows_wig, windows_bed)

    chromosomes = [(f"chr{i}") for i in range(1, 23)]
    chromosomes.append("chrX")
    chromosomes.append("chrY")

    # Create a list of arguments for each sample
    args = [(sample, fragment_folder, windows_bed, bin_dict["cfdna_counter"]) for sample in sample_list]

    # Use a multiprocessing Pool to process each sample in parallel
    with multiprocessing.Pool(processes=num_cpus) as pool:
        pool.starmap(process_sample, args)

    for sample in sample_list:
        fragment_sizes_txt = os.path.join(fragment_folder, f"{sample.name}.fragment.sizes.txt")
        sample.add("fragment_sizes_txt", fragment_sizes_txt)

        fragment_bed = os.path.join(fragment_folder, f"{sample.name}.fragmentation.data.bed")
        fragment_wig = os.path.join(fragment_folder, f"{sample.name}.fragmentation.wig")
        sample.add("fragment_data", fragment_bed)
        sample.add("fragment_wig", fragment_wig)


        # get data from fragmentation distribution
        if not os.path.isfile(fragment_sizes_txt):
            get_read_size_histogram(sample.bam, fragment_sizes_txt)

    fragment_png = os.path.join(fragment_folder, "fragmentation.histogram.png")
    if not os.path.isfile(fragment_png):
        plot_fragment_distribution(sample_list, fragment_png)

    return sample_list
        

def plot_fragmentation_ratio(sample_name, input_bed, output_png):
    """ """

    df = pd.read_csv(input_bed, sep="\t", header=0, names=["chr", "pos", "end", "read_count", 
        "ultra_short_fragments", "short_fragments", "long_fragments", "fragment_size_ratio"])

    df["fsr_zscore"] = ((df["fragment_size_ratio"]-df["fragment_size_ratio"].mean())/df["fragment_size_ratio"].std())
    df["fsr_zscore"] = savgol_filter(df["fsr_zscore"], 15, 3)

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

    plt.figure(figsize=(20, 5))
    ax = sns.lineplot( x=df.index, y=df["fsr_zscore"])

    # ax.set_xticklabels(unique_chromosomes, rotation=45)
    ax.set_yticks([-3, 0, 3], ["-3", "0", "3"], fontsize=12)

    ax.set_xticks(ticks, unique_chromosomes, rotation=45, fontsize=12)


    # Set titles and labels
    plt.title(f"Fragmentation ratio for sample {sample_name}", fontsize=16, weight='bold')
    plt.ylabel("Frag Size Ratio (z-score)", fontsize=14)
    plt.ylim(-3.2, 3.2)

    for chrom in chr_limits:
        plt.axvline(x=chr_limits[chrom], ymin=0, ymax=3, color="grey", linestyle="--")

    plt.legend([],[], frameon=False)

    # Save the plot
    plt.savefig(output_png)
    plt.close()


def plot_fragment_distribution(sample_list, fragment_png):
    """ """

    df_list = []
    for sample in sample_list:
        df = pd.read_csv(sample.fragment_sizes_txt,  names=[sample.name], header=None,  nrows=5000000)
        df_list.append(df)

    result = pd.concat(df_list, axis=1)

    plt.figure(figsize=(10, 6))

    red_palette = [
        "#D10000",
        "#FF0000",  # Pure red
        "#FF3333",
        "#FF6666",
        "#FF9999",
        "#FFCCCC",
        "#FFE5E5"   # Lightest red
    ]

    blue_palette = [
        "#0000CC",
        "#0000FF",  # Pure blue
        "#3333FF",
        "#6666FF",
        "#9999FF",
        "#CCCCFF",
        "#E5E5FF"   # Lightest blue
    ]

    colors = []
    idx = 0
    jdx = 0
    for sample in sample_list:
        if sample.origin == "tumor":
            colors.append(red_palette[idx])
            idx+=1
        else:
            colors.append(blue_palette[jdx])
            jdx+=1

    # colors = ["red", "darkred", "blue", "darkblue"]
    idx = 0
    mode_vals = []
    max_val = 0
    max_vals = []

    means_list = []
    covariances_list = []
    weights_list = []

    interval_stats = []

    # Define the intervals
    intervals = [(50, 220), (220, 400), (400, 800)]

    for col in result.columns:
        if col == "":
            continue
        mode_val = result[col].mode()
        mode_vals.append(mode_val)
        max_val = max(result[col])
        max_vals.append(max_val)

        color = colors[idx]
        col_interval_stats = {}
        for interval in intervals:
            lower, upper = interval
            interval_data = result[col][(result[col] >= lower) & (result[col] < upper)].dropna()
            interval_median = interval_data.median()
            interval_mean = interval_data.mean()
            interval_std = interval_data.std()
            interval_mode = interval_data.mode().iloc[0]
            col_interval_stats[f"{lower}-{upper}"] = {
                'sample': col,
                'mean': interval_mean, 
                'std': interval_std, 
                'median': interval_median, 
                'mode': interval_mode,
            }
        interval_stats.append(col_interval_stats)

        # sns.displot(data=result, x=col, kind="kde",  color=color)
        # sns.kdeplot(
        #     data=result, x=col,
        #     fill=True, common_norm=False, color=color,
        #     alpha=.5, linewidth=0,
        # )

        sns.histplot(result[col], bins=8000, label=col, fill=False, stat='density', color=color, alpha=.5)
        idx+=1

    summary_name = os.path.join(os.path.dirname(fragment_png), "fragmentation.summary.csv")
    o = open(summary_name, "w")
    o.write("Name\tInterval\tMode\tMedian\tMean\tStd\n")

    print(" INFO: Fragment interval stats:")
    for stats in interval_stats:
        for interval in intervals:
            interval_str = str(interval).replace(", ", "-").replace("(", "").replace(")", "")
            print(f"\t({interval_str} bp) ",
                f"mode: {stats[interval_str]['mode']}",
                f"median: {stats[interval_str]['median']}", 
                f"mean: {stats[interval_str]['mean']}",
                f"std: {stats[interval_str]['std']}")

            out_str = f"{stats[interval_str]['sample']}\t{interval_str}\t{stats[interval_str]['mode']}\t{stats[interval_str]['median']}\t{stats[interval_str]['mean']}\t{stats[interval_str]['std']}"
            o.write(out_str+"\n")
    o.close()
            
    plt.title("cfDNA fragment distribution", fontsize=16, weight='bold')
    plt.xlabel("Fragment Size (bp)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(0, 800)

    for idx,mode in enumerate(mode_vals):
        max_val = max_vals[idx]
        plt.axvline(x=mode.iloc[0], ymin=0, ymax=max_val, color="white", linestyle="--")

    plt.tight_layout()
    plt.legend(title="Samples")

    # Save the plot
    plt.savefig(fragment_png)
    msg = f" INFO: Histogram of read sizes saved to {fragment_png}"
    print(msg)





