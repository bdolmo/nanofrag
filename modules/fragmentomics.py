import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import bisect
from sklearn.cluster import KMeans
import numpy as np


def create_windows(ann_dict, window_size, windows_bed):
    """
    Create windows of given size across the genome using bedtools.
    """
    windows_tmp = windows_bed.replace(".bed", ".tmp.bed")

    cmd = f'bedtools makewindows -g {ann_dict["chromosomes"]} -w {window_size} > {windows_tmp}'
    subprocess.run(cmd, shell=True, check=True)

    cmd = f'bedtools subtract -a {windows_tmp} -b {ann_dict["blacklist"]} > {windows_bed}'
    subprocess.run(cmd, shell=True, check=True)

    print(f" INFO: Windows created and saved to {windows_bed}")


def load_windows(windows_bed):
    """
    Load windows from the BED file and store them in a sorted list for efficient lookup.
    """
    windows = []
    with open(windows_bed, 'r') as bed_file:
        for line in bed_file:
            chrom, start, end = line.strip().split()[:3]
            windows.append((chrom, int(start), int(end)))
    return windows


def find_window_for_read(chrom, read_start, windows, window_starts, window_ends_by_chrom):
    """
    Find the window that a given read belongs to using binary search for efficiency.
    """
    if chrom not in window_starts:
        return None

    idx = bisect.bisect_right(window_starts[chrom], read_start) - 1

    if idx >= 0 and window_starts[chrom][idx] <= read_start <= window_ends_by_chrom[chrom][idx]:
        return windows[chrom][idx]

    return None


def parse_samtools_output(output, windows):
    """
    Parse the output from samtools to compute read counts and fragment sizes for each window.
    """
    metrics = {window: {"read_count": 0, "short_fragments": 0, "long_fragments": 0} for window in windows}

    # Precompute window start positions for efficient binary search by chromosome
    window_starts = {}
    window_ends_by_chrom = {}
    windows_by_chrom = {}

    for chrom, start, end in windows:
        if chrom not in window_starts:
            window_starts[chrom] = []
            window_ends_by_chrom[chrom] = []
            windows_by_chrom[chrom] = []
        window_starts[chrom].append(start)
        window_ends_by_chrom[chrom].append(end)
        windows_by_chrom[chrom].append((chrom, start, end))

    # Process each read from the samtools output
    for line in output:
        fields = line.split("\t")
        chrom = fields[2]
        read_start = int(fields[3])
        fragment_size = abs(int(fields[8]))  # Fragment length (TLEN field)

        # Find the window that this read belongs to
        window = find_window_for_read(chrom, read_start, windows_by_chrom, window_starts, window_ends_by_chrom)

        if window:
            # Update metrics for the window
            metrics[window]["read_count"] += 1

            # Count short vs. long fragments (DELFI-like fragment size separation)
            if 0 < fragment_size <= 150:
                metrics[window]["short_fragments"] += 1
            elif fragment_size > 150:
                metrics[window]["long_fragments"] += 1

    # Now calculate the DELFI ratio for each window
    for window_key, data in metrics.items():
        short_fragments = data["short_fragments"]
        long_fragments = data["long_fragments"]
        total_fragments = short_fragments + long_fragments

        if total_fragments > 0:
            fragment_ratio = short_fragments / total_fragments
        else:
            fragment_ratio = 0

        metrics[window_key]["fragment_size_ratio"] = fragment_ratio

    return metrics


def get_fragmentation_metrics_for_chrom(bam_file, windows_bed, chrom):
    """
    Use samtools to count reads and compute fragment sizes for a specific chromosome.
    """
    # Filter windows for the specific chromosome
    windows = load_windows(windows_bed)
    windows = [window for window in windows if window[0] == chrom]

    # Use samtools to stream the BAM file reads for the chromosome
    cmd = f"samtools view -L {windows_bed} {bam_file} {chrom}"
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True) as proc:
        window_metrics = parse_samtools_output(proc.stdout, windows)

    return window_metrics


def save_metrics_to_file(metrics, output_file):
    """
    Save the computed metrics to a file.
    
    Parameters:
    metrics (dict): Dictionary with window metrics (read counts, fragment sizes, and DELFI ratio).
    output_file (str): Path to the output file.
    """
    with open(output_file, 'w') as out:
        out.write("Chromosome\tStart\tEnd\tRead_Count\tShort_Fragments\tLong_Fragments\tFragment_ratio\n")
        for (chrom, start, end), data in metrics.items():  # Unpacking three values
            out.write(f"{chrom}\t{start}\t{end}\t{data['read_count']}\t{data['short_fragments']}\t{data['long_fragments']}\t{data['fragment_size_ratio']:.4f}\n")

    msg = f" INFO: Fragmentation metrics saved to {output_file}"
    print(msg)

import re

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
            cmd = f"samtools view {bam_file}"
            with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True) as proc:
                count = 0
                for line in proc.stdout:
                    if count >= limit:
                        break

                    if len(fields) < 10:
                        continue

                    fields = line.split("\t")
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

    # Save the plot
    plt.savefig(output_histogram)
    msg = f" INFO: Histogram of read sizes saved to {output_histogram}"
    print(msg)


def run_fragmentomic_analysis(sample_list, ann_dict, genome, output_dir, num_cpus, window_size=5000000):
    """ """
    for sample in sample_list:
        fragment_sizes_txt = os.path.join(output_dir, f"{sample.name}.fragment.sizes.txt")
        sample.add("fragment_sizes_txt", fragment_sizes_txt)

        # get data from fragmentation distribution
        if not os.path.isfile(fragment_sizes_txt):
            get_read_size_histogram(sample.bam, fragment_sizes_txt)
    
    fragment_png = os.path.join(output_dir, "framgmentation.histogram.png")
    plot_fragment_distribution(sample_list, fragment_png)


def plot_fragment_distribution(sample_list, fragment_png):
    """ """

    df_list = []
    for sample in sample_list:
        df = pd.read_csv(sample.fragment_sizes_txt,  names=[sample.name], header=None)
        df_list.append(df)

    result = pd.concat(df_list, axis=1)

    plt.figure(figsize=(10, 6))

    colors = ["red", "green", "blue"]
    idx = 0
    for col in result.columns:
        if col == "":
            continue
        print(col)

        color = colors[idx]

        # kmeans_input = np.array(result[col].values.tolist())
        # kmeans_input = np.reshape(kmeans_input, (-1, 2))
        # kmeans = KMeans(n_clusters=3).fit(kmeans_input)

        sns.histplot(result[col], bins=8000, label=col, kde=False, color=color, alpha=0.5)
        idx+=1

    # # Set titles and labels
    plt.title("cfDNA fragmentation", fontsize=16, weight='bold')
    plt.xlabel("Fragment Size (bp)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(0, 800)
    plt.tight_layout()
    plt.legend(title="Samples")


    # Save the plot
    plt.savefig(fragment_png)
    msg = f" INFO: Histogram of read sizes saved to {fragment_png}"
    print(msg)



        # plot_fragment_histogram(histogram_txt, histogram_png, analysis_type)

    # for bam in bams:

    #     windows_bed = os.path.join(output_dir, f"windows_{str(window_size)}.bed")
    #     sys.exit()
    #     create_windows(ann_dict, window_size, windows_bed)


    #     summary_metrics_tsv = os.path.join(output_dir, f"{bam_name}.summary.fragmentation.tsv")


    #     # Get list of chromosomes from the genome file
    #     chromosomes = [line.strip().split()[0] for line in open(ann_dict["chromosomes"], 'r')]

    #     # Create a pool of workers and run in parallel
    #     with multiprocessing.Pool(num_cpus) as pool:
    #         results = pool.starmap(process_chromosome, [(bam, windows_bed, chrom) for chrom in chromosomes])

    #     # Combine results from all chromosomes
    #     all_window_metrics = {}
    #     for chrom, chrom_window_metrics in results:
    #         all_window_metrics.update(chrom_window_metrics)

    #     save_metrics_to_file(all_window_metrics, summary_metrics_tsv)



