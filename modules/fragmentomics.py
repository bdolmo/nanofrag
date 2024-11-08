import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import bisect
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from modules.utils import create_windows_bed_wig
import numpy as np
import re
from scipy.signal import savgol_filter
import multiprocessing


# Define your function to process a single chromosome
def process_chromosome(chrom, sample_bam, windows_bed):
    """ """
    msg = f" INFO: Processing chromosome {chrom}"
    print(msg)
    return get_fragmentation_metrics_for_chrom(sample_bam, windows_bed, chrom)

# Main function to execute the multiprocessing
def parallel_process_fragments(sample, windows_bed, chromosomes, threads):
    """ """
    # Create a pool of workers
    with multiprocessing.Pool(threads) as pool:
        # Distribute the tasks across the pool of workers
        results = pool.starmap(process_chromosome, [(chrom, sample.bam, windows_bed) for chrom in chromosomes])
    return results


def create_windows(ann_dict, window_size, windows_bed):
    """
    Create windows of given size across the genome using bedtools.
    """

    min_window = 1000000

    windows_100kb = windows_bed.replace(".bed", f"{1000000}.tmp.bed")
    windows_tmp = windows_bed.replace(".bed", ".tmp.bed")

    cmd = f'bedtools makewindows -g {ann_dict["chromosomes"]} -w {min_window} | bedtools intersect -a stdin -b {ann_dict["blacklist"]} -v > {windows_bed}'
    subprocess.run(cmd, shell=True, check=True)

    # cmd = f'bedtools makewindows -g {ann_dict["chromosomes"]} -w {window_size} > {windows_tmp}'
    # subprocess.run(cmd, shell=True, check=True)

    # cmd = f'bedtools subtract -a {windows_tmp} -b {ann_dict["blacklist"]} > {windows_bed}'
    # subprocess.run(cmd, shell=True, check=True)

    msg = f" INFO: Windows created and saved to {windows_bed}"
    print(msg)


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
    metrics = {window: {"read_count": 0, "ultra_short_fragments": 0, "short_fragments": 0, "long_fragments": 0} for window in windows}

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

        if len(fields) < 10:
            continue

        sequence = fields[9]  # Sequence
        mapqual = int(fields[4])  # Mapping quality
        cigar = fields[5]  # CIGAR string
        # Get the actual fragment size by adjusting for soft-clipping
        fragment_size = calculate_actual_fragment_size(cigar, len(sequence))

        # fragment_size = abs(int(fields[8]))  # Fragment length (TLEN field)

        # Find the window that this read belongs to
        window = find_window_for_read(chrom, read_start, windows_by_chrom, window_starts, window_ends_by_chrom)

        if window:
            # Update metrics for the window
            metrics[window]["read_count"] += 1

            # Count short vs. long fragments (DELFI-like fragment size separation)

            if 0 < fragment_size <= 80:
                metrics[window]["ultra_short_fragments"] += 1
            if 80 < fragment_size <= 165:
                metrics[window]["short_fragments"] += 1
            elif fragment_size > 165:
                metrics[window]["long_fragments"] += 1

    # Now calculate the DELFI ratio for each window
    for window_key, data in metrics.items():
        ultra_short_fragments = data["ultra_short_fragments"]
        short_fragments = data["short_fragments"]
        long_fragments = data["long_fragments"]
        total_fragments = ultra_short_fragments+short_fragments + long_fragments

        if total_fragments > 0:
            fragment_ratio = (ultra_short_fragments+short_fragments) / total_fragments
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


    # Check if the output file already exists; if not, run the C++ binary
    # if not os.path.isfile(fragment_bed):
    # Run the C++ binary to generate fragment counts
    # The binary expects: <BAM file> <BED file> <output file>
    bam_file = sample.bam  # Assuming you have a method to get the BAM file for each sample
    command = [cxx_binary_path, bam_file, windows_bed, fragment_bed, fragment_wig]
    print(' '.join(command))
    # Execute the C++ binary
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


    # for sample in sample_list:

    #     fragment_bed = os.path.join(fragment_folder, f"{sample.name}.fragmentation.data.bed")
    #     sample.add("fragment_data", fragment_bed)

    #     calculate_fragment_counts(sample, fragment_folder, windows_bed, bin_dict["cfdna_counter"])

    #     # if not os.path.isfile(fragment_bed):
    #     #     o = open(fragment_bed, "a")
    #     #     o.write("chr\tpos\tend\tread_count\tultra_short_fragments\tshort_fragments\tlong_fragments\tfragment_size_ratio\n")
    #     #     results = parallel_process_fragments(sample, windows_bed, chromosomes, num_cpus)

    #     #     for chromdata in results:
    #     #         for region in chromdata:
    #     #             read_count = chromdata[region]["read_count"]
    #     #             ultra_short_fragments = chromdata[region]["ultra_short_fragments"]
    #     #             short_fragments = chromdata[region]["short_fragments"]
    #     #             long_fragments = chromdata[region]["long_fragments"]
    #     #             fragment_size_ratio = chromdata[region]["fragment_size_ratio"]
    #     #             o.write(f"{region[0]}\t{region[1]}\t{region[2]}\t{read_count}\t{ultra_short_fragments}\t{short_fragments}\t{long_fragments}\t{fragment_size_ratio}\n")
    #     #     o.close()

    #     #     with open(fragment_bed) as f:
    #     #         for line in f:
    #     #             line = line.rstrip("\n")
    #     #             if line.startswith("chr\tpos"):
    #     #                 continue
                        
    #     #     f.close()


    #     fragment_size_ratio_png = os.path.join(fragment_folder, f"{sample.name}.fsr.png")
    #     plot_fragmentation_ratio(sample.name, fragment_bed, fragment_size_ratio_png)


    for sample in sample_list:
        fragment_sizes_txt = os.path.join(fragment_folder, f"{sample.name}.fragment.sizes.txt")
        sample.add("fragment_sizes_txt", fragment_sizes_txt)

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
        df = pd.read_csv(sample.fragment_sizes_txt,  names=[sample.name], header=None)
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
        sns.histplot(result[col], bins=8000, label=col, kde=False, color=color)
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
        # print(mode.iloc[0], max_val)
        plt.axvline(x=mode.iloc[0], ymin=0, ymax=max_val, color="white", linestyle="--")
        # plt.text(mode.iloc[0]+2, max_val+2, mode.iloc[0], rotation=60, va='center')

    plt.tight_layout()
    plt.legend(title="Samples")

    # Save the plot
    plt.savefig(fragment_png)
    msg = f" INFO: Histogram of read sizes saved to {fragment_png}"
    print(msg)





