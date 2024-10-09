import os
import sys
import subprocess
import re
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import bisect
from scipy.signal import savgol_filter


def get_corrected_positions(cigar, pos, sequence_length):
    """
    Calculate the actual fragment size by subtracting soft-clipped portions from the ends.
    """
    # Regex to find soft-clipping at the start (beginning of the string) and end (end of the string)
    softclip_start = re.match(r"^(\d+)S", cigar)  # Matches soft-clipping at the start (e.g., 5S...)
    softclip_end = re.search(r"(\d+)S$", cigar)   #} Matches soft-clipping at the end (e.g., ...5S)

    start = pos
    end = pos+sequence_length

    # Subtract soft-clipping at the start
    if softclip_start:
        s_start = softclip_start.group(1)
        sequence_length -= int(softclip_start.group(1))
        start = start+int(s_start)
    # Subtract soft-clipping at the end
    if softclip_end:
        s_end = softclip_end.group(1)
        sequence_length -= int(softclip_end.group(1))
        end = end-int(s_end)

    out_dict = {
        "start": start,
        "end": end
    }
    return out_dict





def create_wps_windows(chrom_sizes, selected_chromosome):
    wps = defaultdict(dict)
    
    with open(chrom_sizes) as f:
        for line in f:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            
            chrom = tmp[0]
            if chrom != selected_chromosome:
                continue

            end = int(tmp[1])

            if chrom not in wps:
                wps[chrom] = []

            # Creating windows and storing their start positions for binary search
            start_positions = []
            for i in range(0, end, 120):
                window_dict = {
                    "start": i,
                    "end": i + 120,
                    "wps": 0
                }
                wps[chrom].append(window_dict)
                start_positions.append(i)

    return wps, start_positions

# Function to perform a fast intersection search
def find_window_for_position(wps_windows, start_positions, pos, end):
    idx = bisect.bisect_left(start_positions, pos)
        
    if idx < len(start_positions):
        window = wps_windows[idx]
        if window["start"] <= pos < window["end"]:
            return window
        if window["start"] <= end < window["end"]:
            return window
    return None


def windowed_protection_scores(sample_list, ann_dict, output_dir):
    """ """

    w = 80
    for sample in sample_list:
        bam_file = sample.bam
        window_start = ""
        WPS = {}

        # result = create_wps_windows(ann_dict["chromosomes"], "chr1")
        # WPS = result[0]
        # start_positions = result[1]

        cmd = f"samtools view {bam_file} chr1:1000000-1200000"
        with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True) as proc:
            for line in proc.stdout:
                fields = line.split("\t")
                chrom = fields[2]

                if len(fields) < 10:
                    continue
                mapqual = int(fields[4])  # Mapping quality
                if mapqual < 20:
                    continue
                pos = int(fields[3])
                sequence = fields[9]  # Sequence
                cigar = fields[5]  # CIGAR string

                # Get the actual read coordinates by adjusting for soft-clipping
                coords_dict = get_corrected_positions(cigar, pos, len(sequence))

                # window = find_window_for_position(WPS[chrom], start_positions, coords_dict["start"], coords_dict["end"])
                
                if not window_start:
                    window_start = coords_dict["start"]
                    window_end = window_start + w
                
                if coords_dict["start"] > window_end:
                    window_start = coords_dict["start"]
                    window_end = window_start + w

                window_coord = f"{chrom}\t{window_start}\t{window_end}"
                if not window_coord in WPS:
                    WPS[window_coord] = defaultdict(dict)
                    windows_dict = {
                        "spanning_reads": 0,
                        "endpoint_reads": 0
                    }
                
                if coords_dict["start"] <= window_start and coords_dict["end"] >= window_end:
                    windows_dict["spanning_reads"] += 1
                if (coords_dict["start"] >= window_start and coords_dict["start"] < window_end) or (coords_dict["end"] >= window_start and coords_dict["end"] < window_end):
                    windows_dict["endpoint_reads"] += 1
                
                WPS[window_coord] = windows_dict
                # print(line, coords_dict, window_coord, windows_dict)
    
        wps_txt = os.path.join(output_dir, f"{sample.name}.wps.tsv")
        wps_png = wps_txt.replace(".tsv", ".png")
        o = open(wps_txt, "w")
        o.write("chr\tstart\tend\tWPS\n")
        for coordinate in WPS:
            WPS[coordinate]["wps"] = WPS[coordinate]["spanning_reads"]-WPS[coordinate]["endpoint_reads"] 
            print(coordinate, WPS[coordinate])
            o.write(f'{coordinate}\t{WPS[coordinate]["wps"]}'+"\n")
        o.close()

        plot_wps(wps_txt, wps_png)

        sys.exit()

def plot_wps(input_txt, output_png):
    """ """

    df = pd.read_csv(input_txt, sep="\t")
    plt.figure(figsize=(20, 4))
    window_size = 15
    poly_order = 4
    y_smooth = savgol_filter(df["WPS"], window_size, poly_order)
    df["WPS"] = y_smooth

    sns.lineplot(data=df, x="start", y="WPS")

    # Set titles and labels
    plt.xlabel("Position (Mb)", fontsize=14)
    plt.ylabel("WPS", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_png)

    print(df)