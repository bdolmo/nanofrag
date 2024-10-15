import os
import subprocess
import re
from collections import defaultdict
from bx.intervals.intersection import Intersecter, Interval
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from modules.utils import get_chromosome_sizes
import gzip

import math

from scipy.signal import savgol_filter


def call_peaks(z_scores_file, peaks_bed):
    """ """

    first_chrom = ""
    first_pos = False
    first_start = ""
    with open(z_scores_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("chrom\t"):
                continue
            tmp = line.split("\t")
            chrom = tmp[0]
            pos = int(tmp[1])
            score = float(tmp[-1])

            if not first_pos:
                score_list = []
                positions_list = []
                first_pos = True
                first_chrom = chrom
                low_scores = 0
                total_length = 0
                first_start = pos

            if score >= -1 and first_pos:
                score_list.append(score)
                if score < 0:
                    low_scores+=1
                total_length+=1
                positions_list.append(pos)

            if score < -1 or low_scores >= 5:
                first_pos = False
                if total_length >= 80 and total_length < 250:
                    # found cluster
                    end = first_start + total_length
                    mean_score = np.mean(score_list)
                    print(chrom, first_start, end, mean_score, sep="\t")
                # first_start = 0
                # first_chrom = chrom
                # low_scores = 0
                
                
                


            # if chrom != first_chrom:
            #     first_chrom = chrom
            #     o.write(f"fixedStep chrom={chrom} start={start} step=1\n")
            # o.write(tmp[-1]+"\n")
    f.close()
    # o.close()


def create_wig(z_scores_file, wig_file, start=0):
    """ """

    o = open(wig_file, "w")

    first_chrom = ""

    with open(z_scores_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("chrom\t"):
                continue
            tmp = line.split("\t")
            chrom = tmp[0]
            if chrom != first_chrom:
                first_chrom = chrom
                o.write(f"fixedStep chrom={chrom} start={start} step=1\n")
            o.write(tmp[-1]+"\n")
    f.close()
    o.close()


def calculate_z_scores(wps_txt, sample_name, output_dir):
    """Calculates z-scores for 1kb windows directly from the WPS file."""

    # Open the WPS file for reading
    with open(wps_txt, "r") as f:
        lines = f.readlines()

    window_size = 500  # 1kb window
    current_window = []
    z_scores = []

    zscore_txt = os.path.join(output_dir, f"{sample_name}.z_scores.tsv")
    
    # Open the output file for writing z-scores
    with open(zscore_txt, "w") as out:
        out.write("chrom\tstart\tend\tWPS\tz_score\n")

        # Process the WPS values in 1kb windows
        for i in range(0, len(lines), window_size):
            current_window = []

            # Collect WPS values for the current 1kb window
            for j in range(i, min(i + window_size, len(lines))):
                line = lines[j].strip()
                if not line:
                    continue

                try:
                    chrom, start, end, wps = line.split("\t")
                    wps = float(wps)  # Convert WPS to a float value
                    current_window.append((chrom, int(start), int(end), wps))
                except ValueError:
                    # Handle malformed lines
                    continue

            # Calculate mean and standard deviation for the current window
            if current_window:
                mean_wps = sum([x[3] for x in current_window]) / len(current_window)
                variance = sum([(x[3] - mean_wps) ** 2 for x in current_window]) / len(current_window)
                std_wps = math.sqrt(variance)

                w_zscores = []

                # Compute z-scores for the current window
                for chrom, start, end, wps in current_window:
                    if std_wps > 0:
                        z_score = (wps - mean_wps) / std_wps
                    else:
                        z_score = 0  # If std is 0, set z-score to 0
                    w_zscores.append(z_score)
                if len(w_zscores) >= 100:
                    z_scores = savgol_filter(w_zscores, 100, 2)
                else:
                    z_scores = w_zscores
                idx = 0
                for chrom, start, end, wps in current_window:
                    z_score = z_scores[idx]
                    # Write the result to the output file
                    out.write(f"{chrom}\t{start}\t{end}\t{wps}\t{z_score}\n")
                    idx+=1


def plot_z_scores_txt(sample_name, zscore_txt, output_png):
    """Plots the z-scores directly from the z-score file."""

    # Prepare lists to hold the data for plotting
    positions = []
    z_scores = []

    # Read the z-score file and extract positions and z-scores
    with open(zscore_txt, "r") as f:
        for line in f.readlines()[1:]:  # Skip header
            try:
                _, start, _, _, z_score = line.strip().split("\t")
                positions.append(int(start))
                z_scores.append(float(z_score))
            except ValueError:
                # Skip lines that cannot be parsed
                continue
    # z_scores = savgol_filter(z_scores, 200, 2)

    plt.figure(figsize=(20, 4))
    plt.plot(positions, z_scores, linestyle="-")
    plt.title(sample_name, fontsize=18)
    plt.xlabel("Genomic Position", fontsize=15)
    plt.ylabel("WPS z-score", fontsize=15)
    plt.ylim(-5, 5)
    plt.savefig(output_png)
    plt.close()


def get_corrected_positions(cigar, pos, sequence_length):
    """
    Calculate the actual fragment size by subtracting soft-clipped portions from the ends.
    """
    softclip_start = re.match(r"^(\d+)S", cigar)
    softclip_end = re.search(r"(\d+)S$", cigar)

    start = pos
    end = pos + sequence_length

    if softclip_start:
        s_start = softclip_start.group(1)
        sequence_length -= int(s_start)
        start += int(s_start)
    if softclip_end:
        s_end = softclip_end.group(1)
        sequence_length -= int(s_end)
        end -= int(s_end)

    return {"start": start, "end": end}


def windowed_protection_scores(sample_list, ann_dict, bin_dict, output_dir, region, protection=120):
    """Calculate WPS for each base with a 120 bp window and plot the Z-scores."""

    chrom_sizes_dict = get_chromosome_sizes(ann_dict["chromosomes"])

    if region:
        tmp_region = re.split(r"[:-]", region)
        r_chrom = tmp_region[0]
        r_pos = int(tmp_region[1])
        r_end = int(tmp_region[2])

    window_size = protection
    half_window = window_size // 2
    protection_half = protection // 2

    for sample in sample_list:
        bam_file = sample.bam
        posRange = defaultdict(lambda: [0, 0, 0])  # First for coverage, second for starts/ends, third for WPS

        filteredReads = Intersecter()  # Use for interval intersection

        if region:
            cmd = f"samtools view {bam_file} {r_chrom}:{r_pos}-{r_end}"
        else:
            cmd = f"samtools view {bam_file} chr12:34443233-34453733"

        wps_txt = os.path.join(output_dir, f"{sample.name}.windowed_wps.tsv")
        o = open(wps_txt, "w")
        o.write("chr\tstart\tend\tWPS\n")

        with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True) as proc:
            for line in proc.stdout:
                fields = line.split("\t")
                if len(fields) < 10:
                    continue

                chrom = fields[2]
                mapqual = int(fields[4])
                if mapqual < 10:
                    continue

                pos = int(fields[3])
                sequence = fields[9]
                cigar = fields[5]

                coords_dict = get_corrected_positions(cigar, pos, len(sequence))

                rstart = coords_dict["start"]
                rend = coords_dict["end"]

                fragment_size = rend-rstart
                if fragment_size < 120:
                    continue
                if fragment_size > 180:
                    continue
                # Add interval to filteredReads for protection calculation later
                filteredReads.add_interval(Interval(rstart, rend))

                # Process positions and store start/end info in posRange
                for i in range(rstart, rend + 1):
                    if r_pos <= i <= r_end:
                        posRange[i][0] += 1  # Total reads covering the position
                if r_pos <= rstart <= r_end:
                    posRange[rstart][1] += 1  # Start of a read
                if r_pos <= rend <= r_end:
                    posRange[rend][1] += 1  # End of a read

        # Calculate WPS using the protection windows
        for pos in range(r_pos, r_end + 1):
            rstart_window = pos - protection_half
            rend_window = pos + protection_half
            gcount, bcount = 0, 0

            # Check for intervals in the window and calculate WPS
            for read in filteredReads.find(rstart_window, rend_window):
                if read.start > rstart_window or read.end < rend_window:
                    bcount += 1
                else:
                    gcount += 1

            posRange[pos][2] += gcount - bcount  # WPS calculation

            # Write WPS to the output file
            o.write(f"{chrom}\t{pos}\t{pos}\t{posRange[pos][2]}\n")

        o.close()

        # Calculate z-scores and plot them
        calculate_z_scores(wps_txt, sample.name, output_dir)
        zscore_txt = os.path.join(output_dir, f"{sample.name}.z_scores.tsv")
        output_png = zscore_txt.replace(".tsv", ".png")

        plot_z_scores_txt(sample.name, zscore_txt, output_png)

        wig_file = zscore_txt.replace(".tsv", ".wig")
        create_wig(zscore_txt, wig_file, start=r_pos)

        peaks_bed = zscore_txt.replace(".tsv", ".peaks.bed")
        call_peaks(zscore_txt, peaks_bed)


def plot_z_scores(df, output_png):
    """Plots the z-scores from the DataFrame."""
    plt.figure(figsize=(20, 4))
    sns.lineplot(data=df, x="start", y="z_score")
    plt.title("Z-Scores of Windowed Protection Score (WPS)")
    plt.xlabel("Genomic Position")
    plt.ylabel("Z-Score")
    plt.savefig(output_png)
    plt.close()
