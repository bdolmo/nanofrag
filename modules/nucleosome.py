import os
import subprocess
import re
from collections import defaultdict
from bx.intervals.intersection import Intersecter, Interval
import bisect

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from modules.utils import get_chromosome_sizes, bed_to_list
import gzip
import math
import sys
import multiprocessing
from scipy.signal import savgol_filter
import glob
import natsort



#########################################################

def preprocess_tss_bed(tss_bed, ann_dict):
    """ """
    tss_1kb = tss_bed.replace(".bed", "1kb.bed")
    o = open(tss_1kb, "w")

    # First, create 1kb windows from tss
    with open(tss_bed) as f:
        for line in f:
            if line.startswith("chromosome"):
                continue
            line = line.rstrip("\n")
            tmp = line.split("\t")
            size = int(tmp[2])-int(tmp[1])

            if size < 1000:
                tosum = int(round((1000-size)/2))
                prev_start = tmp[1]
                prev_end = tmp[2]
                tmp[1] = str (int(tmp[1])- tosum)
                tmp[2] = str (int(tmp[2])+ tosum)
                if (int(tmp[2])-int(tmp[1])) < 1000:
                    tmp[2] = str(int(tmp[2])+1)

                if int(tmp[1]) < 0:
                    tmp[1] = prev_start
                if int(tmp[2]) < 0:
                    tmp[2] = prev_end                
                line = "\t".join(tmp)
            o.write(line+"\n")

    f.close()
    o.close()

    tss_1kb_clean = tss_1kb.replace(".bed", ".clean.bed")

    cmd = f'bedtools  intersect -a {tss_1kb} -b {ann_dict["blacklist"]} -v > {tss_1kb_clean}'
    subprocess.run(cmd, shell=True, check=True)

    return tss_1kb_clean

def plot_tss_coverage(tss_cov_txt):
    """ """

    df = pd.read_csv(tss_cov_txt, sep="\t", names=["Position", "Coverage"], header=None)
    plt.figure(figsize=(20, 6))
    df["Coverage"]=(df["Coverage"]-df["Coverage"].mean())/df["Coverage"].std()


    # sns.lineplot(data=df, x="Position", y="Coverage", linewidth = 0.1)
    sns.scatterplot(data=df, x="Position", y="Coverage", size = 0.01)

    #     sns.histplot(result[col], bins=8000, label=col, kde=False, color=color, alpha=0.5)

    # # Set titles and labels
    plt.xlabel("Position arround TSS", fontsize=14)
    plt.ylabel("Coverage", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    cov_png = tss_cov_txt.replace(".txt", ".png")
    plt.savefig(cov_png)
    plt.close()


def run_tss_analysis(sample_list, ann_dict, output_dir):
    """ """
    tss_bed = preprocess_tss_bed(ann_dict["tss"], ann_dict)

    regions = []
    with open(tss_bed) as f:
        for line in f:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            size = int(tmp[2])-int(tmp[1])
            if size != 1000:
                continue

            region = f"{tmp[0]}:{tmp[1]}-{tmp[2]}"
            regions.append(region)
    f.close()


    for sample in sample_list:


        cov_dict = {}
        for i in range(0, 1001):
            cov_dict[i] = []

        # print(cov_dict)
        # sys.exit()

        # cov_array = [[0] * 1000 for i in range(1)]
        # # print(cov_array)ç
        # for r in cov_array:
        #     print(r)
        # sys.exit()


        for n,region in enumerate(regions):

            if n > 2000:
                break
            # cmd = f'samtools depth -aa -b {tss_bed} {sample.bam}'
            cmd = f'samtools depth -aa -r {region} {sample.bam}'
            with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True) as proc:
                for idx,line in enumerate(proc.stdout):
                    line = line.rstrip("\n")
                    tmp = line.split("\t")
                    cov_dict[idx].append(int(tmp[-1]))

        plotdata_txt = os.path.join(output_dir, "NUCLEOSOME", f"{sample.name}.tss.txt")
        o = open(plotdata_txt, "w")
        for idx in cov_dict:
                
            mean_cov = np.mean(cov_dict[idx])
            print(mean_cov)
            # if idx <= 500:
            idx = idx-500
            o.write(f"{str(idx)}\t{mean_cov}\n")
            # for value in cov_dict[idx]:
            #     idx = idx-500

            #     o.write(f"{str(idx)}\t{value}\n")
        o.close()
        plot_tss_coverage(plotdata_txt)
             



    # now remove low quality regions that overlap the blacklist


#########################################################
    # Let's try a profiling based on available nucleosome positions


# def find_window_for_read(chrom, read_start, windows, window_starts, window_ends_by_chrom):
#     """
#     Find the window that a given read belongs to using binary search for efficiency.
#     """
#     if chrom not in window_starts:
#         return None

#     idx = bisect.bisect_right(window_starts[chrom], read_start) - 1

#     if idx >= 0 and window_starts[chrom][idx] <= read_start <= window_ends_by_chrom[chrom][idx]:
#         return windows[chrom][idx]

#     return None

# def find_window_for_read(chrom, read_start, read_end, windows, window_starts, window_ends_by_chrom):
#     """
#     Find the window that a given read belongs to using binary search for efficiency.
#     """
#     if chrom not in window_starts:
#         return None

#     # Use binary search to find the appropriate window index
#     idx = bisect.bisect_right(window_starts[chrom], read_start) - 1

#     # Check if the read falls within the bounds of the identified window
#     if idx >= 0 and window_starts[chrom][idx] <= read_start <= window_ends_by_chrom[chrom][idx] and window_starts[chrom][idx] <= read_end <= window_ends_by_chrom[chrom][idx]:
#         return windows[chrom][idx]

#     return None



# def plot_nucleosome_enrichment(input_txt, sample_name):
#     """ """
#     df = pd.read_csv(input_txt, names=[sample_name], header=None)
#     print(df)
#     plt.figure(figsize=(10, 6))

#     sns.histplot(df[sample_name],bins=600,kde=False)

#     #     sns.histplot(result[col], bins=8000, label=col, kde=False, color=color, alpha=0.5)

#     # # Set titles and labels
#     plt.xlabel("Relative position to nucleosome center", fontsize=14)
#     plt.ylabel("Frequency", fontsize=14)
#     plt.xticks(fontsize=12)
#     plt.yticks(fontsize=12)
#     relative_pos_png = input_txt.replace(".txt", ".png")
#     plt.savefig(relative_pos_png)
#     plt.close()
# def run_nucleosome_profiling(sample_list, ann_dict, output_dir):
#     """ """

#     window_starts = {}
#     window_ends_by_chrom = {}
#     windows_by_chrom = {}
    
#     regions = bed_to_list(ann_dict["nucleosomes"], gzipped=True)
#     n_regions = 0
#     tmp_nucleosomes = ann_dict["nucleosomes"].replace(".bed.gz", ".tmp.bed")
#     o = open(tmp_nucleosomes, "w")
#     for r in regions:
#         n_regions +=1
#         if n_regions > 500000:
#             break
#         region = r.split("\t")
#         chrom = region[0]
#         start = int(region[1])-300
#         end = int(region[2])+300
#         if not chrom in windows_by_chrom:
#             window_starts[chrom] = []
#             window_ends_by_chrom[chrom] = []
#             windows_by_chrom[chrom] = []
#         window_starts[chrom].append(start)
#         window_ends_by_chrom[chrom].append(end)
#         windows_by_chrom[chrom].append((chrom, start, end))
#         o.write(f"{chrom}\t{str(start)}\t{str(end)}\n")

#     o.close()


#     n_reads = 0
#     for sample in sample_list:
        
#         distances_to_center = []
#         n_reads = 0

#         cmd = f"samtools view {sample.bam} -L {tmp_nucleosomes}"

#         with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True) as proc:
#             for line in proc.stdout:
#                 fields = line.split("\t")
#                 if len(fields) < 10:
#                     continue

#                 chrom = fields[2]
#                 mapqual = int(fields[4])
#                 if mapqual < 50:
#                     continue
#                 n_reads+=1

#                 pos = int(fields[3])
#                 sequence = fields[9]
#                 cigar = fields[5]
#                 if n_reads > 200000:
#                     break


#                 coords_dict = get_corrected_positions(cigar, pos, len(sequence))
#                 read_start = int(coords_dict["start"])
#                 read_end = int(coords_dict["end"])

#                 window = find_window_for_read(chrom, read_start, read_end, windows_by_chrom, window_starts, window_ends_by_chrom)
#                 if window:
#                     has_supplementary = False
#                     for field in fields:    
#                         if field.startswith("SA:"):
#                             has_supplementary = True
#                     if not has_supplementary:

#                         n_reads+=1

#                         print(n_reads, pos,  line, window, type(window))

#                         distance_to_center = read_start-int(window[1]+300)
#                         distances_to_center.append(distance_to_center)

#     output_file = os.path.join(output_dir, "NUCLEOSOME", f"{sample.name}.dist.to.center.txt")    
#     o = open(output_file, "w")
#     for dist in distances_to_center:
#         o.write(str(dist)+"\n")
#     o.close()

#     plot_nucleosome_enrichment(output_file, sample.name)

#     sys.exit()


#         # for idx,region in enumerate(regions):

#         #     print(region)

#         #     tmp_r = region.split("\t")
#         #     r_chr = tmp_r[0]
#         #     r_pos = str(int(tmp_r[1])-75)
#         #     r_end = str(int(tmp_r[2])+75)
            
#         #     region_str = f"{r_chr}:{r_pos}-{r_end}"
#         #     cmd = f"samtools view {sample.bam} {region_str}"

#         #     with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True) as proc:
#         #         for line in proc.stdout:
#         #             fields = line.split("\t")
#         #             if len(fields) < 10:
#         #                 continue

#         #             chrom = fields[2]
#         #             mapqual = int(fields[4])
#         #             if mapqual < 10:
#         #                 continue
#         #             n_reads+=1
#         #             print(n_reads)
#         #             if n_reads > 10000:
#         #                 break

#         #             pos = int(fields[3])
#         #             sequence = fields[9]
#         #             cigar = fields[5]

#         #             coords_dict = get_corrected_positions(cigar, pos, len(sequence))

#         #             read_start = int(coords_dict["start"])
#         #             read_end = int(coords_dict["end"])

#         #             distance_to_center = read_start-int(r_pos)
#         #             distances_to_center.append(distance_to_center)
    

#         # o = open(output_file, "w")
#         # for dist in distances_to_center:
#         #     o.write(dist+"\n")
#         # o.close()




#########################################################33


def call_peaks(z_scores_file, peaks_bed):
    """ """

    first_chrom = ""
    first_pos = False
    first_start = ""


    o = open(peaks_bed, "w")

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
                    o.write(f"{chrom}\t{first_start}\t{end}\t{mean_score}\n")
                # first_start = 0
                # first_chrom = chrom
                # low_scores = 0           
            # if chrom != first_chrom:
            #     first_chrom = chrom
            #     o.write(f"fixedStep chrom={chrom} start={start} step=1\n")
            # o.write(tmp[-1]+"\n")
    f.close()
    o.close()


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
    """
        Calculates z-scores for 1kb windows directly from the WPS file.
    """

    # Open the WPS file for reading
    with open(wps_txt, "r") as f:
        lines = f.readlines()

    window_size = 1000  # 1kb window
    current_window = []
    z_scores = []

    zscore_txt = os.path.join(output_dir, f"{sample_name}.z_scores.tsv")
    
    # Open the output file for writing z-scores
    with open(zscore_txt, "w") as out:
        out.write("chrom\tstart\tWPS\tz_score\n")

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
    """
        Plots the z-scores directly from the z-score file.
    """

    # Prepare lists to hold the data for plotting
    positions = []
    z_scores = []

    # Read the z-score file and extract positions and z-scores
    with open(zscore_txt, "r") as f:
        for line in f.readlines()[1:]:  # Skip header

            line = line.rstrip("\n")
            tmp = line.split("\t")
            positions.append(int(tmp[1]))
            z_scores.append(float(tmp[2]))

    z_scores = savgol_filter(z_scores, 50, 3)

    plt.figure(figsize=(20, 4))
    plt.plot(positions, z_scores, linestyle="-")
    plt.title(sample_name, fontsize=18)
    plt.xlabel("Genomic Position", fontsize=15)
    plt.ylabel("WPS z-score", fontsize=15)
    # plt.ylim(-5, 5)
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

def run_wps_analysis(sample_list, ann_dict, bin_dict, output_dir, protection=120):
    """ """

    nucleosome_folder = os.path.join(output_dir, "NUCLEOSOME")
    if not os.path.isdir(nucleosome_folder):
        os.mkdir(nucleosome_folder)

    regions = []
    windows_bed = os.path.join(output_dir, "FRAGMENTATION", "windows_5000000.bed")
    with open(windows_bed) as f:
        for line in f:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            regions.append(f"{tmp[0]}:{tmp[1]}-{tmp[2]}")


    regions = []
    regions.append("chr12:34267065-34294065")

    for sample in sample_list:

        threads = 10
        with multiprocessing.Pool(threads) as pool:
            # Distribute the tasks across the pool of workers
            for region in regions:
                name_region = region.replace(":", "_").replace("-", "_")
                # wps_bed = os.path.join(nucleosome_folder, f"{sample.name}.{name_region}.windowed_wps.bed")
                # o = open(wps_bed, "w")
                # o.write("chr\tstart\tend\tWPS\n")
                # o.close()
                results = pool.starmap(windowed_protection_scores, [(sample, ann_dict, bin_dict, os.path.join(nucleosome_folder, 
                    f"{sample.name}.{idx}.windowed_wps.bed"), region, 120) for idx,region in enumerate(regions)] )
        sample_wps_data = f"{nucleosome_folder}/{sample.name}.master.wps.bed"

        o = open(sample_wps_data, "w")
        files = glob.glob(f"{nucleosome_folder}/*windowed_wps.bed")
        files = natsort.natsorted(files)
        for file in files:
            with open(file) as f:
                for line in f:
                    # line = line.rstrip("\n")
                    # tmp = line.split("\t")
                    o.write(line)
            print(file)
        o.close()

        sample_peaks = f"{nucleosome_folder}/{sample.name}.wps.peaks.bed"
        call_peaks(sample_wps_data, sample_peaks)


        output_png = f"{nucleosome_folder}/{sample.name}.wps.png"
        plot_z_scores_txt(sample.name, sample_wps_data, output_png)

        # plot_peak_info(sample_peaks)

        sys.exit()



def windowed_protection_scores(sample, ann_dict, bin_dict, wps_bed, region, protection=120):
    """
        Calculate WPS for each base with a 120 bp window and plot the Z-scores.
    """

    chrom_sizes_dict = get_chromosome_sizes(ann_dict["chromosomes"])

    if region:
        tmp_region = re.split(r"[:-]", region)
        r_chrom = tmp_region[0]
        r_pos = int(tmp_region[1])
        r_end = int(tmp_region[2])
        msg = f" INFO: Calculating WPS on {region}"
        print(msg)

    window_size = protection
    half_window = window_size // 2
    protection_half = protection // 2

    bam_file = sample.bam
    posRange = defaultdict(lambda: [0, 0, 0, 0])  # First for coverage, second for starts/ends, third for WPS

    filteredReads = Intersecter()  # Use for interval intersection

    if region:
        cmd = f"samtools view {bam_file} {r_chrom}:{r_pos}-{r_end}"
    else:
        cmd = f"samtools view {bam_file}"


    wps_raw = []
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
            # if fragment_size < 120:
            #     continue
            # if fragment_size > 180:
            #     continue
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


    o = open(wps_bed, "w")
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
        wps_raw.append(posRange[pos][2])

    mean_wps_raw = np.mean(wps_raw)
    std_wps_raw = np.std(wps_raw)
    for pos in range(r_pos, r_end+1):
        zscore = (posRange[pos][2]-mean_wps_raw)/std_wps_raw

        # Write WPS to the output file
        # o.write(f"{r_chrom}\t{pos}\t{zscore}\n")

        o.write(f"{r_chrom}\t{pos}\t{posRange[pos][0]}\n")


        # o.write(f"{r_chrom}\t{pos}\t{posRange[pos][2]}\n")

    o.close()






    return True

    # calculate_z_scores(wps_txt, sample.name, output_dir)
    # zscore_txt = os.path.join(output_dir, f"{sample.name}.z_scores.tsv")
    # output_png = zscore_txt.replace(".tsv", ".png")

    # plot_z_scores_txt(sample.name, zscore_txt, output_png)

    # wig_file = zscore_txt.replace(".tsv", ".wig")
    # create_wig(zscore_txt, wig_file, start=r_pos)

    # peaks_bed = zscore_txt.replace(".tsv", ".peaks.bed")
    # call_peaks(zscore_txt, peaks_bed)


def plot_z_scores(df, output_png):
    """Plots the z-scores from the DataFrame."""
    plt.figure(figsize=(20, 4))
    sns.lineplot(data=df, x="start", y="z_score")
    plt.title("Z-Scores of Windowed Protection Score (WPS)")
    plt.xlabel("Genomic Position")
    plt.ylabel("Z-Score")
    plt.savefig(output_png)
    plt.close()
