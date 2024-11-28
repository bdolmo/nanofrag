import subprocess
import os
import sys
import gzip
import shutil
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from modules.utils import tissues_dict



def ensure_gzipped(file_path):
    """
    Checks if the file is gzipped; if not, gzips it in place using bgzip.

    Args:
        file_path (str): Path to the file to check and gzip if needed.

    Returns:
        str: Path to the gzipped file.
    """
    if not file_path.endswith(".gz"):
        gzipped_path = file_path + ".gz"
        # Run bgzip to compress the file
        subprocess.run(["bgzip", "-c", file_path], stdout=open(gzipped_path, 'wb'))
        print(f"Gzipped {file_path} to {gzipped_path}")
        return gzipped_path
    return file_path

def index_with_tabix(input_gz):
    """ """
    cmd = [
        "tabix", "-p", "bed", input_gz
    ]
    try:
        subprocess.run(cmd, check=True)
        print("tabix index completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during tabix index: {e}")



def run_modkit_dmr_pair(bin_dict, output_dir, norm_pileup, tumor_pileup, regions, ref, dmr_result=None, base="C", threads=1, log_filepath="dmr.log", use_index_a=False, use_index_b=False):
    """
    Run the modkit dmr pair command for differential methylation analysis.

    Args:
        norm_pileup (str): Path to the normal pileup file (gzipped if not already).
        tumor_pileup (str): Path to the tumor pileup file (gzipped if not already).
        regions (str): Path to the BED file with regions of interest (e.g., CpG islands).
        ref (str): Path to the reference FASTA file.
        dmr_result (str): Path to the output BED file for differential methylation results. If None, outputs to stdout.
        base (str): Base to be used for the modification (default is 'C').
        threads (int): Number of threads to use (default is 1).
        log_filepath (str): Path to the log file.
        use_index_a (bool): Whether to use the index file for norm_pileup.
        use_index_b (bool): Whether to use the index file for tumor_pileup.

    Returns:
        None
    """
    # Ensure norm_pileup and tumor_pileup are gzipped

    norm_pileup_gz = f"{norm_pileup}.gz"
    tumor_pileup_gz = f"{tumor_pileup}.gz"

    output_dmr = os.path.join(output_dir, os.path.basename(tumor_pileup).replace(".bed", ".dmr.bed"))

    # if not os.path.isfile(norm_pileup_gz):
    norm_pileup_gz = ensure_gzipped(norm_pileup)

    index_with_tabix(norm_pileup_gz)

    # if not os.path.isfile(tumor_pileup_gz):
    tumor_pileup_gz = ensure_gzipped(tumor_pileup)
    index_with_tabix(tumor_pileup_gz)

    # Build the command
    cmd = [
        bin_dict["modkit"], "dmr", "pair",
        "-a", norm_pileup_gz
    ]

    # Optionally add index for norm_pileup if use_index_a is True
    if use_index_a:
        cmd.extend(["--index-a", norm_pileup_gz + ".tbi"])

    # Add tumor pileup and optional index
    cmd.extend(["-b", tumor_pileup_gz])
    if use_index_b:
        cmd.extend(["--index-b", tumor_pileup_gz + ".tbi"])

    # Add output, regions, reference, base, threads, and log file path
    if dmr_result:
        cmd.extend(["-o", dmr_result])
    cmd.extend([
        "-r", regions,
        "-o", output_dmr,
        "--ref", ref,
        "--base", base,
        "--threads", "4",
        "--log-filepath", log_filepath
    ])
    if not os.path.isfile(output_dmr):
        # Execute the command
        try:
            subprocess.run(cmd, check=True)
            print("modkit dmr pair completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error occurred during modkit dmr pair: {e}")



def run_modkit_pileup(bin_dict, reference_fasta, threads, input_bam, pileup_bed, ignore='h', combine_strands=True):
    """
    Run the modkit pileup command with specified options.
    """
    # Build the base command
    cmd = [
        bin_dict["modkit"], "pileup", input_bam, pileup_bed, "--cpg", "--ref", reference_fasta, "--threads", str(threads),
        "--ignore", ignore
    ]

    # Add the --combine-strands option if requested
    if combine_strands:
        cmd.append("--combine-strands")

    if not os.path.isfile(pileup_bed):
        msg = f" INFO: Performing methylation analysis for {input_bam}"
        print(msg)
        # Execute the command
        try:
            subprocess.run(cmd, check=True)
            msg = " INFO: modkit pileup completed successfully"
            print(msg)
        except subprocess.CalledProcessError as e:
            msg = f"Error occurred during modkit pileup: {e}"
            print(msg)
    else:
        msg = f" INFO: Skipping methylation analysis for {input_bam}"
        print(msg)

def plot_deconvolution(sample_name, input_file, output_png):
    """
    Plots a stacked bar chart for cell type proportions from a text file, in descending order.
    """
    # Read the data from the text file
    df = pd.read_csv(input_file, sep='\t')  # Adjust delimiter if necessary
    
    # Filter to show only non-zero proportions and sort by ascending order for stacking
    df['proportion'] = df['proportion']*100


    df_filtered = df[df['proportion'] > 0].sort_values(by="proportion", ascending=True)
    
    # Calculate cumulative proportions for stacking
    df_filtered['cumulative'] = df_filtered['proportion'].cumsum() - df_filtered['proportion']
    # Use the Spectral palette

    # palette = sns.color_palette()
    # colors = sns.color_palette("tab20c")
    # # print(colors)
    # palette = sns.color_palette(tab20c)
    # Initialize plot
    # colors = 

    plt.figure(figsize=(10, 6))
    sns.set_theme(style="ticks")
    
    # Plot each bar segment individually, stacking in descending order
    for index, row in df_filtered.iterrows():
        plt.bar(
            sample_name, 
            row['proportion'], 
            bottom=row['cumulative'], 
            color=tissues_dict.get(row['cell_type'], "#333333"),
            label=f"{row['cell_type']} ({str(round(float(row['proportion']), 2))}%)"
        )
    
    # Annotate top 3 stacked elements
    top_elements = df_filtered.nlargest(3, 'proportion')
    # for i, (index, row) in enumerate(df_filtered.iterrows()):
    #     plt.text(
    #         x=sample_name, 
    #         y=row['cumulative'] + row['proportion'] / 2, 
    #         s=f"{row['proportion']:.2f}", 
    #         ha='right', 
    #         va='right', 
    #         color="black"
    #     )

    sns.despine(left=True, bottom=True)

    # Set title and labels
    sample_name = sample_name.replace(".methylation", "").replace(".sorted.aligned", "")
    plt.title(sample_name)
    # plt.xlabel("")
    plt.ylabel("Proportion (%)")
    
    # Create a unique legend
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    
    # Invert the order of the legend
    inverted_labels = list(by_label.keys())[::-1]
    inverted_handles = list(by_label.values())[::-1]
    
    plt.legend(inverted_handles, inverted_labels, title="Cell type", 
        bbox_to_anchor=(1.05, 1), loc="upper left")
    
    plt.tight_layout()
    plt.savefig(output_png)
    plt.close()


def deconvolve_tissues_with_nanomix(bin_dict, sample, docker_output, output_dir, methyl_bed, atlas_bed):
    """
    Run the nanomix deconvolute command using Docker with specified options.
    """
    # Adapt methyl_bed to fit nanomix input format
    methyl_nanomix_bed = methyl_bed.replace(".pileup.bed", ".nanomix.bed")

    if not os.path.isfile(methyl_nanomix_bed):
        with open(methyl_nanomix_bed, "w") as o:
            o.write("chr\tstart\tend\ttotal_calls\tmodified_calls\n")
            with open(methyl_bed) as f:
                for line in f:
                    line = line.rstrip("\n")
                    tmp = line.split("\t")
                    o.write(f"{tmp[0]}\t{tmp[1]}\t{tmp[2]}\t{tmp[9]}\t{tmp[11]}\n")

    # Define paths for Docker command
    atlas_path = os.path.abspath(atlas_bed)
    output_path = os.path.abspath(output_dir)
    methylation_input = os.path.abspath(methyl_nanomix_bed)
    output_file = os.path.join(output_path, f"{sample.name}_tissue-proportions_nanomix_5hmC.txt")

    # Construct the Docker command
    nanomix_command = [
        "docker", "run", "--rm", "-v", f"{atlas_path}:/data/methylation_atlas",
        "-v", f"{docker_output}/METHYLATION:/data/output_dir", "bdolmo/nanomix:1.0.0",
        "/bin/bash", "-c",
        f"nanomix deconvolute -a /data/methylation_atlas /data/output_dir/{os.path.basename(methyl_nanomix_bed)} > /data/output_dir/{sample.name}_tissue-proportions_nanomix_5hmC.txt"
    ]
    if not os.path.isfile(output_file):
        # Run the command
        try:
            subprocess.run(nanomix_command, check=True)
            print(f"Nanomix deconvolution completed. Output saved to {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error running nanomix deconvolution: {e}")
    output_png = os.path.join(output_dir, f"{sample.name}.nanomix.deconvolution.png")
    plot_deconvolution(sample.name, output_file, output_png)


def run_methylation_analysis(sample_list, ann_dict, bin_dict, threads, reference_fasta, docker_output, output_dir):
    """ 
    """
    methylation_folder = os.path.join(output_dir, "METHYLATION")
    if not os.path.isdir(methylation_folder):
        os.mkdir(methylation_folder)

    tumor_samples = []
    normal_samples = []
    for sample in sample_list:       
        pileup_bed = os.path.join(methylation_folder, f"{sample.name}.pileup.bed")
        sample.add("methylation_pileup", pileup_bed)
        run_modkit_pileup(bin_dict, reference_fasta, threads, sample.bam, pileup_bed)

        deconvolve_tissues_with_nanomix(bin_dict, sample, docker_output, methylation_folder, pileup_bed, ann_dict["nanomix_atlas"])
        # sys.exit()

        if sample.origin == "tumor":
            tumor_samples.append(sample)
        else:
            normal_samples.append(sample)
    jdx = 0
    for idx, tumor_sample in enumerate(tumor_samples):
        if jdx < len(normal_samples):
            normal_sample = normal_samples[jdx]
        jdx+=1
        t_pileup = tumor_sample.methylation_pileup
        n_pileup = normal_sample.methylation_pileup

        run_modkit_dmr_pair(bin_dict,methylation_folder, n_pileup, t_pileup, ann_dict["cpg_islands"], 
            reference_fasta, log_filepath=f"{methylation_folder}/dmr.log")

    return sample_list



# Example usage:
# run_modkit_pileup("reference.fasta")
