import subprocess
import os
import sys


import subprocess
import os
import sys
import gzip
import shutil

def ensure_gzipped(file_path):
    """
    Checks if the file is gzipped; if not, gzips it in place.

    Args:
        file_path (str): Path to the file to check and gzip if needed.

    Returns:
        str: Path to the gzipped file.
    """
    if not file_path.endswith(".gz"):
        gzipped_path = file_path + ".gz"
        with open(file_path, 'rb') as f_in, gzip.open(gzipped_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f"Gzipped {file_path} to {gzipped_path}")
        return gzipped_path
    return file_path

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

    if not os.path.isfile(norm_pileup_gz):
        norm_pileup_gz = ensure_gzipped(norm_pileup)
    if not os.path.isfile(tumor_pileup_gz):
        tumor_pileup_gz = ensure_gzipped(tumor_pileup)

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
        # Execute the command
        try:
            subprocess.run(cmd, check=True)
            msg = " INFO: modkit pileup completed successfully"
            print(msg)
        except subprocess.CalledProcessError as e:
            msg = f"Error occurred during modkit pileup: {e}"
            print(msg)



def run_methylation_analysis(sample_list, ann_dict, bin_dict, threads, reference_fasta, output_dir):
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
        if sample.origin == "tumor":
            tumor_samples.append(sample)
        else:
            normal_samples.append(sample)

    for idx, tumor_sample in enumerate(tumor_samples):
        normal_sample = normal_samples[idx]

        t_pileup = tumor_sample.methylation_pileup
        n_pileup = normal_sample.methylation_pileup

        run_modkit_dmr_pair(bin_dict,methylation_folder, n_pileup, t_pileup, ann_dict["cpg_islands"], reference_fasta, log_filepath=f"{methylation_folder}/dmr.log")


    return sample_list



# Example usage:
# run_modkit_pileup("reference.fasta")
