import os
import sys
import subprocess

def run_clairs(output_dir, tumor, normal, ref_fasta, threads, platform="ont_r10_dorado_sup_5khz_ssrs"):
    """ """

    platforms = ["ont_r10_dorado_sup_4khz",
        "ont_r10_dorado_sup_5khz_ssrs",
        "ont_r10_dorado_sup_5khz", 
        "ont_r10_guppy", "ont_r9_guppy", 
        "ilmn", "hifi_sequel2", "hifi_revio"]

    
    tumor_dir = os.path.dirname(tumor.bam)
    normal_dir = os.path.dirname(normal.bam)
    ref_dir = os.path.dirname(ref_fasta)

    tumor_bam = os.path.basename(tumor.bam)
    normal_bam = os.path.basename(normal.bam)
    ref_fasta_name = os.path.basename(ref_fasta)
    tumor_name = tumor_bam.replace(".bam", "")
    normal_name = normal_bam.replace(".bam", "")

    # command = [
    #     "docker", "run", "-it",
    #     "-v", f"{tumor_dir}:{tumor_dir}",
    #     "-v", f"{normal_dir}:{normal_dir}",
    #     "-v", f"{output_dir}:{output_dir}",
    #     "-v", f"{ref_dir}:{ref_dir}",
    #     "hkubal/clairs:latest",
    #     "/opt/bin/run_clairs",
    #     "--tumor_bam_fn", f"{tumor_dir}/{tumor_bam}",
    #     "--normal_bam_fn", f"{normal_dir}/{normal_bam}",
    #     "--ref_fn", f"{ref_dir}/{ref_fasta_name}",
    #     "--threads", str(threads),
    #     "--platform", platform,
    #     "--output_dir", output_dir,
    #     "--min_coverage 2",
    #     f"--output_prefix {tumor_name}_{normal_name}"
    # ]


    command = [
        "docker", "run", "-it",
        "-v", f"{tumor_dir}:{tumor_dir}",
        "-v", f"{normal_dir}:{normal_dir}",
        "-v", f"{output_dir}:{output_dir}",
        "-v", f"{ref_dir}:{ref_dir}",
        "hkubal/clairs:latest",
        "/opt/bin/run_clairs",
        "--tumor_bam_fn", f"{tumor_dir}/{tumor_bam}",
        "--normal_bam_fn", f"{normal_dir}/{normal_bam}",
        "--ref_fn", f"{ref_dir}/{ref_fasta_name}",
        "--threads", str(threads),
        "--platform", platform,
        "--output_dir", output_dir,
        "--min_coverage", " 2",
        "--output_prefix", f"{tumor_name}_{normal_name}"
    ]


    print(" ".join(command))
    try:
        subprocess.run(command, check=True)
        print("Clairs run successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error running Clairs: {e}")


def run_small_variant_detection(sample_list, ann_dict, genome, output_dir, threads):
    """ """

    tumor_samples = []
    normal_samples = []
    for sample in sample_list:
        if sample.origin == "tumor":
            tumor_samples.append(sample)
        else:
            normal_samples.append(sample)
    
    for idx,tumor in enumerate(tumor_samples):
        normal = normal_samples[idx]

        run_clairs(output_dir, tumor, normal, genome, threads, platform="ont_r10_dorado_sup_5khz_ssrs")
        sys.exit()