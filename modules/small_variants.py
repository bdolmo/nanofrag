import os
import sys
import subprocess
from multiprocessing import Pool

def generate_contigs(ref_fasta):
    """
        Generate contig information from reference .fai file.
    """
    fai_file = ref_fasta + ".fai"
    contig_file = "contigs.txt"
    
    with open(contig_file, 'w') as out_f:
        with open(fai_file, 'r') as fai_f:
            for line in fai_f:
                contig, length = line.split()[:2]
                out_f.write(f"##contig=<ID={contig},length={length}>\n")
    
    return contig_file


def generate_vcf_header(contig_file, output_dir):
    """
        Generate a complete VCF header programmatically.
    """
    vcf_header = os.path.join(output_dir, "vcf_header.txt")
    
    with open(vcf_header, 'w') as vcf_h:
        # Write contigs
        with open(contig_file, 'r') as contigs_f:
            for line in contigs_f:
                vcf_h.write(line)
        
        # Write other standard headers
        vcf_h.write("##fileformat=VCFv4.2\n")
        vcf_h.write("##INFO=<ID=MBQ,Number=1,Type=Float,Description=\"Mean base quality\">\n")
        vcf_h.write("##INFO=<ID=MMQ,Number=1,Type=Float,Description=\"Mean mapping quality\">\n")
        vcf_h.write("##INFO=<ID=AB,Number=1,Type=Float,Description=\"Allele balance\">\n")
        vcf_h.write("##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Allele depth forward strand\">\n")
        vcf_h.write("##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Allele depth reverse strand\">\n")
        vcf_h.write("##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Allele depth\">\n")
        vcf_h.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depth\">\n")
        vcf_h.write("##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allele depth forward strand\">\n")
        vcf_h.write("##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allele depth reverse strand\">\n")
        vcf_h.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n")
        vcf_h.write("##FORMAT=<ID=SP,Number=1,Type=Float,Description=\"Strand bias\">\n")
        vcf_h.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    
    return vcf_header


def run_mpileup_for_chromosome(chromosome, sample_name, tumor, bed, ref_fasta, output_dir, prior=0.1):
    """
    Run bcftools mpileup for a specific chromosome, pipe output to bcftools call with higher prior, and apply filters.
    """
    vcf_output = os.path.join(output_dir, f"variants_{chromosome}.vcf")
    filtered_vcf_output = os.path.join(output_dir, f"filtered_variants_{sample_name}_{chromosome}.vcf")

    # Step 1: Run bcftools mpileup and pipe to bcftools call
    mpileup_cmd = [
        "bcftools", "mpileup", "-r", chromosome, "-f", ref_fasta, "-X", "ont",
        "-a", "FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR",
        "-T", bed, tumor
    ]

    # Set bcftools call with `-P` and options
    call_cmd = [
        "bcftools", "call", "-m", "-v", "-P", str(prior), "-Ov"
    ]

    with open(vcf_output, 'w') as vcf_f:
        # Use Popen to connect mpileup and call commands via a pipe
        mpileup_proc = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE)
        call_proc = subprocess.Popen(call_cmd, stdin=mpileup_proc.stdout, stdout=vcf_f)
        
        # Ensure mpileup_proc's stdout closes properly before call_proc
        mpileup_proc.stdout.close()
        call_proc.communicate()

    # Step 2: Apply filters to the called VCF
    apply_filters(vcf_output, filtered_vcf_output)

    return filtered_vcf_output


def concatenate_compress_and_sort_vcfs(output_dir, output_vcf_gz):
    """
    Concatenate multiple VCF files into a single gzipped and sorted VCF file.
    
    Args:
        output_dir (str): Directory containing the VCF files to concatenate.
        output_filename (str): Name of the final gzipped and sorted VCF file.
    """
    # Get all filtered VCF files in output_dir
    vcf_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.startswith("filtered_variants_") and f.endswith(".vcf")]

    # Define the path for the output gzipped VCF file
    
    # Temporary file for concatenated but unsorted VCF
    temp_concat_vcf = os.path.join(output_dir, "merged_variants.vcf.gz")
    
    # Concatenate VCF files and compress
    concat_cmd = ["bcftools", "concat"] + vcf_files
    bgzip_cmd = ["bgzip", "-c"]

    with open(temp_concat_vcf, 'wb') as temp_out:
        concat_proc = subprocess.Popen(concat_cmd, stdout=subprocess.PIPE)
        bgzip_proc = subprocess.Popen(bgzip_cmd, stdin=concat_proc.stdout, stdout=temp_out)
        
        concat_proc.stdout.close()
        bgzip_proc.communicate()

    # Sort the gzipped VCF file
    sort_cmd = ["bcftools", "sort", "-Oz", "-o", output_vcf_gz, temp_concat_vcf]
    subprocess.run(sort_cmd)

    index_cmd = ["tabix", "-p", "vcf", output_vcf_gz]
    subprocess.run(index_cmd)

    os.remove(temp_concat_vcf)

    # print(f"  {sorted_vcf_gz}")
    return output_vcf_gz


def run_bcftools(output_dir, sample_name, tumor, bed, ref_fasta, threads):
    """
    Runs bcftools mpileup in parallel for each chromosome and applies variant filters, then concatenates results.
    
    Args:
        output_dir (str): Directory to store the output.
        tumor (str): Path to BAM file of the tumor sample.
        bed (str): Path to BED file for regions of interest.
        ref_fasta (str): Path to the reference fasta file.
        threads (int): Number of parallel jobs.
    """
    
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Step 1: List of chromosomes to process
    chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 
                   'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    
    msg = f" INFO: SNV detection (BCFtools) for {sample_name}"
    print(msg)
    
    # Step 2: Parallelize bcftools mpileup for each chromosome and apply filters
    pool = Pool(processes=threads)
    pool.starmap(run_mpileup_for_chromosome, 
                 [(chrom, sample_name, tumor, bed, ref_fasta, output_dir) for chrom in chromosomes])
    pool.close()
    pool.join()

    # Step 3: Concatenate and compress the filtered VCFs
    msg = " INFO: Bgzip compression, sorting and indexing"
    print(msg)

    snv_vcf = os.path.join(output_dir, f"{sample_name}.snv.vcf.gz")
    merged_vcf_gz = concatenate_compress_and_sort_vcfs(output_dir, snv_vcf)


# def run_mpileup_for_chromosome(chromosome, tumor, bed, ref_fasta, output_dir, prior=0.1):
#     """
#     Run bcftools mpileup for a specific chromosome, pipe output to bcftools call with higher prior, and apply filters.
#     """
#     vcf_output = os.path.join(output_dir, f"variants_{chromosome}.vcf")
#     filtered_vcf_output = os.path.join(output_dir, f"filtered_variants_{chromosome}.vcf")

#     # Step 1: Run bcftools mpileup and pipe to bcftools call
#     mpileup_cmd = [
#         "bcftools", "mpileup", "-r", chromosome, "-f", ref_fasta, "-X", "ont",
#         "-a", "FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR",
#         "-T", bed, tumor
#     ]

#     # Set bcftools call with `-P` and options
#     call_cmd = [
#         "bcftools", "call", "-m", "-v", "-P", str(prior), "-Ov"
#     ]

#     # with open(vcf_output, 'w') as vcf_f:
#     #     # Use Popen to connect mpileup and call commands via a pipe
#     #     mpileup_proc = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE)
#     #     call_proc = subprocess.Popen(call_cmd, stdin=mpileup_proc.stdout, stdout=vcf_f)
        
#     #     # Ensure mpileup_proc's stdout closes properly before call_proc
#     #     mpileup_proc.stdout.close()
#     #     call_proc.communicate()

#     # Step 2: Apply filters to the called VCF
#     apply_filters(vcf_output, filtered_vcf_output)

#     return filtered_vcf_output



def apply_filters(vcf_input, vcf_output):
    """
    Applies filters to the VCF file based on criteria like depth, allele frequency, etc.
    """
    with open(vcf_input, 'r') as vcf_in, open(vcf_output, 'w') as vcf_out:
        for line in vcf_in:
            if line.startswith("#"):
                # Copy VCF headers directly
                vcf_out.write(line)
            else:
                # Parse VCF line and apply filters
                fields = line.strip().split("\t")
                chrom, pos, id_, ref, alt, qual, filter_, info, format_, sample = fields


                
                # Extract necessary fields from the INFO and FORMAT columns
                info_dict = dict(item.split("=") for item in info.split(";") if "=" in item)
                format_fields = format_.split(":")
                sample_fields = sample.split(":")
                format_dict = dict(zip(format_fields, sample_fields))
                
                # Example: Apply filters on depth (DP), allele frequency (AF), strand bias (SP)
                dp = int(format_dict.get("DP", 0))
                af = float(info_dict.get("AF", 0))
                sp = float(format_dict.get("SP", 0))
                ad = format_dict.get("AD", "")

                mapq = float(info_dict.get("MQ", 0))
                baseq = float(info_dict.get("MBP", 0))


                if ad:
                    tmp_ad = ad.split(",")
                    ref_reads = 0
                    alt_reads = 0
                    if len(tmp_ad)>1:
                        ref_reads = int(tmp_ad[0])
                        alt_reads = int(tmp_ad[1])
                    # print(alt_reads, mapq, baseq)

                
                    # Set thresholds for filtering
                    min_dp = 2  # Minimum depth
                    max_sp = 60  # Maximum strand bias
                    min_af = 0.05  # Minimum allele frequency
                    min_alt = 2
                        
                    # Apply the filters
                    if dp >= min_dp and alt_reads >= min_alt and mapq > 20:
                        # print(line)
                        vcf_out.write(line)  # Keep variant if it passes the filters


def run_deepsomatic(output_dir, tumor, ref_fasta, threads, model_type="ONT_TUMOR_ONLY ", regions=None, use_default_pon_filtering=True):
    """
    Function to run DeepSomatic using Docker with specified tumor, normal, reference files, and additional parameters.
    
    Args:
        output_dir (str): The output directory where results will be saved.
        tumor (str): Path to the tumor BAM file.
        normal (str): Path to the normal BAM file.
        ref_fasta (str): Path to the reference FASTA file.
        threads (int): Number of threads to use.
        model_type (str): Type of model to run (e.g., WGS, WES, PACBIO, ONT).
        regions (str, optional): Region of the genome to restrict analysis (e.g., "chr1"). Default is None (whole genome).
        use_default_pon_filtering (bool): Whether to use default PON filtering for tumor-only variant calling.
    """

    # Define directories and file names
    tumor_dir = os.path.dirname(tumor.bam)
    # normal_dir = os.path.dirname(normal.bam)
    ref_dir = os.path.dirname(ref_fasta)

    tumor_bam = os.path.basename(tumor.bam)
    # normal_bam = os.path.basename(normal.bam)
    ref_fasta_name = os.path.basename(ref_fasta)
    tumor_name = "tumor"
    # normal_name = "normal"

    # Define the command for running DeepSomatic with Docker
    # command = [
    #     "sudo", "docker", "run", "-v", f"{tumor_dir}:{tumor_dir}",
    #     "-v", f"{normal_dir}:{normal_dir}",
    #     "-v", f"{output_dir}:{output_dir}",
    #     "google/deepsomatic:1.7.0",
    #     "run_deepsomatic",
    #     f"--model_type={model_type}",
    #     f"--ref={ref_dir}/{ref_fasta_name}",
    #     f"--reads_normal={normal_dir}/{normal_bam}",
    #     f"--reads_tumor={tumor_dir}/{tumor_bam}",
    #     f"--output_vcf={output_dir}/{tumor_name}_{normal_name}.vcf.gz",
    #     f"--output_gvcf={output_dir}/{tumor_name}_{normal_name}.g.vcf.gz",
    #     f"--sample_name_tumor={tumor_name}",
    #     f"--sample_name_normal={normal_name}",
    #     f"--num_shards={threads}",
    #     f"--logging_dir={output_dir}/logs",
    #     f"--intermediate_results_dir={output_dir}/intermediate_results_dir",
    #     f"--use_default_pon_filtering={'true' if use_default_pon_filtering else 'false'}",
    #     f"--dry_run=false"
    # ]

    command = [
        "docker", "run", "-it",
        "-v", f"{tumor_dir}:{tumor_dir}",
        "-v", f"{output_dir}:{output_dir}",
        "-v", f"{ref_dir}:{ref_dir}",
        "google/deepsomatic:1.7.0",
        "run_deepsomatic",
        f"--model_type={model_type}",
        f"--ref={ref_dir}/{ref_fasta_name}",
        f"--reads_tumor={tumor_dir}/{tumor_bam}",
        f"--output_vcf={output_dir}/{tumor_name}.vcf.gz",
        f"--sample_name_tumor={tumor_name}",
        f"--num_shards={threads}",
        f"--logging_dir={output_dir}/logs",
        f"--intermediate_results_dir={output_dir}/intermediate_results_dir",
        f"--use_default_pon_filtering={'true' if use_default_pon_filtering else 'false'}"
    ]





    # Add regions if provided
    if regions:
        command.append(f"--regions={regions}")

    # Print and run the command
    print("Running DeepSomatic with command: " + " ".join(command))
    try:
        subprocess.run(command, check=True)
        print("DeepSomatic ran successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error running DeepSomatic: {e}")


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

    snv_folder = os.path.join(output_dir, "SNV")
    if not os.path.isdir(snv_folder):
        os.mkdir(snv_folder)


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
        # run_deepsomatic(output_dir, tumor, genome, threads, model_type="ONT_TUMOR_ONLY")

        run_bcftools(snv_folder, tumor.name, tumor.bam, ann_dict["somatic_variants"], genome, threads)


        sys.exit()