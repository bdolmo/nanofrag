import os
import sys
import argparse
from modules.utils import get_input_bams, set_annotation_resources, set_binaries_configuration, set_sample_configuration
from modules.fragmentomics import run_fragmentomic_analysis
from modules.nucleosome import run_wps_analysis, run_tss_analysis
from modules.copy_number import run_cn_workflow
from modules.small_variants import run_small_variant_detection
from modules.methylation import run_methylation_analysis

main_dir = os.path.dirname(os.path.realpath(__file__))

def get_args():

    parser = argparse.ArgumentParser(prog="NanoFrag", 
        description="Analaysis of DNA fragmentation from liquid biopsy using nanopore")

    parser.add_argument("--tumor_list", dest="tumor_list", type=str, required=True, 
        help="Analaysis of DNA fragmentation from liquid biopsy using nanopore")
    parser.add_argument("--normal_list", dest="normal_list", type=str, required=True)
    parser.add_argument("--reference", dest="reference", type=str, required=True)
    parser.add_argument("--output_dir", dest="output_dir", type=str, required=True)
    parser.add_argument("--docker_output",  dest="docker_output", type=str, required=False)

    parser.add_argument("--threads", dest="threads", type=int, required=True)
    parser.add_argument("--skip_fragmentation", dest="skip_fragmentation", 
        help="Skip fragment size histogram calculation", action="store_true")
    parser.add_argument("--skip_wps", dest="skip_wps", 
        help="Skip windowed protection score calculation", action="store_true")
    parser.add_argument("--skip_cn", dest="skip_cn", 
        help="Skip copy number analysis", action="store_true")
    parser.add_argument("--skip_small_variants", dest="skip_small_variants", 
        help="Skip small variant analysis (SNV)", action="store_true")
    parser.add_argument("--skip_methylation", dest="skip_methylation", 
        help="Skip methylation anaylsis", action="store_true")

    args = parser.parse_args()
    return args



if __name__ == "__main__":
    args = get_args()

    tumor_list = args.tumor_list
    normal_list = args.normal_list
    output_dir = args.output_dir
    genome = args.reference
    threads = args.threads
    skip_fragmentation = args.skip_fragmentation
    skip_wps = args.skip_wps
    skip_cn = args.skip_cn
    skip_small_variants = args.skip_small_variants
    skip_methylation = args.skip_methylation

    if not hasattr(args, "docker_output"):
        docker_output = output_dir
    else:
        docker_output = args.docker_output


    tumor_bams = get_input_bams(tumor_list)
    normal_bams = get_input_bams(normal_list)

    ann_dict = set_annotation_resources(main_dir)
    sample_list =  set_sample_configuration(tumor_bams, normal_bams)
    bin_dict = set_binaries_configuration(main_dir)

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    if skip_fragmentation == False:
        sample_list = run_fragmentomic_analysis(sample_list, ann_dict, bin_dict, genome, output_dir, threads, window_size=5000000)
    
    if skip_cn == False:
        sample_list = run_cn_workflow(sample_list, docker_output, ann_dict, output_dir)

    if skip_methylation == False:
        sample_list = run_methylation_analysis(sample_list, ann_dict, bin_dict, threads, genome, docker_output, output_dir)

    sys.exit()


    if skip_small_variants == False:
        run_small_variant_detection(sample_list, ann_dict, genome, output_dir, threads)

    # if skip_wps == False:
    #     # windowed_protection_scores(sample_list, ann_dict, output_dir)
    #     # region = "chr12:34267065-34294065"
    #     # region = ""

    #     run_tss_analysis(sample_list, ann_dict, output_dir)
    #     # run_nucleosome_profiling(sample_list, ann_dict, output_dir)
    #     # run_wps_analysis(sample_list, ann_dict, bin_dict, output_dir)



    sys.exit()


    

