import os
import sys
import gzip
import glob

tissues_dict = {
    # Blood cells - soft yellow shades (control-like)
    "erythrocyte_progenitor": "#FFD77F",  # slightly deeper yellow-orange
    "T-cell": "#FFE29F",  # soft yellow
    "NK-cell": "#ffdf7a",  # slightly lighter yellow
    "monocyte": "#ffecb3",  # very light yellow for higher frequency
    "granulocyte": "#FFF9E6",  # lightest yellow for highest frequency
    "B-cell": "#FFE6A0",  # soft yellow-orange

    # Epithelial cells - softened green shades
    "ovary_epithelial": "#c9ede3",
    "thyroid_epithelial": "#a8dfcf",
    "fallopian_epithelial": "#89d2ba",
    "lung_epithelial": "#6bc4a6",
    "prostate_epithelial": "#5bb599",
    "bladder_epithelial": "#a3e3c0",
    "head_neck_epithelial": "#e4f7f2",
    "epiderminal_keratinocyte": "#308570",
    "colon_epithelial": "#267663",
    "small_int_epithelial": "#3d967e",
    "gastric_epithelial": "#4ba68c",#4ba68c

    # Muscle and connective tissue - orange shades
    "heart_fibroblasts": "#ffab73",
    "smooth_muscle": "#ff8c42",
    "skeletal_muscle": "#ff6f00",
    "heart_cardiomyocyte": "#e65c00",

    # Fibroblasts and connective tissue - blue shades
    "endothelial_cell": "#b3cde0",
    "dermal_fibroblast": "#6497b1",
    "colon_fibroblasts": "#005b96",
    "heart_fibroblasts": "#03396c",
    "bone_osteoblast": "#011f4b",

    # Pancreatic cells - purple shades
    "pancreas_duct": "#d7bde2",
    "pancreas_acinar": "#c39bd3",
    "pancreas_delta": "#af7ac5",
    "pancreas_beta": "#9b59b6",
    "pancreas_alpha": "#76448a",

    # Brain and nervous system - brown shades
    "neuron": "#d7ccc8",
    "oligodendrocyte": "#bcaaa4",

    # Adipose tissue and other - gray and earth tones
    "adipocyte": "#b0bec5",
    "hepatocyte": "#8d6e63",
    "kidney_epithelial": "#a1887f",
    "gallbladder": "#6d4c41",

    # Mammary and reproductive tissues - pink shades
    "breast_luminal": "#f8bbd0",
    "breast_basal": "#f48fb1"
}



def bed_to_list(input_bed, gzipped=False):
    """ """

    regions_list = []

    if gzipped:
        with gzip.open(input_bed,'rt') as f:
            for line in f:
                line = line.rstrip("\n")
                regions_list.append(line)
        f.close()
    else:
        with open(input_bed) as f:
            for line in f:
                line = line.rstrip("\n")
                regions_list.append(line)
        f.close()
    return regions_list



class Sample:
    def __init__(self, name):
        self._name = name

    def add(self, key, value):
        """
        Add a new attribute dinamically
        """

        if not key:
            raise ValueError("Expected key and val to be provided.")
        setattr(self, key, value)

    @property
    def name(self):
        return self._name
    

def get_chromosome_sizes(chrom_sizes):
    """ """
    chrom_sizes_dict = {}
    with open(chrom_sizes) as f:
        for line in f:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            
            chrom = tmp[0]
            end = int(tmp[1])
            chrom_sizes_dict[chrom] = end
    return chrom_sizes_dict


def create_windows_bed_wig(chrom_sizes_file, window_size, output_wig_file, output_bed_file):
    # Open the chromosome sizes file, WIG file, and BED file
    with open(chrom_sizes_file, 'r') as f, open(output_wig_file, 'w') as wig, open(output_bed_file, 'w') as bed:
        # Iterate over each chromosome in the chromosome sizes file
        for line in f:
            refName, refLength = line.strip().split()
            if "M" in refName:
                continue
            refLength = int(refLength)
            
            # Adjust window size if it's greater than the chromosome length
            window = min(window_size, refLength)
            
            # Write WIG header for each chromosome
            wig.write(f"fixedStep chrom={refName} start=1 step={window} span={window}\n")
            
            # Initialize start and end positions for the first window
            start = 0
            while start < refLength:
                # Calculate end of the current window
                end = min(start + window, refLength)
                
                # Write BED entry (0-based) and WIG value (use 0 as placeholder for now)
                bed.write(f"{refName}\t{start}\t{end}\n")
                wig.write("0\n")  # Replace with actual data if available
                
                # Move to the next window
                start += window


def set_binaries_configuration(main_dir):
    """ """
    bin_dict = {
        "modkit": os.path.join(main_dir, "bin", "modkit"),
        "wigToBigWig" : os.path.join(main_dir, "bin", "wigToBigWig"),
        "cfdna_counter": os.path.join(main_dir, "bin", "cfDNA_counts", "src", "cfdna_counter"),
    }
    return bin_dict


def set_sample_configuration(tumor_bams, normal_bams):
    """ """
    sample_list = []

    for bam in tumor_bams:
        sample_name = os.path.basename(bam).replace(".bam", "")
        sample = Sample(sample_name)
        sample.add("origin", "tumor")
        sample.add("bam", bam)
        sample_list.append(sample)

    for bam in normal_bams:
        sample_name = os.path.basename(bam).replace(".bam", "")
        sample = Sample(sample_name)
        sample.add("origin", "normal")
        sample.add("bam", bam)
        sample_list.append(sample)
    return sample_list



def set_annotation_resources(main_dir):
    """ """
    ann_dict = {
        "blacklist" : os.path.join(main_dir, "annotations", "consensusBlacklist.hg38.bed"),
        "chromosomes" : os.path.join(main_dir, "annotations", "hg38.chromosomes.txt"),
        "nucleosomes": os.path.join(main_dir, "annotations", "GSE71378_CH01.hg38.reduced.bed.gz"),
        "tss": os.path.join(main_dir, "annotations", "refTSS_v4.1_human_coordinate.hg38.bed"),
        "somatic_variants": os.path.join(main_dir, "annotations", "somatic.variants.grhc38.bed"),
        "cpg_islands": os.path.join(main_dir, "annotations", "cpg_islands_ucsc_cleaned.bed"),
        "nanomix_atlas": os.path.join(main_dir, "annotations", "39Bisulfite.tsv"),
        # "celfie": os.path.join(main_dir, "annotations", "celfie.bed")
    }
    return ann_dict



def get_input_bams(input_bams):
    """ """

    bams = []


    if os.path.isfile(input_bams):
        with open(input_bams) as f:
            for line in f:
                line = line.rstrip("\n")
                bams.append(line)
    if os.path.isdir(input_bams):
        bams = glob.glob(input_bams+"/*.bam")

    if not bams:
        msg = f" INFO: Missing input bam files at {input_bams}"
        print(msg)
        sys.exit()

    return bams

