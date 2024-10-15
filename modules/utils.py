import os
import sys


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

def set_binaries_configuration(main_dir):
    """ """
    bin_dict = {
        "wigToBigWig" : os.path.join(main_dir, "bin", "wigToBigWig"),
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
        "chromosomes" : os.path.join(main_dir, "annotations", "hg38.chromosomes.txt")
    }
    return ann_dict



def get_bams_from_list(input_list_txt):
    """ """

    bams = []

    if not os.path.isfile(input_list_txt):
        msg = f" INFO: Missing input list txt file {input_list_txt}"
        print(msg)
        sys.exit()

    with open(input_list_txt) as f:
        for line in f:
            line = line.rstrip("\n")
            bams.append(line)

    return bams

