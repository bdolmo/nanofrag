import os
import sys



def get_annotation_resources(main_dir):
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

