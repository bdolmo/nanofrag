# NanoFrag

Analyze DNA fragmentation and other patterns from liquid biopsy using nanopore sequencing.

---

## Features

- **Fragmentation Analysis**: fragment size histogram and ratios.
- **Copy Number Alterations**: CNA detection using ichorCNA.
- **Methylation Analysis**: Get methylation patterns genome-wide, perform deconvolution into cell/tissues.
- **Small Variant Detection**: detect SNVs (not ready", experimental)

---

## Installation

To clone the repository:

```bash
git clone https://github.com/bdolmo/nanofrag.git
cd nanofrag
```

Install dependencies:

- python 3.10+ is required.
- Install python dependencies with:

```bash
pip install -r requirements.txt
```

- **Docker** is required to run specific tools such as [ichorCNV](https://github.com/broadinstitute/ichorCNA) and [nanomix](https://github.com/Jonbroad15/nanomix)
---

## Docker Images

- **`seqeralabs/ichorcna:latest`**: Used for CNV analysis.
- **`bdolmo/nanomix:1.0.0`**: Used for methylation and fragmentation analysis.

You can pull them using the following command:

```bash
docker pull seqeralabs/ichorcna:latest
docker pull bdolmo/nanomix:1.0.0
```

---

## Usage

```bash
python nanofrag.py --tumor_list <tumor_bam_list> --normal_list <normal_bam_list> \
                   --reference <reference_genome> --output_dir <output_dir> \
                   --threads <num_threads> [optional flags]
```

### Required arguments

- `--tumor_list`: Path to a text file containing tumor BAM file paths (one per line).
- `--normal_list`: Path to a text file containing normal BAM file paths (one per line).
- `--reference`: Reference genome in FASTA format.
- `--output_dir`: Directory for output files.
- `--threads`: number of threads to use.

### Optional flags (if you wish to skip any step)

- `--skip_fragmentation`: Skip fragmentation analysis.
- `--skip_cn`: Skip copy number analysis.
- `--skip_methylation`: Skip methylation analysis.
- `--skip_small_variants`: Skip small variant analysis (not ready)

---

## Example

```bash
python nanofrag.py --tumor_list tumor_samples.txt \
                   --normal_list normal_samples.txt \
                   --reference hg38.fasta \
                   --output_dir results/ \
                   --threads 8
```

---
## Important!

NanoFrag assumes that all bam files have been generated with methylation compatible tags (MM, and ML).


## Output

The pipeline produces output files in the specified `output_dir`:

- **Fragmentation Analysis**:
  - Fragment size histograms and metrics.
- **Copy Number Analysis**:
  - CNV profiles in standard formats (e.g., `.seg` or `.bed`).
- **Methylation Analysis**:
  - Methylation profiles in `.bed`, deconvolution results in `.txt`.
- **Small Variant Detection**: (not ready)
  - VCF files with identified SNVs.
---


