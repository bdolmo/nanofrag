#!/bin/env perl 

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Basename;
use Parallel::ForkManager;

# Variables for options
my ($input_dir, $reference, $threads, $output_dir);

# Default values
$threads = 4; # Default number of threads

# Parse command-line arguments
GetOptions(
    "input_dir=s"  => \$input_dir,   # Input directory with unaligned BAM file
    "reference=s"  => \$reference,   # Reference genome for minimap2
    "threads=i"    => \$threads,     # Number of threads for minimap2 (default = 4)
    "output_dir=s" => \$output_dir,  # Output directory for the SAM file
) or die("Error in command line arguments\n");


sub help {
   print "Usage: $0 --input-dir <input_directory> --reference <reference_genome> --output-dir <output_directory> [--threads <num_threads>]\n";
   exit;
}

# Check required arguments
if (!$input_dir) {
    help();
}
if (!$reference) {
    help();
}
if (!$output_dir) {
    help();
}

mkdir $output_dir;

# Get samtools and minimap2 paths
my $samtools = `which samtools`;
chomp $samtools;
if (!$samtools) {
    print "samtools not found in PATH\n";
    exit;
}

my $minimap2 = `which minimap2`;
chomp $minimap2;
if (!$minimap2) {
    print "minimap2 not found in PATH\n";
    exit;
}


my @unaligned_bams = glob("$input_dir/*.bam");
if (!@unaligned_bams) {
    print " ERROR: no bam files were found at $input_dir\n";
    exit;
}

my $sample = basename($unaligned_bams[0]);
my @tmpSample = split("_", $sample);
my $sampleName = $tmpSample[0];

# my $idx = 0;
# foreach my $bam(@unaligned_bams) {
#     my $cmd = "$samtools fastq -@ 4 -T MM,ML $bam  | gzip > $output_dir/$sampleName.$idx.fastq.gz";
#     print $cmd . "\n";
#     system $cmd;
#     $idx++;
# }



# Number of parallel processes you want to run (you can adjust this)
my $max_processes = $threads; 

# Initialize ForkManager with the desired number of processes
my $pm = Parallel::ForkManager->new($max_processes);

my $idx = 0;
foreach my $bam (@unaligned_bams) {
    $idx++;

    # Start a new forked process
    $pm->start and next;  # Fork a process and move to the next iteration
    
    my $cmd = "$samtools fastq -@ 4 -T MM,ML $bam  | gzip > $output_dir/$sampleName.$idx.fastq.gz";
    print $cmd . "\n";
    system $cmd;
    $pm->finish;  # End the forked process
}

# Wait for all child processes to complete
$pm->wait_all_children;


my $cmd = "cat $output_dir/*.fastq.gz > $output_dir/$sampleName.fastq.gz";
if (!-e "$output_dir/$sampleName.fastq.gz") {
    system $cmd;
}

sub align_with_minimap2 {
    my $fastq = shift;
    my $bam = shift;
    my $cmd = "$minimap2 -t $threads -y -ax map-ont $reference $fastq "
}




# # Define the output file paths
# my $fastq_file = "$output_dir/unaligned.fastq";
# my $output_sam = "$output_dir/aligned.sam";

# # Convert BAM to FASTQ with methylation marks (MM, ML)
# print "Converting BAM to FASTQ with methylation marks...\n";
# my $bam_to_fastq_cmd = "$samtools fastq -@ $threads -c -O --output-tags MM,ML $bam_file > $fastq_file";
# system($bam_to_fastq_cmd) == 0 or die "Error converting BAM to FASTQ: $!\n";

# # Map the FASTQ file to the reference genome using minimap2
# print "Mapping FASTQ to the genome using minimap2...\n";
# my $minimap2_cmd = "$minimap2 -t $threads -ax map-ont $reference $fastq_file > $output_sam";
# system($minimap2_cmd) == 0 or die "Error mapping FASTQ to reference: $!\n";

# print "Alignment complete. Output SAM: $output_sam\n";

# # Optional: Uncomment if you want to remove the FASTQ file after mapping
# # unlink $fastq_file;

# exit 0;
