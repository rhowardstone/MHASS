#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import shutil
import tempfile
from pathlib import Path
import multiprocessing

# Import local modules (when running as a package)
try:
    from mhass.template_creator import create_per_sequence_templates
    from mhass.pbsim_runner import run_pbsim
    from mhass.fastq_combiner import combine_fastqs
except ImportError:
    # For development/testing
    from template_creator import create_per_sequence_templates
    from pbsim_runner import run_pbsim
    from fastq_combiner import combine_fastqs

def get_package_path():
    """Get the path to the package resources."""
    return Path(__file__).resolve().parent

def get_resource_path(resource_name):
    """Get the path to a resource file bundled with the package."""
    return get_package_path() / "resources" / resource_name

def run_r_count_script(fasta_path, genome_map, num_samples, num_reads, dispersion, 
                      genome_distribution, output_counts, output_meta):
    """Run the R script to generate count matrix."""
    r_script_path = get_package_path() / "get_counts.R"
    
    cmd = [
        "Rscript", str(r_script_path),
        "-f", str(fasta_path),
        "-n", str(num_samples),
        "-r", str(num_reads),
        "-d", str(dispersion),
        "-G", str(genome_map),
        "--genome-distribution", genome_distribution,
        "-o", str(output_counts),
        "-m", str(output_meta)
    ]
    
    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=True)
    return result.returncode == 0

def run_simulation(args):
    """Run the entire simulation pipeline."""
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set up file paths
    counts_file = output_dir / "counts.tsv"
    meta_file = output_dir / "counts_meta.tsv"
    templates_dir = output_dir / "templates"
    barcode_mapping = output_dir / "sample_barcode_map.tsv"
    sim_reads_dir = output_dir / "sim_reads"
    combined_reads = output_dir / "combined_reads.fastq"
    
    # Get default resources
    default_barcode_file = get_resource_path("barcodes.txt")
    default_np_dist = get_resource_path("np_distribution.tsv")
    
    # Use defaults or provided values
    barcode_file = args.barcode_file or default_barcode_file
    np_distribution = args.np_distribution or default_np_dist
    
    # Step 1: Simulate count matrix
    print("\n==> STEP 1: Simulating count matrix <==\n")
    success = run_r_count_script(
        args.amplicon_fasta, 
        args.amplicon_genome_labels,
        args.num_samples,
        args.num_reads,
        args.dispersion,
        args.genome_distribution,
        counts_file,
        meta_file
    )
    if not success:
        print("Error: Failed to generate count matrix")
        return False
    
    # Step 2: Create templates
    print("\n==> STEP 2: Creating templates <==\n")
    templates_dir.mkdir(parents=True, exist_ok=True)
    
    # Prepare np distribution parameters
    np_params = {
        'distribution_type': args.np_distribution_type,
        'empirical_file': np_distribution if args.np_distribution_type == 'empirical' else None,
        'lognormal_mu': args.lognormal_mu,
        'lognormal_sigma': args.lognormal_sigma,
        'np_min': args.np_min,
        'np_max': args.np_max
    }
    
    create_per_sequence_templates(
        args.amplicon_fasta,
        counts_file,
        templates_dir,
        barcode_file,
        barcode_mapping,
        np_params  # Pass the np_params dict instead of just the file
    )
    
    # Step 3: Simulate and process reads
    print("\n==> STEP 3: Simulating and processing reads in parallel <==\n")
    # Get the path to PBSIM and QSHMM from package resources
    pbsim_path = get_resource_path("pbsim3/src/pbsim")
    qshmm_path = get_resource_path("pbsim3/data/QSHMM-RSII.model")
    
    # If running in development mode, try to find PBSIM in system PATH
    if not pbsim_path.exists():
        pbsim_path = shutil.which("pbsim") or "pbsim3/src/pbsim"
        qshmm_path = Path(os.environ.get("PBSIM_QSHMM", "pbsim3/data/QSHMM-RSII.model"))
    
    sim_reads_dir.mkdir(parents=True, exist_ok=True)
    run_pbsim(
        templates_dir,
        sim_reads_dir,
        pbsim_path,
        qshmm_path,
        args.threads,
        args.subread_accuracy
    )

    # Step 4: Combine and relabel CCS reads
    print("\n==> STEP 4: Combining and relabeling CCS reads <==\n")
    combine_fastqs(sim_reads_dir, combined_reads)

    print("\n==> STEP 5: Cleaning up intermediate files <==\n")
    # Move sequence mapping file up one directory
    sequence_mapping_src = templates_dir / "sequence_file_mapping.tsv"
    sequence_mapping_dest = output_dir / "sequence_file_mapping.tsv"
    
    if sequence_mapping_src.exists():
        shutil.move(str(sequence_mapping_src), str(sequence_mapping_dest))
        print(f"Moved sequence mapping file to: {sequence_mapping_dest}")
    else:
        print("Warning: sequence_file_mapping.tsv not found in templates directory")
    
    # Remove templates directory
    if templates_dir.exists():
        shutil.rmtree(templates_dir)
        print(f"Removed templates directory: {templates_dir}")
    
    # Remove sim_reads directory
    if sim_reads_dir.exists():
        shutil.rmtree(sim_reads_dir)
        print(f"Removed sim_reads directory: {sim_reads_dir}")
    
    print("\n==> FINISHED <==\n")
    print(f"Final output: {combined_reads}")
    return True

def main():
    parser = argparse.ArgumentParser(
        description="Microbiome Sequencing Simulator - A complete pipeline for simulating PacBio amplicon data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("--amplicon-fasta", required=True, 
                        help="Input FASTA file of amplicons")
    parser.add_argument("--amplicon-genome-labels", required=True,
                        help="TSV file mapping amplicons to genomes (asvid<TAB>genomeid)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for all files")
    
    # Optional arguments
    parser.add_argument("--num-samples", type=int, default=10,
                        help="Number of samples to simulate")
    parser.add_argument("--num-reads", type=int, default=10000,
                        help="Number of reads per sample")
    parser.add_argument("--dispersion", type=float, default=0.1,
                        help="Dispersion parameter for count simulation")
    parser.add_argument("--genome-distribution", default="uniform",
                        help="Distribution for genome abundances (uniform, lognormal, powerlaw, or empirical:<file>)")
    parser.add_argument("--barcode-file", default=None,
                        help="TSV file with barcodes (id, forward, reverse)")
    parser.add_argument("--subread-accuracy", type=float, default=0.85,
                        help="Mean subread accuracy used in PBSIM (default: 0.85)")
    parser.add_argument("--np-distribution-type", default="empirical",
                    choices=["empirical", "lognormal"],
                    help="Type of np distribution: empirical (from file) or lognormal")
    parser.add_argument("--lognormal-mu", type=float, default=3.88,
                    help="μ parameter for lognormal distribution of NP (default from empirical fit)")
    parser.add_argument("--lognormal-sigma", type=float, default=1.22,
                    help="σ parameter for lognormal distribution of NP (default from empirical fit)")
    parser.add_argument("--np-min", type=int, default=2,
                        help="Minimum np value when using lognormal distribution")
    parser.add_argument("--np-max", type=int, default=50,
                        help="Maximum np value when using lognormal distribution")
    parser.add_argument("--np-distribution", default=None,
                        help="TSV file with empirical num-passes distribution")
    parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count(),
                        help="Number of threads to use for parallel processing")
    args = parser.parse_args()
    
    success = run_simulation(args)
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
