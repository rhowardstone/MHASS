# MHASS
**Microbiome HiFi Amplicon Sequencing Simulator**

A complete pipeline for simulating PacBio HiFi amplicon sequencing data for microbiome studies.

## Installation

```bash
git clone https://github.com/rhowardstone/MHASS.git
cd MHASS
bash install_dependencies.sh
```

This will take some time to install all dependencies. After installation:

```bash
conda activate mhass
```

## Usage

### Basic Command

```bash
mhass --amplicon-fasta amplicons.fa \
      --amplicon-genome-labels labels.tsv \
      --output-dir output/
```

### Full Command with Options

```bash
mhass --amplicon-fasta amplicons.fa \
      --amplicon-genome-labels labels.tsv \
      --output-dir output/ \
      --num-samples 20 \
      --num-reads 5000 \
      --dispersion 0.15 \
      --genome-distribution lognormal \
      --threads 16
```

## Parameters

### Required Parameters
- `--amplicon-fasta`: Input FASTA file containing amplicon sequences
- `--amplicon-genome-labels`: TSV file mapping amplicons to genomes
- `--output-dir`: Directory where all output files will be written

### Optional Parameters
- `--num-samples`: Number of samples to simulate (default: 10)
- `--num-reads`: Number of reads per sample (default: 10,000)
- `--dispersion`: Dispersion parameter for count simulation (default: 0.1)
- `--genome-distribution`: Distribution for genome abundances (default: uniform)
  - Options: `uniform`, `lognormal`, `powerlaw`, or `empirical:<file.tsv>`
- `--barcode-file`: Custom TSV file with barcodes (default: uses built-in barcodes)
- `--np-distribution`: Custom TSV file with num-passes distribution (default: uses built-in distribution)
- `--threads`: Number of threads for parallel processing (default: all available cores)

## Input File Formats

### 1. Amplicon FASTA File
Standard FASTA format with amplicon sequences:
```
>ASV001
ATCGATCGATCGATCG...
>ASV002
GCTAGCTAGCTAGCTA...
```

### 2. Amplicon-Genome Labels File (TSV)
Tab-separated file mapping amplicons to genomes:
```
asvid	genomeid
ASV001	Genome1
ASV002	Genome1
ASV003	Genome2
ASV004	Genome3
```

### 3. Barcode File (TSV) - Optional
Tab-separated file with barcode sequences:
```
BarcodeID	ForwardBarcode	ReverseBarcode
BC01	CACATATCAGAGTGCG	TGGCGTGCATGATTCGA
BC02	ACACACAGACTGTGAG	ACGCACGACATGGACAT
BC03	ACACATCTCGTGAGAG	ACGAGACACTCACATGA
```

### 4. NP Distribution File (TSV) - Optional
Tab-separated file with num-passes distribution:
```
num_passes	count
1	20
2	50
3	100
4	200
5	250
```

### 5. Empirical Abundance File (TSV) - Optional
For `--genome-distribution empirical:<file>`:
```
genomeid	proportion
Genome1	0.45
Genome2	0.30
Genome3	0.25
```

## Output Files

The simulation creates the following files in the output directory:

- `counts.tsv`: Simulated count matrix (ASVs × Samples)
- `counts_meta.tsv`: Metadata for count matrix (sample library sizes)
- `sample_barcode_map.tsv`: Mapping between samples and barcodes
- `sequence_file_mapping.tsv`: Mapping of sequences to template files
- `combined_reads.fastq`: Final combined FASTQ file with simulated reads

## Workflow

MHASS performs the following steps:

1. **Count Simulation**: Uses metaSPARSim to generate realistic count matrices based on genome abundances
2. **Template Creation**: Creates barcoded template sequences for each ASV×Sample combination
3. **Read Simulation**: Uses PBSIM3 to simulate PacBio reads from templates
4. **CCS Processing**: Generates Circular Consensus Sequences using PacBio CCS
5. **Output Generation**: Combines all reads into a single FASTQ file
6. **Cleanup**: Removes intermediate files, keeping only essential outputs

## Example

```bash
# Create test data structure
mkdir -p test_data

# Run simulation with 5 samples, 1000 reads each
mhass --amplicon-fasta test_data/amplicons.fa \
      --amplicon-genome-labels test_data/labels.tsv \
      --output-dir simulation_output/ \
      --num-samples 5 \
      --num-reads 1000 \
      --threads 8
```

## Dependencies

- Python 3.6+
- R 4.0+
- PBSIM3
- PacBio CCS tools
- Required R packages: metaSPARSim, Biostrings, optparse
- Required Python packages: tqdm

All dependencies are automatically installed by the installation script.

## Citation

If you use MHASS in your research, please cite:

```
[Add citation information when available]
```
