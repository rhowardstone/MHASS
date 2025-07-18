# MHASS

**Microbiome HiFi Amplicon Sequencing Simulator**

A complete pipeline for simulating PacBio HiFi amplicon sequencing data for microbiome studies.

## Installation

```bash
git clone https://github.com/rhowardstone/MHASS.git
cd MHASS
bash install_dependencies.sh
```

After installation:

```bash
conda activate mhass
```

## Usage

### Basic Command

```bash
mhass --amplicon-fasta amplicons.fa --amplicon-genome-labels labels.tsv --output-dir output/
```

### Full Command with Options

```bash
usage: mhass [-h] --amplicon-fasta AMPLICON_FASTA --amplicon-genome-labels AMPLICON_GENOME_LABELS --output-dir OUTPUT_DIR
             [--num-samples NUM_SAMPLES] [--num-reads NUM_READS] [--var-intercept VAR_INTERCEPT] [--var-slope VAR_SLOPE]
             [--genome-distribution GENOME_DISTRIBUTION] [--barcode-file BARCODE_FILE] [--subread-accuracy SUBREAD_ACCURACY]
             [--np-distribution-type {empirical,lognormal}] [--lognormal-mu LOGNORMAL_MU] [--lognormal-sigma LOGNORMAL_SIGMA]
             [--np-min NP_MIN] [--np-max NP_MAX] [--np-distribution NP_DISTRIBUTION] [--threads THREADS]

options:
  -h, --help            show this help message and exit
  --amplicon-fasta AMPLICON_FASTA
                        Input FASTA file of amplicons (default: None)
  --amplicon-genome-labels AMPLICON_GENOME_LABELS
                        TSV file mapping amplicons to genomes (asvid<TAB>genomeid) (default: None)
  --output-dir OUTPUT_DIR
                        Output directory for all files (default: None)
  --num-samples NUM_SAMPLES
                        Number of samples to simulate (default: 10)
  --num-reads NUM_READS
                        Number of reads per sample (default: 10000)
  --var-intercept VAR_INTERCEPT
                        Intercept for ASV variability model (controls baseline variation between samples) (default: 1.47565981333483)
  --var-slope VAR_SLOPE
                        Slope for ASV variability model (how variation changes with abundance) (default: -0.909890963463704)
  --genome-distribution GENOME_DISTRIBUTION
                        Distribution for genome abundances (uniform, lognormal, powerlaw, or empirical:<file>) (default: uniform)
  --barcode-file BARCODE_FILE
                        TSV file with barcodes (id, forward, reverse) (default: None)
  --subread-accuracy SUBREAD_ACCURACY
                        Mean subread accuracy used in PBSIM (default: 0.65) (default: 0.65)
  --np-distribution-type {empirical,lognormal}
                        Type of np distribution: empirical (from file) or lognormal (default: empirical)
  --lognormal-mu LOGNORMAL_MU
                        mu parameter for lognormal distribution of Number of Passes (default from empirical fit) (default: 3.88)
  --lognormal-sigma LOGNORMAL_SIGMA
                        sigma parameter for lognormal distribution of Num Passes (default from empirical fit) (default: 1.22)
  --np-min NP_MIN       Minimum np value when using lognormal distribution (default: 2)
  --np-max NP_MAX       Maximum np value when using lognormal distribution (default: 59)
  --np-distribution NP_DISTRIBUTION
                        TSV file with empirical num-passes distribution (default: None)
  --threads THREADS     Number of threads to use for parallel processing (default: 190)
```

## Parameters

### Required Parameters

* `--amplicon-fasta`: Input FASTA file of amplicon reference sequences
* `--amplicon-genome-labels`: TSV file mapping amplicons to genome IDs (`asvid` and `genomeid` columns)
* `--output-dir`: Directory where output files will be written

### Optional Parameters

* `--num-samples`: Number of simulated samples (default: 10)
* `--num-reads`: Number of reads per sample (default: 10,000)
* `--dispersion`: Dispersion parameter for count simulation (default: 0.1)
* `--genome-distribution`: `uniform`, `lognormal`, `powerlaw`, or `empirical:<file.tsv>`
* `--barcode-file`: TSV file with barcodes (default: bundled file)
* `--np-distribution-type`: `empirical` or `gamma` (default: `empirical`)
* `--np-distribution`: TSV file for empirical number-of-passes distribution
* `--lognormal-mu`: Lognormal mean abundance parameter if using `lognormal` np distribution (default: 3.88)
* `--lognormal-sigma`: Lognormal np variability parameter (default: 1.22)
* `--np-min`: Minimum np value (default: 2)
* `--np-max`: Maximum np value (default: 59)
* `--subread-accuracy`: Subread accuracy for PBSIM (default: 0.65)
* `--threads`: Number of threads for parallel PBSIM + CCS execution (default: all CPUs)

## Input File Formats

### 1. Amplicon FASTA

```
>ASV001
ATCGATCGATCGATCG...
```

### 2. Amplicon-to-Genome Label TSV

```
asvid	genomeid
ASV001	Genome1
ASV002	Genome2
```

### 3. Barcode File (optional)

```
BarcodeID	ForwardBarcode	ReverseBarcode
BC01	CACATATCAGAGTGCG	TGGCGTGCATGATTCGA
```

### 4. Num-Passes Distribution (optional)

```
num_passes	count
2	19
3	23159
...
```

### 5. Empirical Genome Abundance (if used)

```
genomeid	proportion
Genome1	0.45
Genome2	0.30
```

## Output Files

* `counts.tsv`: Simulated ASV counts per sample
* `counts_meta.tsv`: Sample library sizes
* `sample_barcode_map.tsv`: Mapping of samples to barcodes
* `sequence_file_mapping.tsv`: Detailed template mappings
* `combined_reads.fastq`: Final output FASTQ file

## Workflow

1. **Count Simulation**: Uses `metaSPARSim` to simulate ASV abundance across samples based on genome distribution.
2. **Template Creation**: Templates are created per ASV×?Sample with barcodes and sampled np values.
3. **Read Simulation**: PBSIM3 simulates PacBio subreads per template with user-defined `--subread-accuracy`.
4. **CCS Processing**: Subreads are collapsed into CCS reads using `ccs`.
5. **Read Relabeling**: All CCS reads are relabeled and merged into `combined_reads.fastq`.
6. **Cleanup**: Intermediate files are cleaned up, leaving only final outputs.

## Example

```bash
mhass --amplicon-fasta test/amplicons.fa \
      --amplicon-genome-labels test/labels.tsv \
      --output-dir output/ \
      --num-samples 5 \
      --num-reads 1000 \
      --threads 4 \
      --subread-accuracy 0.9
```

## Citation

TBD ?? please cite appropriately once the tool is published.
