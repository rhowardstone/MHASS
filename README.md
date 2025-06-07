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
      --barcode-file my_barcodes.tsv \
      --np-distribution-type gamma \
      --gamma-shape 2.0 \
      --gamma-scale 3.0 \
      --np-min 2 \
      --np-max 50 \
      --subread-accuracy 0.85 \
      --threads 16
```

## Parameters

### Required Parameters

* `--amplicon-fasta`: Input FASTA file of amplicon sequences
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
* `--gamma-shape`: Gamma shape parameter if using `gamma` np distribution (default: 2.0)
* `--gamma-scale`: Gamma scale parameter (default: 3.0)
* `--np-min`: Minimum np value (default: 2)
* `--np-max`: Maximum np value (default: 50)
* `--subread-accuracy`: Subread accuracy for PBSIM (default: 0.85)
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
