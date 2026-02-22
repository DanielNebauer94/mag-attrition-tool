# MAG Attrition Tool

A Python-based pipeline for quantifying and visualising lineage-level read attrition during MAG-focused metagenomic workflows.

The tool maps reads against MAG bins, separates mapped and unmapped read fractions, integrates Kraken2 classifications, and generates interactive visualisations.

---

### Functionality

The workflow consists of two scripts:

#### `MAG_attrition_pipeline.py`

- Concatenates MAG bin FASTA files into a reference
- Maps reads using **BWA**
- Separates mapped and unmapped reads with **samtools**
- Integrates Kraken2 classifications
- Produces lineage-level attrition summaries (including average GC content)

#### `make_attribute_dotplot_html.py`

- Scans pipeline results for `*_attrition_lineage*.csv` files
- Generates GC (%) vs mapped-read proportion (%) dotplots
- Colours points by phylum and scales point size by abundance
- Exports interactive Plotly HTML visualisations (per-sample and combined)

---

### Requirements

#### External software (must be on `PATH`)

- `bwa`
- `samtools`

#### Python packages

- `pandas`
- `plotly`

---

### Usage example

#### Run pipeline

```bash
python MAG_attrition_pipeline.py \
  --sheet samples.csv \
  --outdir results_mag_attrition \
  --threads 16
```

**Optional arguments:**

- `--classified-only` - Exclude unclassified Kraken2 reads
- `--min-reads` - Minimum reads per lineage (default: 5)
- `--cache-dir` - Directory for cached references (default: `.cache_binned_refs`)
- `--keep-mapping` - Retain BAM intermediates

---

#### Generate visualisations

```bash
python make_attribute_dotplot_html.py \
  --results_dir results_mag_attrition
```

**Optional arguments:**

- `--min_reads` - Filter low-abundance lineages
- `--outname_all` - Custom filename for combined plot

---

### Input data

#### Sample sheet (CSV)

The pipeline requires a CSV sample sheet with the following columns:

- `sample_id`
- `bins_dir`
- `r1_fastq`
- `r2_fastq`
- `kraken2_report_tsv`
- `kraken2_tsv`

**Example:**

```csv
sample_id,r1_fastq,r2_fastq,,kraken2_tsv,kraken2_report_tsv,bins_dir
S1,reads_R1.fastq.gz,reads_R2.fastq.gz,classifications.kraken2,kraken2_report.tsv,/path/to/bins
```

#### Bins directory

Directory containing MAG FASTA files.

Supported extensions:

- `.fa`
- `.fasta`
- `.fna`

All FASTA files are concatenated into a reference used for read mapping.  
A cached reference is generated automatically to avoid recomputation.

---

#### FASTQ files

- Paired-end reads (`R1` / `R2`)
- Files may be gzipped (`.fastq.gz`)

---

#### Kraken2 outputs

- Kraken2 report file (`kraken2_report_tsv`)
- Kraken2 classification file (`kraken2_tsv`)

These files are used to assign reads to taxonomic lineages.

---

### Outputs

Per sample:

- `mapping/` (intermediate mapping files)
- `mapped_R1.fq.gz`, `mapped_R2.fq.gz`
- `unmapped_R1.fq.gz`, `unmapped_R2.fq.gz`
- `sample_id_attrition_lineage.csv`
- `interactive_attrition_SAMPLE.html`

Combined:

- `interactive_attrition_ALL.html`

---

### Plot interpretation

- **X-axis:** GC content (%)
- **Y-axis:** mapped read proportion (%)

Visual encoding:

- Colour → Phylum
- Size → Lineage abundance

---
