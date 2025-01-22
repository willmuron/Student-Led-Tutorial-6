# Student-Led-Tutorial-6
# Task: 16S rRNA Amplicon Sequencing Analysis
# Date: April 3rd

## **Objective**
Students will:
1. Perform quality control and trimming of 16S rRNA sequencing data.
2. Use **QIIME 2** for denoising, taxonomic assignment, and diversity analysis.
3. Visualize and interpret diversity metrics and taxonomic composition.

---

## **Software and Manuals**
### **Required Software**
1. **QIIME 2**: End-to-end amplicon sequencing analysis.
   - [QIIME 2 Website](https://qiime2.org/)
   - [QIIME 2 Documentation](https://docs.qiime2.org/)
2. **FastQC**: Quality control of raw reads.
   - [FastQC Website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
3. **Trimmomatic**: Read trimming and quality filtering.
   - [Trimmomatic GitHub](http://www.usadellab.org/cms/?page=trimmomatic)
   - [Trimmomatic Manual](http://www.usadellab.org/cms/uploads/supplementary/TrimmomaticManual_v0.39.pdf)
4. **R and RStudio**: Visualization and statistical analysis.
   - [R Project](https://www.r-project.org/)
   - [RStudio Website](https://posit.co/downloads/)

---

## **Step 1: Download SRR Data**
### Dataset: 16S rRNA Sequencing Data
We will use publicly available 16S rRNA sequencing data from a microbiome study.  
- **Source**: [NCBI BioProject PRJNA759579](https://www.ebi.ac.uk/ena/browser/view/PRJEB43384](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=759579)  
  - Download FASTQ files of paired-end reads.
    -  Create a file named `download_sra.sh` with the following content:
       ```bash
       #!/bin/bash
       # Download all SRR files from PRJNA759579
       esearch -db sra -query PRJNA759579 | efetch -format runinfo > runinfo.csv
       cut -d ',' -f 1 runinfo.csv | grep SRR > srr_accessions.txt

       # Use prefetch to download all files
       while read -r srr; do
              echo "Downloading $srr..."
              prefetch $srr
              fasterq-dump $srr --split-files -O fastq_files/
       done < srr_accessions.txt
        ```
    -  Run the `download_sra.sh` script
        ```bash
        bash download_sra.sh
      -  **Output**: Paired-end FASTQ files saved in the `fastq_files`/ directory.
---
## **Step 2: Generate Metadata File**
1. Metadata File Structure
The metadata file is a tab-separated values (TSV) file with the following columns:
- `#SampleID`: Unique sample identifier.
- `forward-absolute-filepath`: Full path to the forward read FASTQ file.
- `reverse-absolute-filepath`: Full path to the reverse read FASTQ file.
- `Group`: Experimental group or condition.
2. Example Metadata File
- Save the following as `metadata.tsv`:
  ``` Text
  #SampleID    forward-absolute-filepath              reverse-absolute-filepath              Group
  Sample1      /path/to/fastq_files/SRR123456_1.fastq /path/to/fastq_files/SRR123456_2.fastq Control
  Sample2      /path/to/fastq_files/SRR123457_1.fastq /path/to/fastq_files/SRR123457_2.fastq Treatment
  Sample3      /path/to/fastq_files/SRR123458_1.fastq /path/to/fastq_files/SRR123458_2.fastq Control
  Sample4      /path/to/fastq_files/SRR123459_1.fastq /path/to/fastq_files/SRR123459_2.fastq Treatment
- Replace file paths and groups with actual information from the dataset.
## **Step 3: Import Data into QIIME 2**
1. Import the data:
    ```bash
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path metadata.tsv \
        --output-path paired-end-demux.qza \
        --input-format PairedEndFastqManifestPhred33V2
2. Visualize quality plots:
    ```bash
    qiime demux summarize \
        --i-data paired-end-demux.qza \
        --o-visualization demux.qzv
## **Step 4: Denoising with DADA2**
1. Run DADA2 to denoise and generate ASVs:
```bash
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs paired-end-demux.qza \
    --p-trim-left-f 0 --p-trim-left-r 0 \
    --p-trunc-len-f 250 --p-trunc-len-r 250 \
    --o-table feature-table.qza \
    --o-representative-sequences rep-seqs.qza \
    --o-denoising-stats denoising-stats.qza
```
2. Summarize the feature table:
```bash
qiime feature-table summarize \
    --i-table feature-table.qza \
    --o-visualization feature-table.qzv \
    --m-sample-metadata-file metadata.tsv
```
## **Step 5: Taxonomic Assignment**
1. Assign taxonomy using a pre-trained classifier (e.g., SILVA):
```bash
qiime feature-classifier classify-sklearn \
    --i-classifier silva-132-99-515-806-nb-classifier.qza \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza
```
2. Visualize taxonomic composition:
```bash
qiime taxa barplot \
    --i-table feature-table.qza \
    --i-taxonomy taxonomy.qza \
    --m-sample-metadata-file metadata.tsv \
    --o-visualization taxa-bar-plots.qzv
```
## **Step 6: Diversity Analysis**
1. Generate a phylogenetic tree
```bash
qiime alignment mafft \
    --i-sequences rep-seqs.qza \
    --o-alignment aligned-rep-seqs.qza
```
2. Remove/Mask highly variable positions that might introduce noise:
```bash
qiime alignment mask \
    --i-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza
```
3. Build unrooted tree
```bash
qiime phylogeny fasttree \
    --i-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza
```
4. Root the tree
```bash
qiime phylogeny midpoint-root \
    --i-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza
```
5.Calculate alpha and beta diversity metrics:
```bash
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree.qza \
    --i-table feature-table.qza \
    --p-sampling-depth 10000 \
    --m-sample-metadata-file metadata.tsv \
    --output-dir core-metrics-results
```
6. Generate PCoA plots:
```bash
qiime emperor plot \
    --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
    --m-sample-metadata-file metadata.tsv \
    --o-visualization unweighted-unifrac-emperor.qzv
```
## **Step 7: Visualization and Interpretation**

1. **Open `.qzv` Files**  
   Use [QIIME 2 View](https://view.qiime2.org) to explore the following visualizations:
   - Quality control summaries.
   - Taxonomic bar plots.
   - Diversity metrics (alpha and beta diversity).

2. **Interpret Results**  
   - Identify dominant taxa at different taxonomic levels.
   - Compare diversity between experimental groups.
