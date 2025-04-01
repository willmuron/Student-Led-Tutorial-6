# Student-Led-Tutorial-6
# Task: 16S rRNA Amplicon Sequencing Analysis
# Date: April 3rd.

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
1. Go to your ocean folder
```
myocean
```
2. Fork and clone the "Student-Led-Tutorial-6" repository
3. Download FASTQ files of paired-end reads.
    -  Create a file named `download_sra.sh` with the following content:
    -  Type
```
vi download_sra.sh
```
   - Type `I` and enter the following workflow:
```bash
#!/bin/bash
# Download all SRR files from PRJNA759579
# esearch -db sra -query PRJNA759579 | efetch -format runinfo > runinfo.csv
# cut -d ',' -f 1 runinfo.csv | grep SRR > srr_accessions.txt

# Use prefetch to download all files
module load sra-toolkit
while read -r SRR; do
   echo "Downloading $SRR..."
   prefetch --max-size 100G $SRR 
   fastq-dump --gzip $SRR --split-files -O fastq_files/
done < srr_accessions.txt
```
- Exit the vi editor by typing `:wq`
- Run the `download_sra.sh` script
  
```
bash download_sra.sh
```
-  **Output**: Paired-end FASTQ files saved in the `fastq_files`/ directory.
---
## **Step 2: Data Quality Control**

1. Run fastqc
```
module load FastQC
mkdir fastqc_raw_data_out
fastqc -t 4 fastq_files/*.fastq.gz -o fastqc_raw_data_out
```
2. Run multiqc (will not work if you activate conda environment first, in such case `conda deactivate` before running multiqc)
```
module load python/3.8.6
multiqc --dirs fastqc_raw_data_out --filename multiqc_raw_data.html
```
3. Add, commit and push file to github repository. You can also transfer files using GLOBUS connect.
```
git add "multiqc_raw_data.html"
git commit -m "Adding multiqc results"
git push origin main
```
4. Go to web browser github and download `multiqc_raw_data.html`
5. Open file in web browser
---
## **Step 3: Generate Metadata File**
1. Metadata File Structure
The metadata file is a tab-separated values (TSV) file with the following columns:
- `#SampleID`: Unique sample identifier or SRA accession number
- `Group`: Experimental group or condition.
2. Example Metadata File
- Save the following as `metadata.tsv`:
  ``` Text
  #SampleID      Group
  SRRxxxxxx      Control
  SRRxxxxxx      Treatment
  SRRxxxxxx      Control
  SRRxxxxxx      Treatment
- Replace accession numbers and groups with actual information from the dataset. This file is given in the repository.
## **Step 3: Import Data into QIIME 2**
### Import the data:
1. We need to change our raw files names so that they can be read by qiime
  - Create directory for renamed files
```
mkdir -p casava_reads
```
  - Open the vi editor
    ```
    vi rename.sh
    ```
  - Type `I` to edit and write the following:
  ```
#!/bin/bash
for file in fastq_files/*_1.fastq.gz; do
    base=$(basename "$file" _1.fastq.gz)
    
    cp "fastq_files/${base}_1.fastq.gz" "casava_reads/${base}_S1_L001_R1_001.fastq.gz"
    cp "fastq_files/${base}_2.fastq.gz" "casava_reads/${base}_S1_L001_R2_001.fastq.gz"
done
```
- Run rename.sh script using:
```
bash rename.sh
```
2. Let's create the following directory
```         
mkdir reads_qza
```
3. You had installed the latest version of QIIME2, so you can either activate that environment with the command below.

```
module load anaconda3        
conda activate qiime2-amplicon-2024.2
```
- When you are finished with this tutorial make sure you deactivate the environment using `conda deactivate`

4. Feed qiime the raw reads
```         
qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-path casava_reads \
  --output-path reads_qza/reads.qza \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt
```
- cassava format, files follow a pattern `SampleID_L001_R1_001.fastq.gz`
---
5. remove primer sequences from reads, these are the primers used to enrich for a specific locus, e.g.:16S, COI, etc
```
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences reads_qza/reads.qza \
  --p-cores 4 \
  --p-front-f ACGCGHNRAACCTTACC \
  --p-front-r ACGGGCRGTGWGTRCAA \
  --p-discard-untrimmed \
  --p-no-indels \
  --o-trimmed-sequences reads_qza/reads_trimmed.qza
```
6. Visualize your data now
```         
qiime demux summarize \
  --i-data reads_qza/reads_trimmed.qza \
  --o-visualization reads_qza/reads_trimmed_summary.qzv
```
- Add, commit and push `reads_qza/reads_trimmed_summary.qzv` to git
```
git add "reads_qza/reads_trimmed_summary.qzv"
git commit -m "Adding reads summary results"
git push origin main
```
- Download and open `reads_trimmed_summary.qzv` in https://view.qiime2.org/
---

## **Step 4: Denoising with DADA2**
# Denoising reads
1. Join pair-end reds
```         
qiime vsearch merge-pairs \
  --i-demultiplexed-seqs reads_qza/reads_trimmed.qza \
  --output-dir reads_qza/reads_joined
```
2. Filter out low-quality reads
```         
qiime quality-filter q-score \
  --i-demux reads_qza/reads_joined/merged_sequences.qza \
  --o-filter-stats filt_stats.qza \
  --o-filtered-sequences reads_qza/reads_trimmed_joined_filt.qza
```
3. Summarize results
```         
qiime demux summarize \
  --i-data reads_qza/reads_trimmed_joined_filt.qza \
  --o-visualization reads_qza/reads_trimmed_joined_filt_summary.qzv
```
- Add, commit and push `reads_qza/reads_trimmed_joined_filt.qzv` to Git
- Download and open in https://view.qiime2.org/

4. Ask for some computer resources
```
salloc --mem=32G --time=2:00:00 --cpus-per-task=35 
```
6. Run DADA2 to denoise and generate ASVs:
```bash
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs reads_qza/reads.qza \
    --p-trim-left-f 0 --p-trim-left-r 0 \
    --p-trunc-len-f 250 --p-trunc-len-r 250 \
      --p-n-threads 32 \
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
#Download database to compare to
wget https://data.qiime2.org/2023.9/common/silva-138-99-nb-classifier.qza

qiime feature-classifier classify-sklearn \
   --i-classifier silva-132-99-515-806-nb-classifier.qza \
   --i-reads rep-seqs.qza \
   --p-n-jobs 32 \
   --output-dir taxa
```
2. If the above takes too long, please run:
```
ln -s /ocean/projects/agr250001p/shared/week-7-data/taxa .
```
---
### Filtering resultant table
1. Filter out rare ASVs
```         
qiime feature-table filter-features \
  --i-table feature-table.qza \
  --p-min-frequency 2 \
  --p-min-samples 1 \
  --o-filtered-table dada_table_filt.qza
```
- replace 2 and 1 with appropiate frequencies and observe the changes.
2. Filter out contaminant and unclassified ASVs
```         
qiime taxa filter-table \
  --i-table dada_table_filt.qza \
  --i-taxonomy taxa/classification.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table dada_table_filt_contam.qza
```
3. Subset and summarize filtered table
```         
qiime feature-table summarize \
  --i-table dada_table_filt_contam.qza \
  --o-visualization dada_table_filt_contam_summary.qzv
```
- Git add, commit and push `dada_table_filt_contam_summary.qzv` Download and open in https://view.qiime2.org/

4. Copy final table to current directory
```         
cp dada_table_filt_contam.qza dada_table_final.qza
```
5. Produce final table summary
```         
qiime feature-table filter-seqs \
  --i-data rep-seqs.qza \
  --i-table dada_table_final.qza  \
  --o-filtered-data rep_seqs_final.qza
```
```         
qiime feature-table summarize \
  --i-table dada_table_final.qza \
  --o-visualization dada_table_final_summary.qzv
```

- Git add, commit and push `dada_table_final_summary.qzv` Download and open in https://view.qiime2.org/

## **Visualization**
2. Visualize taxonomic composition:
```bash
qiime taxa barplot \
    --i-table feature-table.qza \
    --i-taxonomy taxa/classification.qza \
    --m-sample-metadata-file metadata.tsv \
    --o-visualization taxa-bar-plots.qzv
```
## **Step 6: Diversity Analysis** From 1-5 may take too long, think about producing files ahead of time.
1. Generate a phylogenetic tree using mafft. This software will be installed during tutorial 5
```
qiime tools export \
  --input-path rep-seqs.qza \
  --output-path exported_seqs
```
```bash
conda deactivate
conda create -n mafft-env -c bioconda mafft
conda activate mafft-env
```
```
mafft --thread 32 exported_seqs/dna-sequences.fasta > aligned-rep-seqs.fasta
```
```
conda deactivate
conda activate qiime2-amplicon-2024.2
```
```
qiime tools import \
  --type 'FeatureData[AlignedSequence]' \
  --input-path aligned_seqs.fasta \
  --output-path aligned-rep-seqs.qza
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
