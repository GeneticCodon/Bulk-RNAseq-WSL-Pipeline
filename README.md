# Bulk RNA-seq Pipeline (Windows + WSL)

A comprehensive guide for installing, configuring, and running a bulk RNA-seq workflow on **Windows** (with or without WSL) using open-source tools. This documentation covers tool installation, genome index setup, typical workflow commands, and downstream R analysis.

---

## Table of Contents

1. [Overview & Workflow](#overview--workflow)
2. [Pipeline Tools](#pipeline-tools)
3. [Installation Instructions](#installation-instructions)
    - [a. WSL (Ubuntu)](#a-wsl-ubuntu)
    - [b. Native Windows](#b-native-windows)
4. [Genome Downloading & Indexing](#genome-downloading--indexing)
5. [Typical Workflow Commands](#typical-workflow-commands)
6. [Downstream R Analysis](#downstream-r-analysis)
7. [References & Resources](#references--resources)

---

## Overview & Workflow

This pipeline is designed for robust, reproducible RNA-seq data analysis, from quality control to differential expression. It is suitable for both beginners and experienced users, and can be run entirely on Windows (via WSL or natively).

**Workflow Summary:**

1. **Quality Control:** FastQC → MultiQC
2. **Read Trimming:** Fastp
3. **Alignment:** HISAT2
4. **Post-processing:** Samtools / Sambamba
5. **Quantification:** StringTie / featureCounts
6. **Differential Expression:** R (DESeq2, EnhancedVolcano)

---

## Pipeline Tools

| Tool                      | Purpose                              | Resource Link                                                                              |
|---------------------------|--------------------------------------|-------------------------------------------------------------------------------------------|
| FastQC                    | Raw read quality check               | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)                      |
| MultiQC                   | Aggregate QC reports                 | [MultiQC](https://multiqc.info/)                                                          |
| Fastp                     | Read trimming & quality control      | [fastp GitHub](https://github.com/OpenGene/fastp)                                         |
| HISAT2                    | Indexing & alignment                 | [HISAT2 GitHub](https://github.com/DaehwanKimLab/hisat2)                                  |
| Sambamba                  | BAM deduplication/manipulation       | [Sambamba GitHub](https://github.com/biod/sambamba)                                       |
| StringTie                 | Transcript assembly/quantification   | [StringTie](https://ccb.jhu.edu/software/stringtie/)                                      |
| featureCounts (Subread)   | Count matrix generation              | [Subread](https://subread.sourceforge.net/)                                               |
| R (DESeq2, EnhancedVolcano)| Differential expression & visualization | [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [EnhancedVolcano](https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html) |

---

## Installation Instructions

### a. WSL (Ubuntu Recommended)

1. **Enable WSL** and install Ubuntu from the Microsoft Store.
2. **Update system:**
    ```sh
    sudo apt update && sudo apt upgrade -y
    ```

3. **Install pipeline tools:**
    ```sh
    # FastQC
    sudo apt install fastqc -y

    # Python3 + pip (for MultiQC)
    sudo apt install python3-pip -y
    pip install multiqc

    # Fastp
    sudo apt install fastp -y

    # HISAT2
    sudo apt install hisat2 -y

    # Sambamba
    sudo apt install sambamba -y

    # StringTie
    sudo apt install stringtie -y

    # featureCounts (Subread)
    sudo apt install subread -y

    # R
    sudo apt install r-base -y
    ```

4. **Install R packages:**  
    Inside the R console:
    ```r
    install.packages("BiocManager")
    BiocManager::install("DESeq2")
    BiocManager::install("EnhancedVolcano")
    ```

---

### b. Native Windows

> **Note:** Add each tool's binary folder to your Windows `PATH` ([Microsoft guide](https://learn.microsoft.com/en-us/previous-versions/office/developer/sharepoint-2010/ee537574(v=office.14))).

- **FastQC:** Download `.zip` from [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), unzip, add to PATH.
- **MultiQC:** Install via Python (after installing [Python for Windows](https://www.python.org/downloads/)):
    ```sh
    pip install multiqc
    ```
- **Fastp:** Download precompiled `.exe` from [releases](https://github.com/OpenGene/fastp/releases), add to PATH.
- **HISAT2:** Download `.zip` from [releases](https://github.com/DaehwanKimLab/hisat2/releases), unzip, add to PATH.
- **Sambamba:** Download `.exe` from [releases](https://github.com/biod/sambamba/releases).
- **StringTie:** Download binaries from [StringTie](https://ccb.jhu.edu/software/stringtie/).
- **Subread (featureCounts):** Download Windows binary from [Subread](https://sourceforge.net/projects/subread/files/).
- **R:** Install from [CRAN](https://cran.r-project.org/).  
  Then install DESeq2 + EnhancedVolcano as above.

---

## Genome Downloading & Indexing

Example: **Human genome (GRCh38)**

```sh
# Create genome folder
mkdir -p ~/genomes/hg38 && cd ~/genomes/hg38

# Download genome FASTA
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Download GTF annotation
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz

# Indexing with HISAT2
hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa hg38_index
```
> This generates `hg38_index.*.ht2` files for alignment.
> “Alignment with HISAT2 may need ~8GB RAM for the human genome.”
> “If you run into memory errors, try mouse genome first.”

---

## Typical Workflow Commands

```sh
# 1. Quality Control
fastqc sample_R1.fastq.gz sample_R2.fastq.gz -o fastqc_reports
multiqc fastqc_reports -o multiqc_report

# 2. Read Trimming
fastp -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
      -o clean_R1.fastq.gz -O clean_R2.fastq.gz \
      -h fastp.html -j fastp.json

# 3. Alignment
hisat2 -x ~/genomes/hg38/hg38_index \
       -1 clean_R1.fastq.gz -2 clean_R2.fastq.gz \
       -S sample.sam

# 4. Convert SAM to BAM, sort, deduplicate
samtools view -bS sample.sam | samtools sort -o sample_sorted.bam
sambamba markdup sample_sorted.bam sample_dedup.bam

# 5. Quantification
stringtie sample_dedup.bam -G Homo_sapiens.GRCh38.115.gtf -o sample.gtf -A sample_gene_abundance.tab
# or using featureCounts
featureCounts -a Homo_sapiens.GRCh38.115.gtf -o counts.txt sample_dedup.bam
```

---

## Downstream R Analysis

```r
# Load libraries
library(DESeq2)
library(EnhancedVolcano)

# Import count matrix (edit path if needed)
cts <- as.matrix(read.csv("counts.txt", sep="\t", row.names=1))
coldata <- data.frame(condition=factor(c("treated","control","treated","control")))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds)

# Volcano plot
EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue')
```

---

## References & Resources

- [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC Documentation](https://multiqc.info/docs/)
- [fastp Documentation](https://github.com/OpenGene/fastp)
- [HISAT2 Manual](https://daehwankimlab.github.io/hisat2/)
- [Sambamba Manual](https://lomereiter.github.io/sambamba/)
- [StringTie Manual](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)
- [Subread/featureCounts User Guide](https://subread.sourceforge.net/)
- [DESeq2 Bioconductor](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [EnhancedVolcano Bioconductor](https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html)
- [Adding to Windows PATH](https://learn.microsoft.com/en-us/previous-versions/office/developer/sharepoint-2010/ee537574(v=office.14))

---

**Maintainer**: [GeneticCodon](https://github.com/GeneticCodon)  
**License**: MIT 

---

## Hall of Testers (Windows + WSL)
Want your name here? Comment on this issue: https://github.com/GeneticCodon/Bulk-RNAseq-WSL-Pipeline/issues/1

- @username — Windows 11 + WSL2 Ubuntu 22.04 — “fastp worked out of the box”
- @username — Windows 10 + WSL2 Ubuntu 20.04 — “needed to fix PATH for HISAT2”


<a href="https://github.com/GeneticCodon/Bulk-RNAseq-WSL-Pipeline/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=GeneticCodon/Bulk-RNAseq-WSL-Pipeline" />
</a>

