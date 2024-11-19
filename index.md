# Bulk RNA-seq Quantification Pipeline 2023

## Overview

Bulk RNA sequencing (RNA-Seq) is a highly sensitive and accurate tool for meansuring expression across the transcriptome. In addition to the transcriptome quantification, RNA-Seq also allows researchers to detect new splicing junctions (e.g. TOPHAP/TOPHAP2-regtools), novel transcripts (e.g. Cufflinks), gene fusion (e.g. STAR-Fusion, Arriba), single nucleotide variants (e.g. STAR-GATK), and other features. **<u>This pipeline is for transcriptome quantification purpose only</u>.**

The current bulk RNA-Seq quantification methods can be grouped into two categories, **alignment-based** and **alignment-free**, as summarized in the table below. 

|            | Alignment-based methods                                      | Alignment-free methods                                       |
| ---------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Definition | Methods that quanfity from the **alignments to either transcritpome or genome** (i.e. BAM/SAM files). The coordinates of mapped reads are provided. | Methdos that quantify from **read-transcript matches**. The coordinates of mapped reads are **NOT** provided. |
| Principle  | "Seed and extend" for aligners.<br />Quantifiers vary in rules and weighting methods to count reads. | "*k-mer*s" based indexing;<br />Multiple models to identify read-transcript matches,<br /> e.g. SEMEs for  Salmon, T-DBG for Kallisto |
| Examples   | **Aligner**: Bowtie/Bowtie2, STAR, BWA, HISAT2, TopHat2, et. al.<br />**Quantifier**: RSEM, HTSeq, featureCounts, IsoEM, Cufflinks et. al. | Salmon, Kallisto, Sailfish, Fleximer, RNA-Skim, RapMap, et. al. |
| Accuracy   | High                                                         | a little bit lower or equal                                  |
| Speed      | Slow, a few hours for a typical run                          | Super-fast, a few minutes for a type run                     |

<u>**To ensure the accuracy of quantification, we employ one signature method from each of these two categories for cross-validation**:</u> **1)** **RSEM**, the most highly cited alignment-based method which shows the highest accuracy in most benchmarks; **2)** **Salmon**, one wicked-fast and highly-accurate alignment-free method which is recently further enhanced by integrating selective alignment and decoy sequences. We also introduce **3) STAR**, another alignment-based method recomended by GDC, as an optional method.

Below is an overview of the pipelines, which contains three sections: 1) **Preprocessing**: to prepare the standard inputs for quantification analysis; 2) **Quantification**: to **quantify** the gene and transcript expression in both alignment-based and alignment-free methods; 3) **Summarization**: to compile the **expression matrices** at both gene and transcript levels, and generate the **quanlity control report**.

![Picture1](/Users/qpan/Desktop/Picture1.png)

To serve better, we have:

* Complied all tools required into one single conda environment, which can be easily launched by:

  ```bash
  module load conda3/202210
  conda activate /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023
  ```

* Updated all the index libraries to the latest version

  | Genome | GENCODE | Path                                                         |
  | ------ | ------- | ------------------------------------------------------------ |
  | hg38   | v43     | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43 |
  | hg19   | v43     | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43 |
  | mm39   | vM32    | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32 |
  | mm10   | vM25    | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25 |

* Prepared the test cases of multi-format inputs

  | Cases   | Library Type | File Format           | Path                                                         |
  | ------- | ------------ | --------------------- | ------------------------------------------------------------ |
  | Sample1 | Paired-end   | FASTQ                 | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample1 |
  | Sample2 | Single-end   | FASTQ                 | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample2 |
  | Sample3 | Paired-end   | FASTQ, multiple lanes | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample3 |
  | Sample4 | Single-end   | FASTQ, multiple lanes | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample4 |
  | Sample5 | Paired-end   | BAM                   | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample5 |
  | Sample6 | Single-end   | BAM                   | /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample6 |

For traning purpose, we will go through the pipelines in a step-by-step way with real cases as listed above. 

## Preprocessing

**Why is preprocessing required?** The RNA-Seq data you start with could **vary in format** (e.g., FASTQ, BAM, FASTA) and usually **contain noisy sequences** (e.g., adapters leftovers, poor quality bases and other contaminations). So, we need to pre-process these data to generate the standard-in-format, clean-in-sequence FASTQ files which can be directly proceed to quantification analysis.

### 1. Prepare raw FASTQ files

Usually, you have two FASTQ files (R1, R2) for each paired-end sequencing sample (e.g. Sample1), or one single FASTQ file for each single-end sequencing sample (e.g. Sample2). If so, you are good to move forward.

However, if this is not the case, you will need to generate the raw FASTQ files by yourself:

* If you start with multiple FASTQ files for each mate, usually generated in different lanes, you need to **merge** them:

  ```bash
  ## For paired-end sequencing
  dir_sample3=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample3
  cat $dir_sample3/fqRaw_R1_L001.fq.gz $dir_sample3/fqRaw_R1_L002.fq.gz $dir_sample3/fqRaw_R1_L003.fq.gz $dir_sample3/fqRaw_R1_L004.fq.gz > /your_path/fqRaw_R1.fq.gz
  cat $dir_sample3/fqRaw_R2_L001.fq.gz $dir_sample3/fqRaw_R2_L002.fq.gz $dir_sample3/fqRaw_R2_L003.fq.gz $dir_sample3/fqRaw_R2_L004.fq.gz > /your_path/fqRaw_R2.fq.gz
  
  ## For single-end sequencing
  dir_sample4=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample4
  cat $dir_sample4/fqRaw_L001.fq.gz $dir_sample4/fqRaw_L002.fq.gz $dir_sample4/fqRaw_L003.fq.gz $dir_sample4/fqRaw_L004.fq.gz > /your_path/fqRaw.fq.gz
  ```

  > **NOTE:** The files of lanes **MUST BE IN THE SAME ORDER** between mate 1 and mate 2.

* If you start with BAM/SAM files, usually collected from other sources, you need to **convert** them:

  We previously used `**bedtools bamtofastq**` to convert BAM/SAM to FASTQ files, in which the BAM/SAM files MUST BE SORTED BY NAME. As recommended by GDC, we now move to `**bamtofastq**` integrated in the toolset named biobambam2. This command does not need the input BAM/SAM files sorted. The only thing you need to figure out before running it is that the input BAM/SAM file was generated from single-end or paired-end sequencing. This could be done by `samtools view -c -f 1 input.bam` which counts the matching records. It returns 0 for single-end sequencing or a non-zero integer for paired-end sequencing.
  
  ```bash
  ## to tell the BAM/SAM files are single- or paired-end
  samtools view -c -f 1 input.bam
  # This command counts the matching records in the bam/sam file, and returns 0 for single-end sequeing. Otherwise, the input bam/sam file is paired-end.
  
  ## For paired-end sequencing
  dir_sample5=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample5
  bamtofastq filename=$dir_sample5/rawBAM.toGenome.bam inputformat=bam gz=1 F=/your_path/fqRaw_R1.fq.gz F2=/your_path/fqRaw_R2.fq.gz # If the input file is in SAM format, change inputformat argument from bam to sam
  
  ## For single-end sequencing
  dir_sample6=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/testdata/sample6
  bamtofastq filename=$dir_sample6/rawBAM.toTranscriptome.bam inputformat=bam gz=1 S=/your_path/fqRaw.fq.gz # If the input file is in SAM format, change inputformat argument from bam to sam
  ```
  

After the step, you should have two FASTQ files for each paired-end sequencing sample, or one single FASTQ file for each single-end sequencing sample. If so, you are good to move forward.

### 2. Quality control

With the raw FASTQ file(s) prepared, a quality control analysis by FastQC is suggested to: 1) generate a html report which will give you a whole view of the FASTQ file quality; and 2) figure out some key informatioins which will be used in sebsequent analysis, e.g. encoding method of Phred quality score.

```bash
## For paired-end sequencing
fastqc /your_path/fqRaw_R1.fq.gz -o /your_path
fastqc /your_path/fqRaw_R2.fq.gz -o /your_path

## For single-end sequencing
fastqc /your_path/fqRaw.fq.gz -o /your_path
```

This analysis will generate a .html file with all metrics integerated. For more details of this report, please check out here: https://rtsf.natsci.msu.edu/sites/_rtsf/assets/File/FastQC_TutorialAndFAQ_080717.pdf

Here are some key infomations you may pay attention to:

* **Encoding of Phred Quality Score**:  it defines how the quality score of base calling is encoded. There are two options: Phred+33 or **Phred+64**. Phred+33 is denoted as **Sanger / Illumina 1.9** in the html report. This should always be the case for those sequenced in recent years. You will only find Phred+64 encoding (denoted as **Illumina 1.5 or lower**) on older data, which was sequenced several years ago. For more details: https://sequencing.qcfail.com/articles/incorrect-encoding-of-phred-scores/)
* **Sequenc Length**: The most common values are 46/45, 76/75, 101/100 or 151/150. The alignment of reads shorter than 50 could be tricky.
* **Adapter Content**: FastQC could detect several widely-used (**not all**) adapters, e.g. Illumina Universal Adapter (AGATCGGAAGAG), Illumina Small RNA 3' Adapter (TGGAATTCTCGG), Illumina Small RNA 5' Adapter (GATCGTCGGACT), Nextera Transposase Sequence (CTGTCTCTTATA) and SOLID Small RNA Adapter(CGCCTTGGCCGT). 

### 3. Adapter Trimming

Adapter trimming analysis trims not only the **adapter sequences**, but also the **sequences of unknown or low-quality bases**. It also discards the reads of **too-short length**. So, even though no significant adapter content was found in quality control analysis, it is still highly recommended to perform this analysis to  remove the low-quality reads.

We previous employed **Cutadapt** (https://cutadapt.readthedocs.io/en/stable/guide.html) for adapter trimming. It requires the sequence of adapter(s) and takes ~ two hours for a regular run. Now, we move to the **fastp** (https://github.com/OpenGene/fastp#adapters) which can automatically detect the adapter sequence(s) and trim much faster.

``` bash
## For paired-end sequencing
fastp -w 8 -l 30 -q 20 -n 5 -i /your_path/fqRaw_R1.fq.gz -I /your_path/fqRaw_R2.fq.gz -o /your_path/fqClean_R1.fq.gz -O /your_path/fqClean_R2.fq.gz

## For single-end sequencing
fastp -w 8 -l 30 -q 20 -n 5 -i /your_path/fqRaw.fq.gz -o /your_path/fqClean.fq.gz
```

**<u>Key arguments:</u>**

* <u>**-6/--phred64**: enable it if the input is using Phred+64 encoding. If enabled, fastp will automatically convert the Thread scores from Phread+64 to Phread+33. So the outputs of fastp are always encoded by Phread+33. </u>
* **-w/--thread**: number of threads to use concurrently.
* **-l/--length_required**: the trimmed reads shorter than this value will be discarded. The deault is 15, but 30 is recommended. The shorter reads tend to have multiple alignments, which may affect the quantification accuracy.
* **-q/--qualified_quality_phred:** the quality value that a base is qualified. The default is 15, but 20 is recommended.
* **-n/--n_base_limit**: the read/pair with more N bases will be discarded. The default is 5.

<u>**Key outputs:**</u>

* **fqClean.fq.gz**: FASTQ files with trimmed reads.
* **fastp.html**: an htmal report which summarizes some key matrices, e.g. read counts before and after trimming, insert size, base quality distribution, et. al.
* **fastp.json**: an json report with same matrices. This file provides the number of reads before and after adapter trimming used in the final QC report.

After the adapter timing, your FASTQ files are clean-in-sequence and can be directly proceed to quantification analysis.



## Alignment-free method

Salmon is a bulk RNA-Seq quantifier with **wicked-fast speed** and **comparable accuracy** (https://salmon.readthedocs.io/en/latest/salmon.html). It provides two working modes:

* **Mapping-based mode**: this is the feature mode that makes Salmon famous. Samlon employes a SEME (super maximal exact match)-based chaining algorithm to find and score potential mapping loci, making it super fast (because no alignment is needed) and comprably accurate. From version 1.0.0 and higher, Salmon introduced the **Selective Alignment** and **Decoy Sequences** to further improve its quantification accuracy (for more details: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8).
* **Alignment-based mode**: Salmon is not an aligner, so this mode, though called alignment-based, does NOT take FASTQ files. Instead, you can simpley provide: 1) BAM/SAM files of alignments generated by other aligners, e.g., STAR, Bowtie; and 2) FASTA file of reference transcriptome. No indexing is required for this mode.

In this pipeline, we suggest to use the mapping-based mode ONLY.

```bash
## Indexing directories
index_hg38v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43/bulkRNAseq/Salmon/index_decoy
index_hg19v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43/bulkRNAseq/Salmon/index_decoy
index_mm39vM32=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32/bulkRNAseq/Salmon/index_decoy
index_mm10vM25=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25/bulkRNAseq/Salmon/index_decoy

## transcript-to-gene mapping files
tr2gene_hg38v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43/gencode.v43.primary_assembly.annotation.gtf
tr2gene_hg19v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43/gencode.v43lift37.annotation.gtf
tr2gene_mm39vM32=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32/gencode.vM32.primary_assembly.annotation.gtf
tr2gene_mm10vM25=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25/gencode.vM25.primary_assembly.annotation.gtf

## Quantification by Salmon
# For paired-end sequencing
salmon quant -i $hg38v43_index -l A -p 8 -g $tr2gene_hg38v43 -1 /your_path/fqClean_R1.fq.gz -2 /your_path/fqClean_R2.fq.gz --validateMappings -o /your_path/quantSalmon
# For single-end sequencing
salmon quant -i $hg38v43_index -l A -p 8 -g $tr2gene_hg38v43 -r /your_path/fqClean.fq.gz --validateMappings -o /your_path/quantSalmon
```

**<u>Key arguments:</u>**

Salmon is very smart, It could learn most of the arguments by itseft. So there are only a few arguments you need to specify, making it super easy to use.

* **-l/--libType**: library type. Salmon employs a three-letter string to denote the library type (https://salmon.readthedocs.io/en/latest/library_type.html#fraglibtype): 1) the relative orientation of two matched mates, including I = inward, O = outward and M = matching; 2) the protocol is stranded or unstranded, including S = stranded and U = unstranded; 3) the strand from which the read originates in a strand-specific protocol, including F = read 1 (or single-end read) from the forward stand and R =  read 1 (or single-end read) from the reverse stand. As a result, there are 9 library types: **3 for unstranded library (IU, MU, OU)** and **6 for stranded library (ISF, ISR, MSF, MSR, OSR, OSF)**, as shown below. Unless you know the exact library type of your samples, just use "A", which asks Salmon to figure it out automatically.

  ![_images/ReadLibraryIllustration.png](https://salmon.readthedocs.io/en/latest/_images/ReadLibraryIllustration.png)

* -i/--index: salmon index files. We have generated them for hg38, hg19, mm39 and mm10. You don't need to regenerate them.

* -g/--geneMap: file contaning a mapping of transcript to gene. **This must be provided to genereate the gene-level abundance estimates.**

* -p/--threads: number of threads to use concurrently.

<u>**Key outputs**</u>: (for more details: https://salmon.readthedocs.io/en/latest/file_formats.html#fileformats)

* **quant.sf:** contains transcript-level abundance estimates, TPM and raw counts.

* **quant.genes.sf**: contains the aggregated gene-level abundance estimates, TPM and raw counts.

* **lib_format_counts.json**: this file contains a count of the number of mappings that matched each possible library type. You can find the infered library type from the line of "expected_format". Here is a summary of library types used in different tools:

  | Salmon (--libType, SE/PE) | RSEM (--strandedness) | TopHap (--library-type) | HTSeq (--stranded) |
  | ------------------------- | --------------------- | ----------------------- | ------------------ |
  | U/IU                      | none                  | -fr-unstranded          | no                 |
  | SR/ISR                    | reverse               | -fr-firststrand         | reverse            |
  | SF/ISF                    | forward               | -fr-secondstrand        | yes                |

  

## Alignment-based methods

RSEM (RNA-Seq by Expectation Maximization) is an accurate and user-friendly software tool for quantifying transcript and gene abundances from RNA-seq data. Here is the detailed introduction of RSEM: https://github.com/bli25/RSEM_tutorial.

* RSEM is a quantifier, not an aligner. RSEM can directly take the FASTQ files, but it does not align the reads by itself. Instead, it empolys the Bowtie2 (by default) or STAR/HISAT2 (optional) for reads alignment.
* RSEM doesn't rely on the existence of a reference genome. It quantifies from the alignments to reference transcriptome.
* RSEM is famouse for its ability to **effectively use ambiguously-mapping reads**. This is the main reason for its high quantification accuracy.

### 1. Bowtie2-RSEM (default)

In this pipeline, Bowtie2-RSEM is set as the default. 

```bash
## Indexing directories
index_hg38v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43/bulkRNAseq/RSEM
index_hg19v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43/bulkRNAseq/RSEM
index_mm39vM32=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32/bulkRNAseq/RSEM
index_mm10vM25=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25/bulkRNAseq/RSEM

## Transcript bins
bins_hg38v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43/bulkRNAseq/genebodyBins
bins_hg19v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43/bulkRNAseq/genebodyBins
bins_mm39vM32=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32/bulkRNAseq/genebodyBins
bins_mm10vM25=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25/bulkRNAseq/genebodyBins

## For paired-end sequencing
rsem-calculate-expression --num-threads 8 \
--bowtie2 --bowtie2-path /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/bin \
--bowtie2-sensitivity-level sensitive \
--strandedness reverse --phred33-quals \
--sort-bam-by-coordinate \
--paired-end /your_path/fqClean_R1.fq.gz /your_path/fqClean_R2.fq.gz \
$index_hg38v43/index_bowtie2 /your_path/quantRSEM

## For single-end sequencing
rsem-calculate-expression --num-threads 8 \
--bowtie2 --bowtie2-path /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/bin \
--bowtie2-sensitivity-level sensitive \
--strandedness reverse --phred33-quals \
--sort-bam-by-coordinate \
/your_path/fqClean.fq.gz \
$index_hg38v43/index_bowtie2 /your_path/quantRSEM

## Prepare gene body coverage statistics
bedtools multicov -bams /your_path/quantRSEM.transcript.sorted.bam -bed $bins_hg38v43/genebodyBins_HouseKeepingTranscripts.txt > /your_path/genebodyCoverage.txt
```

**<u>Key arguments:</u>**

There is only one key argument you need to manually specify, `**--strandedness [none|forward|reverse]**`. You can figure it out from the **lib_format_counts.json** returned by Salmon, as mentioned before.

**<u>Key outputs</u>**: (for more details: https://salmon.readthedocs.io/en/latest/file_formats.html#fileformats)

* **quantRSEM.isoforms.results**: contains transcript-level abundance estimates, raw counts, TPM and FPKM.
* **quantRSEM.genes.results**: contains the aggregated gene-level abundance estimates, raw counts, TPM and FPKM.
* **quantRSEM.transcript.bam/quantRSEM.transcript.sorted.bam**: BAM files with alignments to reference transcriptome. These files will be used to calculate gene body coverage and others.
* **quantRSEM.stat/quantRSEM.cnt**: contains alignment statistics based purely on the alignment results obtained from aligners. This file provides the number of totally-mapped and uniquely-mapped reads used in the final QC report.

### 2. STAR-RSEM (optional)

Following the RSEM tutorial, we also provide the codes for STAR-RSEM **as optional**. STAR (Spliced Transcripts Alignment to a Reference) was uniquely designed for RNA-Seq alignment. It's ultrafast and it does support the spliced alignment. For more details: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/.

```bash
## Indexing directories
index_hg38v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43/bulkRNAseq/RSEM
index_hg19v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43/bulkRNAseq/RSEM
index_mm39vM32=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32/bulkRNAseq/RSEM
index_mm10vM25=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25/bulkRNAseq/RSEM

## For paired-end sequencing
rsem-calculate-expression -num--threads 8 \
--star --star-path /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/bin \
--star-gzipped-read-file \
--strandedness reverse --phred33-quals \
--sort-bam-by-coordinate \
--paired-end /your_path/fqClean_R1.fq.gz /your_path/fqClean_R2.fq.gz \
$index_hg38v43/index_star /your_path/quantRSEM

## For single-end sequencing
rsem-calculate-expression -num--threads 8 \
--star --star-path /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/bin \
--star-gzipped-read-file \
--strandedness reverse --phred33-quals \
--sort-bam-by-coordinate \
/your_path/fqClean.fq.gz \
$index_hg38v43/index_star /your_path/quantRSEM
```

There are two key arguments you need to specify: 1) `**--strandedness [none|forward|reverse]**`, which is same with Bowtie2-RSEM; and 2) `**--star-gzipped-read-file**`, which should be enabled if the input FASTQ files are gzipped.

The outputs are exactly same with those from Bowtie2-RSEM.

### 3. STAR-HTSeq (optional)

In addtion RSEM, HTSeq (https://htseq.readthedocs.io/en/master/) is another popular quantifier for RNA-Seq data. Compared with RSEM which employed the Expectaton Maximization-algorithsm to handle the ambiguously-mapped reads, the way HTSeq counts reads is simpler and more transparent: as shown in the figure below, **HTSeq does NOT count ambiguously-mapped reads and it does NOT weight the reads by either their overlaps with refence genes or mapping scores.**

![_images/count_modes.png](https://htseq.readthedocs.io/en/master/_images/count_modes.png)

```bash
## Indexing directories
index_hg38v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43/bulkRNAseq/STAR/index_overhang100
index_hg19v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43/bulkRNAseq/STAR/index_overhang100
index_mm39vM32=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32/bulkRNAseq/STAR/index_overhang100
index_mm10vM25=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25/bulkRNAseq/STAR/index_overhang100

## Alignment by STAR: for paired-end sequencing
STAR --genomeDir $index_hg38oh100 --readFilesIn /your_path/fqClean_R1.fq.gz /your_path/fqClean_R2.fq.gz --readFilesCommand zcat \
--runThreadN 8 --twopassMode Basic --quantMode TranscriptomeSAM \
--outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif \
--outSAMattributes All --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --outSAMattrRGline ID:SampleName SM:SampleName LB:Illumina PL:Illumina PU:Illumina --outFileNamePrefix /your_path/quantSTAR

## Alignment by STAR: for single-end sequencing
STAR --genomeDir $index_hg38oh100 --readFilesIn /your_path/fqClean.fq.gz --readFilesCommand zcat \
--runThreadN 8 --twopassMode Basic --quantMode TranscriptomeSAM \
--outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif \
--outSAMattributes All --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --outSAMattrRGline ID:SampleName SM:SampleName LB:Illumina PL:Illumina PU:Illumina --outFileNamePrefix /your_path/quantSTAR

## Quantification by HTSeq
htseq-count -f bam -r pos -s reverse -i gene_id -m intersection-nonempty -n 8 /your_path/quantSTAR/Aligned.out.bam $gtf_hg38 > /your_path/quantSTAR/htseq_counts.txt
```

**<u>Key arguments</u>**:

In this pipeline, STAR was employed as the aligner, with the arguments recommended by GDC (https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#mrna-expression-transformation). The only key argument for STAR you need to specify is the `**--outSAMattrRGline**`.

For HTSeq, there are two key arguments to be specified: 1) `**-r/--order [pos|name]**`, which is the way BAM files were sorted; and 2) `**-s/--stranded [no|yes|reverse]**`, which is the library type.

**<u>Key outputs:</u>**

* **htseq_counts.txt**: contains gene-level abundance estimates, raw counts only.
* **Aligned.sortedByCoord.out.bam**: BAM file with alignments to reference genome, sorted by coordinates.
* **Aligned.toTranscriptome.out.bam**: BAM file with alignments to reference transcriptome.
* **Log.final.out**: contains alignment statistics generated by STAR. This file provides the number of totally-mapped and uniquely-mapped reads used in the final QC report.



## Quantification Summary

For each single sample, we will generate a QC report from the quantification results by Salmon and RSEM. There are five QC metrics in each report:

* Alignment statistics: showing read counts and mapping rates et. al.;
* Quantification statistics: showing the number of genes/transcripts identified by two methods and the correlations of quantification results of them;
* Biotype distribution: showiing the composition of types of identified transcripts and genes;
* Quantification accuracy: correlations of gene expression measurements by Salmon and RSEM;
* Genebody coverage statistics: showing if the RNA samples were degraded or not.

![image-20230901163554962](/Users/qpan/Library/Application Support/typora-user-images/image-20230901163554962.png)

```sh
Rscript -e "rmarkdown::render(input = '/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/bin/summaryIndividual.Rmd', clean = TRUE, quiet = F, output_format = 'html_document', output_file = 'quantSummary.html', output_dir = '/your_path/quantSummary', params = list(sampleName = 'sample1', quant_dir = '/your_path'))"
```

As shown above, we built a R markdown script, **summaryIndividual.Rmd**, to generate this report from the outputs of quantification pipelines. The only argument you need to pay attention to is the **quant_dir**. The script will call these files to generate the QC report:

* quant_dir/preProcessing/adapterTrimming.json
* quant_dir/quantRSEM/quant.isoforms.results
* quant_dir/quantRSEM/quant.genes.results
* quant_dir/quantRSEM/quant.stat/quant.cnt
* quant_dir/quantRSEM/geneCoverage.txt
* quant_dir/quantSalmon/quant.sf
* quant_dir/quantSalmon/quant.genes.sf

As for the outputs, the script will generate a .html file named quantSummary.html which summarizes all five QC metrics. And a few .txt files and .pdf files will be generated as well to provide the source data of the html QC report. These files can be also used to generate the multi-sample QC report if multiple samples were sequenced.



## Best Practice

The step-by-step tutorial shown above is set to get you familiar with this pipeline. In practice, it could be tedious to run this pipeline in a step-by-step manor, especially in the cases with multiple samples. For the best practice, we designed **a sample table-centered strategy** to run this pipeline.

### 1. Sample Table

As shown below, sample table is a tab-delimited text file with 6 columns:

* **sampleID**: name of samples, which can only contain letters, numbers and underscores. They should NOT start with numbers.
* **libraryType**: [PE | SE]. For the alignment files, like BAM or SAM, please use samtools view -c -f 1 input.bam to tell. See above for how.
* **phredMethod**: [Phred33 | Phred64]. Those data generated in recent 5 years should be Phred33. FastQC can tell it. See above for how.
* **reference**: [hg38 | hg19 | mm39 | mm10]. Reference genomes supported by this pipeline.
* **input**: files for quantification. This pipeline supports FASTQ, BAM and SAM.
* **output**: path to save the output files.

![image-20230901180916287](/Users/qpan/Library/Application Support/typora-user-images/image-20230901180916287.png)

That's it! This table is the only one you need to prepare. You can generate it in either of these ways:

* any coding language you prefer, e.g. BASH, R, Python, Perl et. al.
* Excel or VIM. And for you convince, we have a templete avaible on HPC: /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/scripts/temp_bulkRNAseq.txt. You can simply copy it to your folder and edit it in VIM. (To insert TAB in VIM: on INSERT mode press control + v + TAB)

### 2. To run the pipeline

```bash
## 0. activate the conda env
module load conda3/202210
conda activate /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023

## 1. preprocessing
/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/scripts/preProcessing.pl sampleTable.txt

## 2. adapterTrimming
/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/scripts/adapterTrimming.pl sampleTable.txt

## 3. quantify by Salmon
/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/scripts/quantSalmon.pl sampleTable.txt

## 4. genebody coverage
/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/scripts/geneCoverage.pl sampleTable.txt

## 5. quantification summary
/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq_2023/git_repo/scripts/quantSummary.pl sampleTable.txt

```







 



