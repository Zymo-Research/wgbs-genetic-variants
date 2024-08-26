# Variant Calling from Whole Genome Bisulfite Sequencing Data (WGBS): From FASTQ to VCF

This tutorial provides a step-by-step guide on how to perform variant calling on Whole Genome Bisulfite Sequencing (WGBS) data, from BAM files to obtaining a VCF file. We will use the following tools:

1. **Revelio** to convert WGBS BAM files to genomics-compatible BAM files.
2. **FreeBayes** for variant calling.

For the purposes of this tutorial, we already assume your WGBS data has been appropriately trimmed, aligned (for example, using a tool such as Bismark), and QC-checked. We will begin with an aligned BAM file and demonstrate how to obtain called variants in a VCF (Variant Call Format) file.

## Conceptual Introduction

We have elected to use a tool called Revelio, which transforms a BAM file from WGBS data into one that is compatible with the same tools as any genomics analysis. The advantage is that standard variant calling tools can now be used on WGBS data. Also, if a researcher is performing both WGBS and whole genome sequencing (WGS, or DNA-seq), the same variant calling pipeline could be used to compare between the datasets. Revelio achieves this by transforming quality scores based on the likelihood that a C->T conversion is due to a genetic change or a bisulfite conversion. It is able to do this, for stranded (directional), WGBS libraries by leveraging unconverted information present on the second strand. The quality score conversion will effectively mask any bisulfite conversion events, but will allow the variant caller to analyze any other mutation events. See the Revelio publication for more information on how the algorithm works.

## Prerequisites

* Compute node suitable for bioinformatics analysis of large next generation sequencing (NGS) datasets. While the exact system specifications will depend on the size of your data, we generally use compute nodes with the Linux operating system and 16 CPU cores, 128GB of RAM, and 1TB of disk space.
* Data from a directional (stranded) whole-genome bisulfite sequencing (WGBS) run that has been properly trimmed, aligned, and quality controlled. In our testing, we used the [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) aligner on our WGBS data to create our input BAM file. For an example of a whole pipeline that can pre-process your WGBS data, please consult the [`nf-core methylseq` project](https://nf-co.re/methylseq/).
* A FASTA file of the reference genome that your BAM file was aligned to.

## Conventions

In the example below, we will call the input BAM file `sample.bam`. We will continue to label the resulting files as simply "sample". For your purposes, please replace this with the name of your actual BAM file name and sample name that you are processing. We will also call the reference genome this BAM file was aligned to `genome.fa`. Please replace `genome.fa` with the FASTA file for the reference genome you are using.

## Step 1: Software Installation

In our example, we use [Miniconda](https://docs.anaconda.com/miniconda/) and [Bioconda](https://bioconda.github.io/) to create a virtual environment where we can run the needed software for this tutorial. Feel free to install the software manually or use another dependency management solution if you would prefer.

### Option to Use a Conda Environment

Install Miniconda first, if you have not already, following the [official installation directions](https://docs.anaconda.com/miniconda/).

Creation of Conda environment.
```bash
conda create -n bioinfo
```

Add channels (bioconda needed to install `samtools`).
```bash
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

Activation of Conda environment.
```bash
conda activate bioinfo
```

Installation of dependencies.
```bash
conda install samtools
conda install -c bioconda pysam
conda install -c bioconda freebayes
```

Also, fetch the source code for Revelio from the official GitHub repository. Make sure you are in the working directory you want use for your analysis when you do this.
```bash
git clone https://github.com/bio15anu/revelio.git
```

### Option to Install Software Directly

* Refer to [Revelio's repository](https://github.com/bio15anu/revelio) for installation instructions.

* Refer to [Samtools's website](https://www.htslib.org/) for installation instructions.

* Refer to [Freebayes's website](https://github.com/freebayes/freebayes) for installation instructions.

## Step 2: Running Revelio

### Prepare BAM File

First, make sure your `genome.fa` has been indexed.
```bash
samtools faidx genome.fa
```

Then, generate MD tags based on sample BAM and reference genome FASTA. This pre-processing is necessary for Revelio to run most efficiently, as it no longer requires access to the reference genome FASTA file once MD tags have been generated.
```bash
samtools calmd -b sample.bam genome.fa 1> sample_calmd.bam 2> /dev/null
samtools index sample_calmd.bam
```

### Generated Masked BAM

Replace `<num_threads>` with the number of CPU cores available on your machine (or a lower number of your choosing). We also create a `tmp` directory in the current working directly in this example, and set Revelio to use this rather than the default, as we encountered errors with the `--temp` path being on a different device than where the output BAM file is saved.

```bash
mkdir tmp
./revelio/revelio.py --threads <num_threads> --temp tmp/ sample_calmd.bam sample_masked.bam
samtools index sample_masked.bam
```

## Step 3: Variant Calling using Freebayes

For this tutorial, we demonstrate variant calling using Freebayes, but note that you can choose another variant caller.

### Running Freebayes
Now we can call variants from the Revelio-masked BAM using the Freebayes variant caller. See the [Freebayes documentation](https://github.com/freebayes/freebayes) for additional options, as variant calling procedures often need to be fine-tuned for specific datasets and needs.
```bash
freebayes -f genome.fa sample_masked.bam > variants.vcf
```

### Quality Filtering
It is generally desirable to perform filtering on the VCF file to retain high-quality variants. Adjust this filtering based on your needs. In this case, we used QUAL > 30 and DP > 10 for quality filtering filtering.

```bash
bcftools filter -i 'QUAL > 30 && FORMAT/DP > 10' -O z -o "$output_dir/filtered1.vcf.gz" "$input_file1"`
vcffilter -f "QUAL > 20" variants.vcf > filtered_variants.vcf
```

## Notes

- Ensure all tools are properly installed and accessible in your PATH.
- Adjust file paths and parameters as necessary for your specific data and computational environment.
- The quality of the input data and the parameters used for alignment and variant calling can significantly affect the results.