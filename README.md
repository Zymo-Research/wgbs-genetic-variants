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
conda install -c bioconda freebayes=1.3.6
```

Note we have requested a specific version of `freebayes` as we encountered runtime errors with version 1.3.7 at time of writing.

Also, fetch the source code for Revelio from the official GitHub repository. Make sure you are in the working directory you want use for your analysis when you do this.
```bash
git clone https://github.com/bio15anu/revelio.git
```

### Option to Install Software Directly

* Refer to [Revelio's repository](https://github.com/bio15anu/revelio) for installation instructions.

* Refer to [Samtools's website](https://www.htslib.org/) for installation instructions.

* Refer to [Freebayes's website](https://github.com/freebayes/freebayes) for installation instructions.

### Set Number of CPU Cores

Various tools will require an argument for the number of CPU cores you wish to use. We will use an environment variable to pass this value in our example. If you wish to do the same, please set this variable, replacing `<num_cores>` with the number of CPU cores you wish to use.

```bash
export NUM_CORES=<num_cores>
```

## Step 2: Running Revelio

### Prepare BAM File

First, make sure your `genome.fa` has been indexed.
```bash
samtools faidx genome.fa
```

Then, generate MD tags based on sample BAM and reference genome FASTA. This pre-processing is necessary for Revelio to run most efficiently, as it no longer requires access to the reference genome FASTA file once MD tags have been generated.
```bash
samtools calmd -@ $NUM_CORES -b sample.bam genome.fa 1> sample_calmd.bam 2> /dev/null
samtools index -@ $NUM_CORES sample_calmd.bam
```

### Generated Masked BAM

Replace `<num_threads>` with the number of CPU cores available on your machine (or a lower number of your choosing). We also create a `tmp` directory in the current working directly in this example, and set Revelio to use this rather than the default, as we encountered errors with the `--temp` path being on a different device than where the output BAM file is saved.

```bash
mkdir tmp
./revelio/revelio.py --threads $NUM_CORES --temp tmp/ sample_calmd.bam sample_masked.bam
samtools index -@ $NUM_CORES sample_masked.bam
```

### Sort the BAM file

Variant callers (including Freebayes) often require BAMs that have been coordinate sorted. To ensure our BAM is properly sorted, we will perform sorting first. Again, replace <num_threads> with the number of CPU cores you wish to use.

```bash
samtools sort -@ $NUM_CORES -o sample_masked_sorted.bam sample.bam
```

## Step 3: Variant Calling using Freebayes

For this tutorial, we demonstrate variant calling using Freebayes, but note that you can choose another variant caller.

### Running Freebayes
Now we can call variants from the Revelio-masked BAM using the Freebayes variant caller. See the [Freebayes documentation](https://github.com/freebayes/freebayes) for additional options, as variant calling procedures often need to be fine-tuned for specific datasets and needs. In this case, we specifically compress and index the VCF file as well to save disk space and be compatible with downstream tools.
```bash
freebayes -f genome.fa sample_masked_sorted.bam | bgzip -c > variants.vcf.gz
tabix -@ $NUM_CORES -p vcf variants.vcf.gz
```

### Quality Filtering
It is generally desirable to perform filtering on the VCF file to retain high-quality variants. Adjust this filtering based on your needs. In this case, we used QUAL > 30 (quality score from the variant caller, PHRED-scaled) and DP > 10 (sequencing depth at the position, in reads) for quality filtering filtering. Please see the [`bcftools` documentation](https://samtools.github.io/bcftools/bcftools.html#filter) for more information about filtering options. What we have provided here is just an illustrative example. Note there is a difference between "hard filtering" (removing variants from the VCF itself) and "soft filtering" (flagging but not removing variants that do not pass the filter). For simplicity, we've show hard filtering here.

```bash
bcftools filter -i 'QUAL>30 && FORMAT/DP>10' -Oz -o variants_filtered.vcf.gz variants.vcf.gz
bcftools index --threads $NUM_CORES variants_filtered.vcf.gz

```

## Notes

- Ensure all tools are properly installed and accessible in your PATH.
- Adjust file paths and parameters as necessary for your specific data and computational environment.
- The quality of the input data and the parameters used for alignment and variant calling can significantly affect the results.

## Citations

### Tools Used or Mentioned in This How To

**VCF** file format used as the output:

* Danecek, Petr, Adam Auton, Goncalo Abecasis, Cornelis A. Albers, Eric Banks, Mark A. DePristo, Robert E. Handsaker, et al. 2011. “The Variant Call Format and VCFtools.” Bioinformatics (Oxford, England) 27 (15): 2156–58. [[link]](https://doi.org/10.1093/bioinformatics/btr330)

**nf-core** Project recommended for a pre-processing pipeline for WGBS data:

* Ewels, Philip A., Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso, and Sven Nahnsen. 2020. “The Nf-Core Framework for Community-Curated Bioinformatics Pipelines.” Nature Biotechnology 38 (3): 276–78. [[link]](https://doi.org/10.1038/s41587-020-0439-x)

**FreeBayes** used for variant calling on the Revelio-processed BAM file:

* Garrison, Erik, and Gabor Marth. 2012. “Haplotype-Based Variant Detection from Short-Read Sequencing.” arXiv. [[link]](https://doi.org/10.48550/arXiv.1207.3907)

**Bioconda** Project recommended for dependency management:

* Grüning, Björn, Ryan Dale, Andreas Sjödin, Brad A. Chapman, Jillian Rowe, Christopher H. Tomkins-Tinch, Renan Valieris, Johannes Köster, and Bioconda Team. 2018. “Bioconda: Sustainable and Comprehensive Software Distribution for the Life Sciences.” Nature Methods 15 (7): 475–76. [[link]](https://doi.org/10.1038/s41592-018-0046-7)

**bismark** Aligner recommended for WGBS alignment:

* Krueger, Felix, and Simon R. Andrews. 2011. “Bismark: A Flexible Aligner and Methylation Caller for Bisulfite-Seq Applications.” Bioinformatics (Oxford, England) 27 (11): 1571–72. [[link]](https://doi.org/10.1093/bioinformatics/btr167)

**samtools** used for general BAM file processing:

* Li, Heng, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, and 1000 Genome Project Data Processing Subgroup. 2009. “The Sequence Alignment/Map Format and SAMtools.” Bioinformatics (Oxford, England) 25 (16): 2078–79. [[link]](https://doi.org/10.1093/bioinformatics/btp352)

**Revelio** used to process WGBS BAM file to be compatible with variant calling software:

* Nunn, Adam, Christian Otto, Mario Fasold, Peter F. Stadler, and David Langenberger. 2022. “Manipulating Base Quality Scores Enables Variant Calling from Bisulfite Sequencing Alignments Using Conventional Bayesian Approaches.” BMC Genomics 23 (1): 477. [[link]](https://doi.org/10.1186/s12864-022-08691-6)

### Related Papers and Alternative Tools

This list includes other tools that can extract variants from WGBS data and papers comparing among such tools. These tools could be considered as alternatives to the Revelio+FreeBayes approach used in this How To.

**MethylExtract** contains a function for WGBS variant calling:

* Barturen, Guillermo, Antonio Rueda, José L. Oliver, and Michael Hackenberg. 2013. “MethylExtract: High-Quality Methylation Maps and SNV Calling from Whole Genome Bisulfite Sequencing Data.” F1000Research 2:217. [[link]](https://doi.org/10.12688/f1000research.2-217.v2)

**BS-SNPer** tool for WGBS variant calling:

* Gao, Shengjie, Dan Zou, Likai Mao, Huayu Liu, Pengfei Song, Youguo Chen, Shancen Zhao, et al. 2015. “BS-SNPer: SNP Calling in Bisulfite-Seq Data.” Bioinformatics (Oxford, England) 31 (24): 4006–8. [[link]](https://doi.org/10.1093/bioinformatics/btv507)

**CGmapTools** contains a function for WGBS variant calling:

* Guo, Weilong, Ping Zhu, Matteo Pellegrini, Michael Q. Zhang, Xiangfeng Wang, and Zhongfu Ni. 2018. “CGmapTools Improves the Precision of Heterozygous SNV Calls and Supports Allele-Specific Methylation Detection and Visualization in Bisulfite-Sequencing Data.” Bioinformatics (Oxford, England) 34 (3): 381–87. [[link]](https://doi.org/10.1093/bioinformatics/btx595)

**MethGo** contains a function for WGBS variant calling:

* Liao, Wen-Wei, Ming-Ren Yen, Evaline Ju, Fei-Man Hsu, Larry Lam, and Pao-Yang Chen. 2015. “MethGo: A Comprehensive Tool for Analyzing Whole-Genome Bisulfite Sequencing Data.” BMC Genomics 16 Suppl 12 (Suppl 12): S11. [[link]](https://doi.org/10.1186/1471-2164-16-S12-S11)

Comparison study of various methods for WGBS variant calling:

* Lindner, Melanie, Fleur Gawehns, Sebastiaan Te Molder, Marcel E. Visser, Kees van Oers, and Veronika N. Laine. 2022. “Performance of Methods to Detect Genetic Variants from Bisulphite Sequencing Data in a Non-Model Species.” Molecular Ecology Resources 22 (2): 834–46. [[link]](https://doi.org/10.1111/1755-0998.13493)

**Bis-SNP** tool for WGBS variant calling:

* Liu, Yaping, Kimberly D. Siegmund, Peter W. Laird, and Benjamin P. Berman. 2012. “Bis-SNP: Combined DNA Methylation and SNP Calling for Bisulfite-Seq Data.” Genome Biology 13 (7): R61. [[link]](https://doi.org/10.1186/gb-2012-13-7-r61)

**gemBS** contains a function for WGBS variant calling:

* Merkel, Angelika, Marcos Fernández-Callejo, Eloi Casals, Santiago Marco-Sola, Ronald Schuyler, Ivo G. Gut, and Simon C. Heath. 2019. “gemBS: High Throughput Processing for DNA Methylation Data from Bisulfite Sequencing.” Bioinformatics (Oxford, England) 35 (5): 737–42. [[link]](https://doi.org/10.1093/bioinformatics/bty690)

The EpiDiverse Toolkit implements a larger pipeline around the Revelio tool:

* Nunn, Adam, Sultan Nilay Can, Christian Otto, Mario Fasold, Bárbara Díez Rodríguez, Noé Fernández-Pozo, Stefan A. Rensing, Peter F. Stadler, and David Langenberger. 2021. “EpiDiverse Toolkit: A Pipeline Suite for the Analysis of Bisulfite Sequencing Data in Ecological Plant Epigenetics.” NAR Genomics and Bioinformatics 3 (4): lqab106. [[link]](https://doi.org/10.1093/nargab/lqab106)

**BISCUIT** contains a function for WGBS variant calling:

* Zhou, Wanding, Benjamin K. Johnson, Jacob Morrison, Ian Beddows, James Eapen, Efrat Katsman, Ayush Semwal, et al. 2024. “BISCUIT: An Efficient, Standards-Compliant Tool Suite for Simultaneous Genetic and Epigenetic Inference in Bulk and Single-Cell Studies.” Nucleic Acids Research 52 (6): e32. [[link]](https://doi.org/10.1093/nar/gkae097)
