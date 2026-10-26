# RNAChrom_ALLvsALL_data_processing

Package for processing data of several "ALL-to-ALL" RNA-Chromatin interactions capturing experiments, e.g. [RADICL](https://pubmed.ncbi.nlm.nih.gov/32094342/), [GRID](https://pubmed.ncbi.nlm.nih.gov/31175345/), [iMARGI](https://pubmed.ncbi.nlm.nih.gov/30718424/), [CHAR](https://pubmed.ncbi.nlm.nih.gov/29648534/), etc.

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
    1. [Required software](#software1)
    2. [Setting up virtual environment](#venv)
    3. [Actual installation](#install1)
    4. [Data preparation](#genomeprep)
3. [List of rnachromprocessing tools](#tools)
    1. [rnachromprocessing](#rnachromprocessing)
    2. [detect-strand](#detect-strand)
    3. [X-RNA inference](#xrna)

<a name="overview"></a>
##  Overview
This package is desined to unify the data processing of different RNA-Chromatin interactions capturing experiments.

After installation it provides access to several command-line tools:
* rnachromprocessing -- pipeline for uniform processing of "All-versus-All" RNA-Chromatin interactions data. 
* detect-strand -- program that checks whether orientation of RNA parts of contacts was inverted during sequencing. 
* infer-xrna -- pipeline that searches for most promising not-annotated transcripts in the experimental data.

<a name="software1"></a>
## Installation
### Required software

Basic dependencies:

* Python (v 3.12 or higher)
* Common unix command-line text-processing utilities: grep, sed, awk, cut, uniq etc.

This package also requires several bioinformatic tools to run.

* [samtools](http://www.htslib.org/) ver 1.11 or higher (should be available via PATH variable)
* [bedtools](https://bedtools.readthedocs.io/en/latest/) (should be available via PATH variable)
* [hisat2](http://daehwankimlab.github.io/hisat2/)
* any of the supported trimming tools ([trimmomatic](https://academic.oup.com/bioinformatics/article/30/15/2114/2390096), ...)
* any of the supported deduplication tools ([fastuniq](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0052249), 
                                            [fastq-dupaway](https://github.com/AndrewSigorskih/fastq-dupaway), ...)

<a name="venv"></a>
### Setting up virtual environment (recommended)
It is highly recommended to install package using in a virtual environment. This can be achieved via [python's venv](https://docs.python.org/3/library/venv.html) or [conda envs](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).<br> Moreover, all of the required bioinformatic tools can be installed via conda as well (but mind the tricky [conda samtools installation](https://www.biostars.org/p/455593/)!) 


<a name="install1"></a>
### Clone and install package
Clone repository from github and install. All non-python dependencies are listed in the <b>environment.yml</b> file. We highly recommend using [mamba](https://anaconda.org/conda-forge/mamba) for faster environment setup, however installation via vanilla conda is also possible:

```bash
git clone https://github.com/AndrewSigorskih/RNAChrom_ALLvsALL_data_processing.git
cd RNAChrom_ALLvsALL_data_processing
# set up virtual environment with all bioinformatics tools
# (provided you use conda as your package manager and have already 
# installed mamba and/or conda-libmamba-solver in base environment)
mamba env create -f environment/py312_conda_env.yml
conda activate rnachromproc_py12
# install package with all its python dependincies
pip install .
```

<a name="genomeprep"></a>
### Data preparation

#### Paired-end reads order
The *rnachromprocessing* pipeline expects reads in paired files to be synchronized, i.e. for every index i read number i from rna.fastq file and read number i from dna.fastq file have same id. If you download data from SRA archive using fasterq-dump make sure to use `--split-3` flag.

#### Filtering genome fasta
Chromosome names in genome and gene annotation should be the same. We also recommend to remove mitochondrial and "Unplaced/Unlocalized" sequences from genome fasta. This can be easily achieved using tools like [seqkit](https://bioinf.shenwei.me/seqkit/): for example, command `seqkit grep -vrp "^chrUn" file.fa > clean.fa` will remove all such chromosome fragments from genome file.

#### Genome indexing
Hisat2 requires genome to be indexed. In order to build the index, use command `hisat2-build -p 16 genome.fa prefix`, where **prefix** is the prefix name for genome index files. In order to properly map the RNA parts you will also need the splicecite file, that can be obtained from gene annotation file using hisat's script: `hisat2_extract_splice_sites.py genome.gtf > genome.ss`. See [hisat2 manual](http://daehwankimlab.github.io/hisat2/howto/#building-indexes) for more information.

<a name="tools"></a>
## List of rnachromprocessing tools

### [rnachromprocessing](docs/rnachromprocessing/README.md)

This program runs the pipeline of RNA-Chromatin interactions data processing. It takes raw fastq files as input and produces tables of RNA-DNA contacts as main output.

---

### [detect-strand](docs/detect-strand/README.md)

This program checks whether orientation of RNA parts of contacts was inverted or lost during sequencing.

---

<a name="xrna"></a>
### [infer-xrna](docs/x-rna/README.md)

This program tries to infer novel transcripts that do not correspond to any gene from user-provided annotation.
