<a name="head"></a>
# infer-xrna

This program tries to infer novel transcripts that do not correspond to any gene from user-provided annotation. It heavily relies on the results from 2 previous tools: *rnachromprocessing* and *detect-strand*.

## Table of Contents
1. [Inputs](#inputs)
    * [Requirements](#requirements)
2. [Usage](#usage)
3. [Config contents](#config)
4. [Outputs](#outputs)
5. [Procedure overview](#procedure)

<a name="inputs"></a>
## Inputs

* List of ids corresponding to RNA bam files to be processed.
* A path to a directory containing *bam* files, preferably obtained by using the [rnachromprocessing](../rnachromprocessing/README.md) tool from this package. Ids of the files should match with those listed as input ids.
* Reference gene annotation in **GTF** format. Reads that overlap with *any* interval from this annotation will be filtered out and not used for transcript inference.
* Main output of the [detect-strand](../detect-strand/README.md) tool from this package. This file should contain information about files to be processed.

**!NB!** Your reference gene annotation **must** contain both "gene" and "transcript" records corresponding to the genomic intervals that have a known / proven product. The "transcript" records are used by Stringtie tool directly, while the "gene" records will be used by this tool to determine which reads to exclude from analysis.

<a name="requirements"></a>
#### Requirements

Anongside the Python requirements used by the package, this module requires the following programs to be present:

* bedtools (buit and tested with v2.31)
* samtools (buit and tested with v1.24)
* Stringtie (buit and tested with v3.0.3)
* awk


<a name="usage"></a>
## Usage

```bash
infer-xrna [-h] -c CONFIG [-v/-vv]
```

Supported arguments:
```
-h, --help          show help message and exit
-c CONFIG, --config CONFIG
                    Configuration file in json or yaml format (required).
-v, --verbose       Verbosity level. By default little to none information is printed.
                    Use -v once to increase information logs about each step, and -vv to 
                    print every command that is being run.
```

<a name="config"></a>
## Config contents

Field name|Required|Description
---|---|---
file_ids|Yes|non-empty array of strings. **Must** contain base names of files with RNA parts of reads, without directory names or file extensions. For example, for the file named <i>SRR9201799_1.fastq.gz</i> id is SRR9201799_1.
bam_input_dir|Yes|A path to a directory containing bam input files.
output_dir|Yes|Path to directory to store results in. Will be created if doesn't exist. Resulting output files and folders will be copied into this path.
base_dir|No|string. Directory to run analysis in (a temporal directory will be created in this directory). Path should exist if provided. If not provided, current working directory will be used instead.
cpus|No|int > 0. Number of tasks to run simultaneously on each step. Default is 1.
ouputs_prefix|No|string. Prefix for main output files and plots. Default value is "xrna".
keep_extras|No|boolean. Whether or not save all intermediate files alongside the main output. Mianly used for debug.
annotation|Yes|A sub-config containing paths to input annotation files (see below).
stringtie|No|A sub-config containing information required to run stringtie tool (see below).

#### annotation sub-config:

Field name|Required|Description
---|---|---
gene_annotation|Yes|Path to reference gene annotation in **GTF** format.
strand_info|Yes|Path to main output of the `detect-strand` tool.


#### stringtie sub-config:
Field name|Required|Description
---|---|---
tool_path|No|Path to stringtie executable if not available via PATH.
cpus|No|int > 0. Number of tasks to run simultaneously for all parallel stringtie calls. If not set, the "global" value will be used.
stringtie_threads|No|int > 0. How many threads stringtie will use **for every** simultaneously-run task. (Value to pass to stringtie's `-p` option). Default is 1.

[Back to top](#head)
