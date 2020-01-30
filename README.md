# ASAP to kite
A script to process fastqs from ASAP-seq for downstream processing with kite (kallisto | bustools). 

The options are designed to mirror that of CellRanger/CellRanger-ATAC for convenience in processing

## Sample use cases

### One sample, one directory

The most basic use case is when we have one library sequenced one. From the demultiplexing,
we should see files that look like this:

```
test/data1/test1_S1_L001_R1_001.fastq.gz
test/data1/test1_S1_L001_R2_001.fastq.gz
test/data1/test1_S1_L001_R3_001.fastq.gz
test/data1/test1_S1_L002_R1_001.fastq.gz
test/data1/test1_S1_L002_R2_001.fastq.gz
test/data1/test1_S1_L002_R3_001.fastq.gz
test/data1/test1_S1_L003_R1_001.fastq.gz
test/data1/test1_S1_L003_R2_001.fastq.gz
test/data1/test1_S1_L003_R3_001.fastq.gz
test/data1/test1_S1_L004_R1_001.fastq.gz
test/data1/test1_S1_L004_R2_001.fastq.gz
test/data1/test1_S1_L004_R3_001.fastq.gz
```

Here, the sequencing run is in the folder `test/data1` and we are interested in the `test1` sample. 

We can process these fastqs:

```
python asap_to_kite_v1.py -f test/data1 -s test1 -o one_one
```

Here, the `-s` specifies the `sample name`; `-f` specifies the `fastq folder`; `-o` specifies the `output` naming convention.

### One sample, multiple directories 

If multiple sequencing rounds are performed, we can supply all sequencing libraries as a comma-separated list:

```
python asap_to_kite_v1.py -f test/data1,test/data2 -s test1 -o one_many
```

### One sample, named multiple ways, in multiple directories

Suppose that the sequencing library is named two different ways over the two sequencing runs. 
We can stack the comma-separated nature of the sample names and the sequencing runs to 
synthesize the libraries

```
python asap_to_kite_v1.py -f test/data1,test/data2 -s test1,test2 -o many_many
```

### Write to a different output destination

Finally, just to showcase that we can write these files out to a different path:

```
python asap_to_kite_v1.py -f test/data1,test/data2 -s test1,test2 -o test/many_many
```

## Important

This code works for one biological sample at a time. If multiple samples are supplied in the 
command line execution, then they will be merged (under the assumption that they were
called different things). Execute the code sequentially for each sample in the event of 
multiple biological samples. 


## Options

```
python asap_to_kite_v1.py --help
```

yields

```
Usage: asap_to_kite_v1.py [options] [inputs] Script to reformat raw sequencing 
data from CellRanger-ATAC demultiplexing to a format 
compatible with kite (kallisto|bustools)

Options:
  -h, --help            show this help message and exit
  -f FASTQS, --fastqs=FASTQS
                        Path of folder created by mkfastq or bcl2fastq; can be
                        comma separated that will be collapsed into one
                        output.
  -s SAMPLE, --sample=SAMPLE
                        Prefix of the filenames of FASTQs to select; can be
                        comma separated that will be collapsed into one output
  -o ID, --id=ID        A unique run id, used to name output.
  -c CORES, --cores=CORES
                        Number of cores for parallel processing. Default = 4.
  -n NREADS, --nreads=NREADS
                        Maximum number of reads to process in one iteration.
                        Decrease this if in a low memory environment (e.g.
                        laptop). Default = 10,000,000.
  -r, --no-rc-R2        By default, the reverse complement of R2 (barcode) is
                        performed (when sequencing with, for example, the
                        NextSeq). Throw this flag to keep R2 as is-- no
                        reverse complement (rc).
```

<br><br>

