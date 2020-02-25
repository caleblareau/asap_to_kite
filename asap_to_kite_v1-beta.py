import sys
import re
import shutil
import os 
import glob
import gzip
import string

from optparse import OptionParser
from os import path


opts = OptionParser()
usage = "usage: %prog [options] [inputs] Script to reformat raw sequencing \ndata from CellRanger-ATAC demultiplexing to a format \ncompatible with kite (kallisto|bustools)"
opts = OptionParser(usage=usage)

opts.add_option("--fastqs", "-f", help="Path of folder created by mkfastq or bcl2fastq; can be comma separated that will be collapsed into one output.")
opts.add_option("--sample", "-s", help="Prefix of the filenames of FASTQs to select; can be comma separated that will be collapsed into one output")
opts.add_option("--id", "-o", default = "asap2kite", help="A unique run id, used to name output.")
opts.add_option("--cores", '-c', default = 4, help="Number of cores for parallel processing. Default = 4.")
opts.add_option("--nreads", '-n', default = 10000000, help="Maximum number of reads to process in one iteration. Decrease this if in a low memory environment (e.g. laptop). Default = 10,000,000.")
opts.add_option("--no-rc-R2", '-r', action="store_true", default = False, help="By default, the reverse complement of R2 (barcode) is performed (when sequencing with, for example, the NextSeq). Throw this flag to keep R2 as is-- no reverse complement (rc).")
options, arguments = opts.parse_args()

folder_of_fastqs = options.fastqs
sample_name = options.sample
out = options.id
n_cpu = int(options.cores)
n_reads= int(options.nreads)
rc_R2 = not (options.no_rc_R2)
print("\nASAP-to-kite, version 1\n")

print("User specified options: ")
print(options)

#----------
# Functions - I/O
#----------

# Function to verify R1/R2/R3 are present for nominated samples
def verify_sample_from_R1(list_of_R1s):
	verified_R1s = []
	for R1file in list_of_R1s:
		R2file = R1file.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
		R3file = R1file.replace("_R1_001.fastq.gz", "_R3_001.fastq.gz")
		if(path.exists(R2file) and path.exists(R3file)):
			verified_R1s.append(R1file)
	return(verified_R1s)

# identify all sequencing data that should be parsed for conversion
def parse_directories(folder_of_fastqs):
	list_folders = folder_of_fastqs.split(",")
	list_samples = sample_name.split(",")

	all_R1s = []
	
	# Look into all supplied folders for specific files:
	for path_to_folder in list_folders:
		
		# Look at all of the possible sample names
		for sample_name_one in list_samples:
			matching_R1s = glob.glob(path_to_folder+"/*" + sample_name_one + "*" + "_R1_001.fastq.gz")
			for file in matching_R1s:
				all_R1s.append(file)
	verified_R1s = verify_sample_from_R1(all_R1s)
	return(verified_R1s)

# Import files
R1s_for_analysis = parse_directories(folder_of_fastqs)

tab = str.maketrans("ACTG", "TGAC")
def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]

# Fast fastq generator function
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break
                
# Reformat read for export
def formatRead(title, sequence, quality):
	return("@%s\n%s\n+\n%s\n" % (title, sequence, quality))

def asap_to_kite_v2(title1, title2, title3, sequence1, sequence2, sequence3, quality1, quality2, quality3):	

	# process R2
	if(rc_R2):
	
		# Update sequence with reverse complement
		sequence2 = reverse_complement_table(sequence2)
		
		# update the quality
		quality2 = quality2[::-1]
	
	# Recombine attributes
	new_sequence1 = sequence2 + sequence1[0:10]
	new_sequence2 = sequence3
	
	new_quality1 = quality2 + quality1[0:10]
	new_quality2 = quality3
	
	out_fq1 = formatRead(title1, new_sequence1, new_quality1)
	out_fq2 = formatRead(title2, new_sequence2, new_quality2)
	return([out_fq1, out_fq2])

# Main loop -- process input reads and write out the processed fastq files
print("\nProcessing these fastq samples: ")
for r in R1s_for_analysis:
	print(r.replace("_R1_001.fastq.gz", ""))

outfq1file = out + "_R1.fastq.gz"
outfq2file = out + "_R2.fastq.gz"
with gzip.open(outfq1file, "wt") as out_f1:
	with gzip.open(outfq2file, "wt") as out_f2:
		for R1file in R1s_for_analysis:
			R2file = R1file.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
			R3file = R1file.replace("_R1_001.fastq.gz", "_R3_001.fastq.gz")
		
			# Read in fastq in chunks the size of the maximum user tolerated number
			fh1 = gzip.open(R1file, "rt")
			fh2 = gzip.open(R2file, "rt")
			fh3 = gzip.open(R3file, "rt")
			for name1, seq1, qual1 in readfq(fh1):				
				name2, seq2, qual2 = readfq(fh2).__next__()
				name3, seq3, qual3 = readfq(fh3).__next__()
				
				# process and write out
				fq_data = asap_to_kite_v2(name1, name2, name3, seq1, seq2, seq3, qual1, qual2, qual3)
				out_f1.writelines(fq_data[0])
				out_f2.writelines(fq_data[1])
				
print("\nDone!\n")
