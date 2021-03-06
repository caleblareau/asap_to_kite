import sys
import re
import shutil
import os 
import glob
import gzip

from optparse import OptionParser
from os import path
from multiprocessing import Pool, freeze_support
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from rapidfuzz import process

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Script to reformat raw sequencing \ndata from CellRanger-ATAC demultiplexing to a format \ncompatible with kite (kallisto|bustools)"
opts = OptionParser(usage=usage)

opts.add_option("--fastqs", "-f", help="Path of folder created by mkfastq or bcl2fastq; can be comma separated that will be collapsed into one output.")
opts.add_option("--TSA", "-a", help="File containing list of TotalSeq™-A barcodes")
opts.add_option("--TSB", "-b", help="File containing list of TotalSeq™-B barcodes")
opts.add_option("--sample", "-s", help="Prefix of the filenames of FASTQs to select; can be comma separated that will be collapsed into one output")
opts.add_option("--id", "-o", default = "asap2kite", help="A unique run id, used to name output.")

opts.add_option("--cores", '-c', default = 4, help="Number of cores for parallel processing. Default = 4.")
opts.add_option("--nmismatch", '-m', default = 2, help="Maximum number of mismatches to tolerate when looking for TSA/TSB assignments.")
opts.add_option("--nreads", '-n', default = 10000000, help="Maximum number of reads to process in one iteration. Decrease this if in a low memory environment (e.g. laptop). Default = 10,000,000.")

opts.add_option("--no-rc-R2", '-r', action="store_true", default = False, help="By default, the reverse complement of R2 (barcode) is performed (when sequencing with, for example, the NextSeq). Throw this flag to keep R2 as is-- no reverse complement (rc).")
options, arguments = opts.parse_args()

folder_of_fastqs = options.fastqs
sample_name = options.sample
out = options.id
n_cpu = int(options.cores)
n_reads= int(options.nreads)
n_mismatch= int(options.nmismatch)

rc_R2 = not (options.no_rc_R2)

# import barcodes
A_barcodes = [line.rstrip('\n') for line in open(options.TSA)]
B_barcodes = [line.rstrip('\n') for line in open(options.TSB)]

# Reformat read for export
def formatRead(title, sequence, quality):
	return("@%s\n%s\n+\n%s\n" % (title, sequence, quality))


def asap_to_kite_vDetect(trio):
	listRead1 = trio[0]; listRead2 = trio[1]; listRead3 = trio[2]
	
	# Parse aspects of existing read
	title1 = listRead1[0]; sequence1 = listRead1[1]; quality1 = listRead1[2]
	title2 = listRead2[0]; sequence2 = listRead2[1]; quality2 = listRead2[2]
	title3 = listRead3[0]; sequence3 = listRead3[1]; quality3 = listRead3[2]
	
	# process R2
	if(rc_R2): # rc_R2
		# Update sequence with reverse complement
		bio_r2 = str(Seq(sequence2).reverse_complement())
		sequence2 = bio_r2
		
		# update the quality
		quality2 = quality2[::-1]
	
	conjugation = ""
	maybeA = sequence3[0:14]
	maybeB = sequence3[10:24]
	A_match = process.extractOne(maybeA, A_barcodes)
	B_match = process.extractOne(maybeB, B_barcodes)

	if(A_match[1] >= (len(maybeA)-n_mismatch)/len(maybeA)*100): 
		conjugation = "TotalSeqA"

	if(B_match[1] >= (len(maybeB)-n_mismatch)/len(maybeB)*100): 
		conjugation = "TotalSeqB"

	new_sequence1 = ""; new_sequence2 = ""
	new_quality1 = ""; new_quality2 = ""
	
	# Recombine attributes based on conjugation logic
	if(conjugation == "TotalSeqA"):
		new_sequence1 = sequence2 + sequence1[0:10]
		new_sequence2 = sequence3
	
		new_quality1 = quality2 + quality1[0:10]
		new_quality2 = quality3
		
	elif(conjugation == "TotalSeqB"):
		new_sequence1 = sequence2 + sequence3[0:10] + sequence3[25:34]
		new_sequence2 = sequence3[10:25]
	
		new_quality1 = quality2 + quality3[0:10] + quality3[25:34]
		new_quality2 = quality3[10:25]
	
	# Prepare reads for exporting
	out_fq1 = formatRead(title1, new_sequence1, new_quality1)
	out_fq2 = formatRead(title2, new_sequence2, new_quality2)
	if(conjugation == "TotalSeqA"):
		return([[out_fq1, out_fq2, "", ""]])
	elif(conjugation == "TotalSeqB"):
		return([["", "", out_fq1, out_fq2]])
	else:
		return([['','','','']])

if __name__ == '__main__':
	
	print("\nASAP-to-kite, version A-B\nEnables simultaneous parsing and detection of TotalSeqA and B tags in one sample.\n")

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

	# Process through iterator
	def batch_iterator(iterator, batch_size):
		entry = True  # Make sure we loop once
		while entry:
			batch = []
			while len(batch) < batch_size:
				try:
					entry = iterator.__next__()
				except StopIteration:
					entry = None
				if entry is None:
					# End of file
					break
				batch.append(entry)
			if batch:
				yield batch



	# Main loop -- process input reads and write out the processed fastq files
	print("\nProcessing these fastq samples: ")
	for r in R1s_for_analysis:
		print(r.replace("_R1_001.fastq.gz", ""))

	outfq1fileA = out + "_TSA_R1.fastq.gz"
	outfq2fileA = out + "_TSA_R2.fastq.gz"
	outfq1fileB = out + "_TSB_R1.fastq.gz"
	outfq2fileB = out + "_TSB_R2.fastq.gz"

	with gzip.open(outfq1fileA, "wt") as out_Af1:
		with gzip.open(outfq2fileA, "wt") as out_Af2:
			with gzip.open(outfq1fileB, "wt") as out_Bf1:
				with gzip.open(outfq2fileB, "wt") as out_Bf2:
					for R1file in R1s_for_analysis:
						R2file = R1file.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
						R3file = R1file.replace("_R1_001.fastq.gz", "_R3_001.fastq.gz")
		
						# Read in fastq in chunks the size of the maximum user tolerated number
						it1 = batch_iterator(FastqGeneralIterator(gzip.open(R1file, "rt")), n_reads)
						it2 = batch_iterator(FastqGeneralIterator(gzip.open(R2file, "rt")), n_reads)
						it3 = batch_iterator(FastqGeneralIterator(gzip.open(R3file, "rt")), n_reads)
		
						for i, batch_R1 in enumerate(it1):
							batch_R2 = it2.__next__()
							batch_R3 = it3.__next__()

							pool = Pool(processes=n_cpu)
							pm = pool.map(asap_to_kite_vDetect, zip(batch_R1, batch_R2, batch_R3))
							pool.close()
				
							# process and write out
							fq_data = list(map(''.join, zip(*[item.pop(0) for item in pm])))
							out_Af1.writelines(fq_data[0])
							out_Af2.writelines(fq_data[1])
							out_Bf1.writelines(fq_data[2])
							out_Bf2.writelines(fq_data[3])
						
			print("\nDone!\n")
		
		
