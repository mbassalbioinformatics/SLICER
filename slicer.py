#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################

####     READ CONFIGURATION 1     ####     BBS-CoreSeq-LFS-Barcode-RFS     ####

###    BackBone  |  BBS end  |        Core sequence        | LFS start |      LFS      |  LFS end  | Barcode | RFS start |           RFS           |  BackBone   
### <------------|===========|-----------------------------|===========|---------------|===========|---------|===========|-------------------------|------------>

### On READ CONFIGURATION 1, this script searches for === sequences to locate barcode and core sequence


###    BackBone  |  BBS end  |        Core sequence        | LFS start |      LFS      |  LFS end  | Barcode | RFS start |           RFS           |  BackBone   
### <------------|-----------|+++++++++++++++++++++++++++++|-----------|---------------|+++++++++++|+++++++++|+++++++++++|-------------------------|------------>

### On READ CONFIGURATION 1, this script extracts +++ sequences for alignment and discards the rest of the read

###############################################################################

####     READ CONFIGURATION 2     ####     LFS-Barcode-RFS-CoreSeq-BBS     ####

###    BackBone  |           LFS           |  LFS end  | Barcode | RFS start |      RFS      |  RFS end  |        Core sequence        | BBS start |  BackBone   
### <------------|-------------------------|===========|---------|===========|---------------|===========|-----------------------------|===========|------------>

### On READ CONFIGURATION 2, this script searches for === sequences to locate barcode and core sequence


###    BackBone  |           LFS           |  LFS end  | Barcode | RFS start |      RFS      |  RFS end  |        Core sequence        | BBS start |  BackBone   
### <------------|-------------------------|+++++++++++|+++++++++|+++++++++++|---------------|-----------|+++++++++++++++++++++++++++++|-----------|------------>

### On READ CONFIGURATION 2, this script extracts +++ sequences for alignment and discards the rest of the read

###############################################################################

from genericpath import isfile
import os
import pysam # required for plotting
import seaborn as sns # required for plotting
from matplotlib import pyplot as plt # required for plotting
import argparse
import pyfastx
from os import listdir
import sys
import numpy as np
import pandas as pd
import multiprocess
import threading
import math
import time
import csv
from venn import venn
from scipy import stats
from scipy.spatial.distance import hamming
import Levenshtein as lev
from colorama import Fore, Back, Style

###############################################################################

def chunks(lst, n):
    # Yield successive n-sized chunks from lst
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def fastq_filtering_function(fastq_filtering_function_arg, fastq_filtering_queue):

	read_big_chunk = fastq_filtering_function_arg[0]
	process_number = fastq_filtering_function_arg[1]

	start_index = time.time()
	raw_fq = pyfastx.Fastq(step_01_output_filename)
	end_index = time.time()
	# print('Child process #{} FASTQ object successfully generated in {} seconds'.format(process_number + 1, (end_index - start_index)))

	number_of_threads = 10
	read_small_chunk_list = list(chunks(read_big_chunk, number_of_threads))

	for read_small_chunk in read_small_chunk_list:
		
		current_read_name_list = [('@{}'.format(raw_fq[read_num].name)) for read_num in read_small_chunk]
		current_read_seq_list = [(raw_fq[read_num].seq) for read_num in read_small_chunk]
		current_read_qual_list = [(raw_fq[read_num].qual) for read_num in read_small_chunk]
		current_read_record_list = list(zip(current_read_name_list, current_read_seq_list, current_read_qual_list))

		if sys.platform == 'darwin': # Single thread if running on MacOS
			for current_read_record in current_read_record_list:
				fastq_filtering_thread(current_read_record, fastq_filtering_queue)

		if sys.platform == 'linux' or sys.platform == 'linux2': # Multithread if running on Linux
			# define thread
			try:
				t1 = threading.Thread(target = fastq_filtering_thread, args = (current_read_record_list[0], fastq_filtering_queue))
				t2 = threading.Thread(target = fastq_filtering_thread, args = (current_read_record_list[1], fastq_filtering_queue))
				t3 = threading.Thread(target = fastq_filtering_thread, args = (current_read_record_list[2], fastq_filtering_queue))
				t4 = threading.Thread(target = fastq_filtering_thread, args = (current_read_record_list[3], fastq_filtering_queue))
				t5 = threading.Thread(target = fastq_filtering_thread, args = (current_read_record_list[4], fastq_filtering_queue))
				t6 = threading.Thread(target = fastq_filtering_thread, args = (current_read_record_list[5], fastq_filtering_queue))
				t7 = threading.Thread(target = fastq_filtering_thread, args = (current_read_record_list[6], fastq_filtering_queue))
				t8 = threading.Thread(target = fastq_filtering_thread, args = (current_read_record_list[7], fastq_filtering_queue))
				t9 = threading.Thread(target = fastq_filtering_thread, args = (current_read_record_list[8], fastq_filtering_queue))
				t10 = threading.Thread(target = fastq_filtering_thread, args = (current_read_record_list[9], fastq_filtering_queue))
			
			except:
				pass
				
			# starting thread
			try:
				t1.start()
				t2.start()
				t3.start()
				t4.start()
				t5.start()
				t6.start()
				t7.start()
				t8.start()
				t9.start()
				t10.start()
			
			except:
				pass

			# wait until all threads completely executed
			try:
				t1.join()
				t2.join()
				t3.join()
				t4.join()
				t5.join()
				t6.join()
				t7.join()
				t8.join()
				t9.join()
				t10.join()
			
			except:
				pass


def fastq_filtering_thread(fastq_filtering_thread_arg, fastq_filtering_queue):

	current_read_name = fastq_filtering_thread_arg[0]
	current_read_seq = fastq_filtering_thread_arg[1]
	current_read_qual = fastq_filtering_thread_arg[2]

	try:
		marker_1_index = current_read_seq.index(marker_1_seq)
		marker_1_error_state = 0
	except:
		marker_1_error_state = 1

	try:
		marker_2_index = current_read_seq.index(marker_2_seq)
		marker_2_error_state = 0
	except:
		marker_2_error_state = 1

	try:
		marker_3_index = current_read_seq.index(marker_3_seq)
		marker_3_error_state = 0
	except:
		marker_3_error_state = 1

	try:
		marker_4_index = current_read_seq.index(marker_4_seq)
		marker_4_error_state = 0
	except:
		marker_4_error_state = 1


	if marker_1_error_state == 0 and marker_2_error_state == 0 and marker_3_error_state == 0 and marker_4_error_state == 0:

		current_barcode = current_read_seq[(marker_1_index + search_length) : marker_2_index]

		current_coreseq = current_read_seq[(marker_3_index + search_length) : marker_4_index]

		if read_config == 1:
			
			current_trimmed_read_seq = current_read_seq[(marker_3_index + search_length) : marker_4_index] + \
									current_read_seq[marker_1_index : (marker_2_index + search_length)]

			current_trimmed_read_qual = current_read_qual[(marker_3_index + search_length) : marker_4_index] + \
									current_read_qual[marker_1_index : (marker_2_index + search_length)]
			
		if read_config == 2:
			
			current_trimmed_read_seq = current_read_seq[marker_1_index : (marker_2_index + search_length)] + \
									current_read_seq[(marker_3_index + search_length) : marker_4_index]

			current_trimmed_read_qual = current_read_qual[marker_1_index : (marker_2_index + search_length)] + \
									current_read_qual[(marker_3_index + search_length) : marker_4_index]
	
		fastq_filtering_queue.put([
			current_read_name,
			current_trimmed_read_seq,
			current_trimmed_read_qual,
			marker_1_error_state,
			marker_2_error_state,
			marker_3_error_state,
			marker_4_error_state,
			current_barcode, 
			current_coreseq])

	else:
		fastq_filtering_queue.put([
			current_read_name,
			'NA',
			'NA',
			marker_1_error_state,
			marker_2_error_state,
			marker_3_error_state,
			marker_4_error_state,
			'NA',
			'NA'])


def fastq_filtering_listener(fastq_filtering_queue):
	'''listens for messages on the q, writes to file. '''

	if read_config == 1:
		discarded_read_dict = {
			'LFS end missing': set([]),
			'RFS start missing': set([]),
			'BBS end missing': set([]),
			'LFS start missing': set([])
			}

	if read_config == 2:
		discarded_read_dict = {
			'LFS end missing': set([]),
			'RFS start missing': set([]),
			'RFS end missing': set([]),
			'BBS start missing': set([])
			}

	no_marker_1_count = 0
	no_marker_2_count = 0
	no_marker_3_count = 0
	no_marker_4_count = 0
	invalid_sized_barcode_count = 0
	invalid_sized_coreseq_count = 0
	valid_read_count = 0
	discarded_read_count = 0
	total_read_count = 0
	barcode_sequence_list = []
	barcode_length_list = []
	first_read_list = []
	
	with open(step_02_output_filename, 'w') as filtered_fq_output, open(barcode_coreseq_list_filename, 'w') as barcode_coreseq_list_output:

		while 1:
			filtered_read = fastq_filtering_queue.get()

			if filtered_read == 'kill':
				
				with open(step_02_read_stats_filename, 'w') as raw_read_stats:

					raw_read_stats.write('no_lfs_end_count\t{}\n'.format(no_marker_1_count))
					raw_read_stats.write('no_rfs_start_count\t{}\n'.format(no_marker_2_count))
					
					if read_config == 1:
						raw_read_stats.write('no_bbs_end_count\t{}\n'.format(no_marker_3_count))
						raw_read_stats.write('no_lfs_start_count\t{}\n'.format(no_marker_4_count))
					
					if read_config == 2:
						raw_read_stats.write('no_rfs_end_count\t{}\n'.format(no_marker_3_count))
						raw_read_stats.write('no_bbs_start_count\t{}\n'.format(no_marker_4_count))

					raw_read_stats.write('invalid_sized_barcode_count\t{}\n'.format(invalid_sized_barcode_count))
					raw_read_stats.write('invalid_sized_coreseq_count\t{}\n'.format(invalid_sized_coreseq_count))
					raw_read_stats.write('valid_read_count\t{}\n'.format(valid_read_count))
					raw_read_stats.write('discarded_read_count\t{}\n'.format(discarded_read_count))
					raw_read_stats.write('total_read_count\t{}\n'.format(total_read_count))
					raw_read_stats.write('barcode_length_min\t{}\n'.format(np.min(barcode_length_list)))
					raw_read_stats.write('barcode_length_max\t{}\n'.format(np.max(barcode_length_list)))
					raw_read_stats.write('barcode_length_mean\t{}\n'.format(np.mean(barcode_length_list)))
					raw_read_stats.write('barcode_length_median\t{}\n'.format(np.median(barcode_length_list)))
					raw_read_stats.write('barcode_length_mode\t{}\n'.format(stats.mode(barcode_length_list)[0]))
				
					venn(discarded_read_dict, fmt = '{percentage:.1f}%', cmap = list('rgby'), fontsize = 10, legend_loc = 'upper left', figsize = (10,10))
					plt.savefig(step_02_discarded_reads_venn_filename, dpi = 600)
					plt.clf()

				break

			else:

				if filtered_read[3] == 1:
					no_marker_1_count += 1
					discarded_read_dict['LFS end missing'].add(filtered_read[0])
					discarded_read_dict['LFS end missing'].add(filtered_read[0])

				if filtered_read[4] == 1:
					no_marker_2_count += 1
					discarded_read_dict['RFS start missing'].add(filtered_read[0])
					discarded_read_dict['RFS start missing'].add(filtered_read[0])

				if filtered_read[5] == 1:
					no_marker_3_count += 1
					if read_config == 1:
						discarded_read_dict['BBS end missing'].add(filtered_read[0])
					if read_config == 2:
						discarded_read_dict['RFS end missing'].add(filtered_read[0])

				if filtered_read[6] == 1:
					no_marker_4_count += 1
					if read_config == 1:
						discarded_read_dict['LFS start missing'].add(filtered_read[0])
					if read_config == 2:
						discarded_read_dict['BBS start missing'].add(filtered_read[0])

				if filtered_read[7] != 'NA':
					barcode_sequence = filtered_read[7]
					barcode_length = len(filtered_read[7])
					barcode_length_list.append(barcode_length)

				if filtered_read[8] != 'NA':
					coreseq_sequence = filtered_read[8]
					coreseq_length = len(filtered_read[8])

				total_read_count += 1

				if all(error_state == 0 for error_state in filtered_read[3:7]) and barcode_length != 'NA':

					if len(barcode_length_list) == 1:
						print('Screening for the first 1000 reads with detectable barcodes to estimate a range of acceptable barcode sizes...')

					if len(barcode_length_list) <= 1000:
						first_read_list.append(filtered_read)

					if len(barcode_length_list) == 1000:
						mode_barcode_length = int(stats.mode(barcode_length_list)[0])
						barcode_length_floor = int(round(mode_barcode_length * 0.75))
						barcode_length_ceiling = int(round(mode_barcode_length * 1.5))
						
						print('Screening is finished. Mode barcode length is {}. Acceptable barcode sizes are between {} and {} base pairs.'.format(
							mode_barcode_length,
							barcode_length_floor, 
							barcode_length_ceiling))
							
						print('Reads with barcode size outside this range are listed below and subsequently discarded.')

						for first_read in first_read_list:
							invalid_size_error_state = 0

							if len(first_read[7]) < barcode_length_floor or len(first_read[7]) > barcode_length_ceiling:
								invalid_size_error_state = 1
								invalid_sized_barcode_count += 1

							if len(first_read[8]) < mode_barcode_length:
								invalid_size_error_state = 1
								invalid_sized_coreseq_count += 1

							if invalid_size_error_state == 0:
								barcode_sequence_list.append(first_read[7])
								if reference_exist == False:
									barcode_coreseq_list_output.write('{}\t{}\n'.format(first_read[7], first_read[8]))
								filtered_fq_output.write('{}\n{}\n+\n{}\n'.format(first_read[0], first_read[1], first_read[2]))
								valid_read_count += 1
							else:
								discarded_read_count += 1

					if len(barcode_length_list) > 1000:
						invalid_size_error_state = 0

						if barcode_length < barcode_length_floor or barcode_length > barcode_length_ceiling:
							invalid_size_error_state = 1
							invalid_sized_barcode_count += 1

						if coreseq_length < mode_barcode_length:
							invalid_size_error_state = 1
							invalid_sized_coreseq_count += 1

						if invalid_size_error_state == 0:
							barcode_sequence_list.append(barcode_sequence)
							if reference_exist == False:
								barcode_coreseq_list_output.write('{}\t{}\n'.format(barcode_sequence, coreseq_sequence))
							filtered_fq_output.write('{}\n{}\n+\n{}\n'.format(filtered_read[0], filtered_read[1], filtered_read[2]))
							valid_read_count += 1
						else:
							discarded_read_count += 1

				else:
					discarded_read_count += 1

			if reference_exist == False:
				barcode_coreseq_list_output.flush()
			filtered_fq_output.flush()

		return barcode_sequence_list, stats.mode(barcode_length_list)[0]


def reference_generator_function(collapsed_barcode): # The function to generate reference sequences for the automatically generated reference file
	coreseq_sequence_list = []
	barcode_instance_count = 0
	
	with open(barcode_coreseq_list_filename, 'r') as barcode_coreseq_list_input:
	
		while True:

			if barcode_instance_count < 100000: # Sample up to 100000 first valid reads with the exact same barcode, more than this is an unnecessary waste of memory
				
				barcode_coreseq_line = barcode_coreseq_list_input.readline()

				if not barcode_coreseq_line:
					break

				barcode_coreseq_combo = barcode_coreseq_line.strip('\n').split('\t')

				if barcode_coreseq_combo[0] == collapsed_barcode[0]:

					barcode_instance_count += 1 # Count up to 100000
					coreseq_sequence_list.append(barcode_coreseq_combo[1])

			else:
				break


	unique_coreseq_distance_average_list = []

	unique_coreseq_sequence_list, unique_coreseq_sequence_count_list = np.unique(coreseq_sequence_list, return_counts = True) # For reads with the same barcode, list all unique coreseq sequence and count

	for coreseq_sequence_a in unique_coreseq_sequence_list: # Among those with the right barcode length...
		unique_coreseq_distance_value_list = []

		for unique_coreseq_index, coreseq_sequence_b in enumerate(unique_coreseq_sequence_list): # ...iterate to check amongst themselves...
			unique_coreseq_distance = lev.distance(coreseq_sequence_a, coreseq_sequence_b) / np.mean([len(coreseq_sequence_a), len(coreseq_sequence_b)])

			for duplicate_counter in range(unique_coreseq_sequence_count_list[unique_coreseq_index]):
				unique_coreseq_distance_value_list.append(unique_coreseq_distance)

		unique_coreseq_distance_average_list.append(np.mean(unique_coreseq_distance_value_list))

	top_unique_coreseq_sequence_by_distance = unique_coreseq_sequence_list[np.argsort(unique_coreseq_distance_average_list)[0]]

	top_unique_coreseq_count_by_distance = unique_coreseq_sequence_count_list[unique_coreseq_sequence_list.tolist().index(top_unique_coreseq_sequence_by_distance)]

	top_unique_coreseq_sequence = top_unique_coreseq_sequence_by_distance

	top_unique_coreseq_count = top_unique_coreseq_count_by_distance
	
	if read_config == 1:
		reference_sequence = top_unique_coreseq_sequence + marker_1_seq + collapsed_barcode[0] + marker_2_seq # With the top coreseq sequence at hand, we can potentially create full reference sequence
		
	if read_config == 2:
		reference_sequence = marker_1_seq + collapsed_barcode[0] + marker_2_seq + top_unique_coreseq_sequence # With the top coreseq sequence at hand, we can potentially create full reference sequence
	
	# But! Given that one barcode exclusively marks one specific coreseq, we may be able to identify a variant barcode which is derived from the actual barcode but contains too much SNP to be correlatable by hamming distance 

	return [reference_sequence, collapsed_barcode[0], top_unique_coreseq_sequence, collapsed_barcode[1], top_unique_coreseq_count] # So return everything necessary for that purpose


def barcode_extracting_function(file_chunk, barcode_extracting_queue):
	for file_num in file_chunk:
		
		current_file_name = dir_files[file_num]

		# make a list of the extracted barcodes
		barcode_sequence_list = []

		# open the fastq file # first time it will genderate index which will take some time
		current_file = pyfastx.Fastq(current_file_name)
		
		# print(current_file_name, len(current_file))

		for read_num in range(len(current_file)):
			
			# pull out sequence
			current_read_seq = current_file[read_num].seq
			# figure out where the barcode is with respect to the flanking sequence
			# with pacbio can have indels and snps so need to be more lax to find barcode sequence
			# if the barcode lengths is way way too big (as seen it happen), then ignore it

			marker_1_index = current_read_seq.index(marker_1_seq)

			marker_2_index = current_read_seq.index(marker_2_seq)
			
			current_barcode = current_read_seq[(marker_1_index + search_length) : marker_2_index]
	
			barcode_sequence_list.append(current_barcode)

		# determine the unique barcodes and their counts
		unique_barcodes, counts = np.unique(barcode_sequence_list, return_counts=True)

		# make a df of the counts and report them to the reference barcode name
		extracted_data = {'Reference': (current_file_name.split('_')[-2] + '_' + current_file_name.split('_')[-1].split('.')[0]), 
						'Barcode': unique_barcodes, 
						'Barcode Count': counts}
		barcodes_dataframe = pd.DataFrame(data = extracted_data)
		
		# construct the filename for the output
		output_filename = '{}.barcode_counts.tsv'.format(current_file_name[:-5])

		# output the sorted barcodes out to file
		barcodes_dataframe.sort_values('Barcode Count', ascending = False).to_csv(output_filename, header = True, index = None, sep = '\t', mode = 'w')

		barcode_extracting_queue.put(barcodes_dataframe.values.tolist())

	return


def barcode_extracting_listener(barcode_extracting_queue):
	'''listens for messages on the q, writes to file. '''

	with open(step_05_global_barcode_list_filename, 'w') as global_barcode_output, open(step_05_global_collapsed_barcode_list_filename, 'w') as global_collapsed_barcode_output:
		global_writer = csv.writer(global_barcode_output, delimiter = '\t')
		global_collapsed_writer = csv.writer(global_collapsed_barcode_output, delimiter = '\t')

		while 1:
			barcodes_list = barcode_extracting_queue.get()

			if barcodes_list == 'kill':
				break

			barcodes_list_transposed = list(map(list, zip(*barcodes_list)))

			global_writer.writerows(barcodes_list)
			global_collapsed_writer.writerow([barcodes_list_transposed[0][np.argmax(barcodes_list_transposed[2])], sum(barcodes_list_transposed[2])])

			global_barcode_output.flush()
			global_collapsed_barcode_output.flush()

###############################################################################

parser = argparse.ArgumentParser(description = 'Filter, Trim, and Split PacBio uBAM or FASTQ.GZ based on Barcode and Identifier Sequences.')

parser.add_argument('--arglist', required = False, 
					help = 'Path to argument list for the command line. The each argument in this list will be overridden if the flag argument is provided in the command line.', action = 'store')
parser.add_argument('--input', required = False, 
					help = 'Path to input uBAM or fq.gz file.', action = 'store')
parser.add_argument('--outdir', required = False, 
					help = 'Path to output directory.', action = 'store')
parser.add_argument('--setname', required = False, 
					help = 'The name to be used for the processed sample.', action = 'store')
parser.add_argument('--mapref', required = False,
					help = 'Path to the reference file used by minimap2 for reads mapping. Mapping reference will be auto-generated if this file is not provided.', action = 'store')
parser.add_argument('--autoref', required = False, 
					help = 'The method that will be used for auto referencing when reference file (mapref) is not provided,', action = 'store', choices = ['slope','distance'])
parser.add_argument('--thread', required = False, type = int,
					help = 'Number of cores to use for multi-threaded applications. Will be set to half the maximum number of cores available if not specified.', action = 'store',
					choices = range(1, (multiprocess.cpu_count() + 1), 1),
                    metavar = '[1-{}]'.format(multiprocess.cpu_count()))
parser.add_argument('--config', required = False, 
					help = 'Your long reads construct configuration. Choose 1 if BBS-CoreSeq-LFS-Barcode-RFS. Choose 2 if it is LFS-Barcode-RFS-CoreSeq-BBS.', action = 'store', choices = [1,2])
parser.add_argument('--bbs', required = False, type = str,
					help = 'Left (if config 1) or right (if config 2) BackBone Sequence (BBS). Any letter case is accepted.', action = 'store')
parser.add_argument('--lfs', required = False, type = str,
					help = 'Left barcode Flanking Sequence (LFS). Any letter case is accepted.', action = 'store')
parser.add_argument('--rfs', required = False, type = str,
					help = 'Right barcode Flanking Sequence (RFS). Any letter case is accepted.', action = 'store')
parser.add_argument('--slen', required = False, type = int,
					help = 'Search length for BBS, LFS, and RFS. Recommended value is 20 or above. Cannot be longer than BBS, LFS, or RFS.', action = 'store')
parser.add_argument('--mapq', required = False, type = int,
					help = 'The minimum MAPQ score of reads to keep after minimap2 alignment. Will be set to default value of 5 if not specified.', action = 'store',
					choices = range(0, 40, 1),
                    metavar = '[0-40]')

args=parser.parse_args()

###############################################################################

if args.arglist:
	arglist = os.path.abspath(args.arglist)

	argument_df = pd.read_csv(arglist, delimiter = '\t', header =  None)

	argument_list = ['input',
					'outdir',
					'setname',
					'mapref',
					'autoref',
					'thread',
					'config',
					'bbs',
					'lfs',
					'rfs',
					'slen',
					'mapq']

	argument_array = argument_df.values.tolist()

	argument_dict = {}

	for argument in argument_list:
		argument_dict[argument] = []

	for argument_array_row in argument_array:
		if argument_array_row[0] in argument_list:
			current_flag = argument_array_row[0]
			current_argument = argument_array_row[1]
			argument_dict[current_flag].append(current_argument)

	try:
		input             		= str(argument_dict['input'][0])
	except:
		input					= 'nan'
	try:
		outdir        			= str(argument_dict['outdir'][0])
	except:
		outdir					= 'nan'
	try:
		setname          		= str(argument_dict['setname'][0])
	except:
		setname					= 'nan'
	try:
		minimap_reference_path	= str(argument_dict['mapref'][0])
	except:
		minimap_reference_path 	= 'nan'
	try:
		auto_reference_method 	= str(argument_dict['autoref'][0])
	except:
		auto_reference_method	= 'nan'
	try:
		cpu_count    			= int(argument_dict['thread'][0])
	except:
		cpu_count				= 'nan'
	try:
		read_config				= int(argument_dict['config'][0])
	except:
		read_config				= 'nan'
	try:
		bbs_sequence     		= str(argument_dict['bbs'][0])
	except:
		bbs_sequence			= 'nan'
	try:
		lfs_sequence        	= str(argument_dict['lfs'][0])
	except:
		lfs_sequence			= 'nan'
	try:
		rfs_sequence     		= str(argument_dict['rfs'][0])
	except:
		rfs_sequence			= 'nan'	
	try:
		search_length    		= int(argument_dict['slen'][0])
	except:
		search_length    		= 'nan'
	try:
		mapq_score    			= int(argument_dict['mapq'][0])
	except:
		mapq_score				= 'nan'

###############################################################################

missing_input_status = 0

print('\n\n||| {} - SLICER - Program Started |||\n\n'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())))

if args.thread:
	cpu_count = args.thread
	print('User-determined number of processes: {}.'.format(cpu_count))
else:
	if cpu_count != 'nan':
		print('User-determined number of processes: {}.'.format(cpu_count))
	else:
		print(Fore.GREEN + 'The value for variable "thread" is not provided by user. Setting it to {} cores.'.format(math.ceil(multiprocess.cpu_count() / 2)) + Style.RESET_ALL)
		cpu_count = math.ceil(multiprocess.cpu_count() / 2)

if cpu_count > multiprocess.cpu_count():
	print(Fore.GREEN + 'The value for variable "thread" exceeds maximum available cores. Setting it to {} cores.'.format(multiprocess.cpu_count()) + Style.RESET_ALL)
	cpu_count = math.ceil(multiprocess.cpu_count())

print('')

if args.input:
	input = os.path.abspath(args.input)
	if os.path.isfile(input):
		print('Path to the input uBAM file: {}.'.format(input))
	else:
		print(Fore.RED + 'Unable to find {}. Please check if the path is correct!'.format(input) + Style.RESET_ALL)
		missing_input_status = 1
else:
	if input != 'nan':
		if os.path.isfile(input):
			print('Path to the input uBAM file: {}.'.format(input))
		else:
			print(Fore.RED + 'Unable to find {}. Please check if the path is correct!'.format(input) + Style.RESET_ALL)
			missing_input_status = 1
	else:
		print(Fore.RED + 'Please provide the value for variable "input" either through the argument list or command line flag argument!' + Style.RESET_ALL)
		missing_input_status = 1

print('')

if args.outdir:
	outdir = os.path.abspath(args.outdir)
	outdir = outdir.rstrip('/')
	if os.path.isdir(outdir):
		print('Path to the output directory: {}.'.format(outdir))
	else:
		print(Fore.RED + 'Unable to find {}. Please check if the path is correct!'.format(outdir) + Style.RESET_ALL)
		missing_input_status = 1
else:
	if outdir != 'nan':
		outdir = outdir.rstrip('/')
		if os.path.isdir(outdir):
			print('Path to the output directory: {}.'.format(outdir))
		else:
			print(Fore.RED + 'Unable to find {}. Please check if the path is correct!'.format(outdir) + Style.RESET_ALL)
			missing_input_status = 1
	else:
		print(Fore.RED + 'Please provide the value for variable "outdir" either through the argument list or command line flag argument!' + Style.RESET_ALL)
		missing_input_status = 1

print('')

if args.setname:
	setname = args.setname
	print('Dataset name: {}.'.format(setname))
else:
	if setname != 'nan':
		print('Dataset name: {}.'.format(setname))
	else:
		print(Fore.RED + 'Please provide the value for variable "setname" either through the argument list or command line flag argument!' + Style.RESET_ALL)
		missing_input_status = 1

print('')

if args.mapref:
	minimap_reference_path = os.path.abspath(args.mapref)
	if os.path.isfile(minimap_reference_path):
		reference_exist = True
		print('Path to the mapping reference FASTA file: {}.'.format(minimap_reference_path))
	else:
		print(Fore.RED + 'Unable to find {}. Please check if the path is correct!'.format(minimap_reference_path) + Style.RESET_ALL)
		missing_input_status = 1
else:
	if minimap_reference_path != 'nan':
		if os.path.isfile(minimap_reference_path):
			reference_exist = True
			print('Path to the mapping reference FASTA file: {}.'.format(minimap_reference_path))
		else:
			print(Fore.RED + 'Unable to find {}. Please check if the path is correct!'.format(minimap_reference_path) + Style.RESET_ALL)
			missing_input_status = 1
	else:
		if outdir != 'nan':
			minimap_reference_path = '{}/02_automatic_reference_{}.fa'.format(outdir, setname)
			reference_exist = False
			print(Fore.GREEN + 'Reads alignment reference path not assigned, the program will automatically generate a reference file.' + Style.RESET_ALL)
			print(Fore.GREEN + 'You can later check the automatically generated reference here: "{}".'.format(minimap_reference_path) + Style.RESET_ALL)
			if args.autoref:
				auto_reference_method = args.autoref
				print('Selected reference prediction method: "{}"'.format(auto_reference_method))
			else:
				if auto_reference_method != 'nan':
					if auto_reference_method == 'slope' or auto_reference_method == 'distance':
						print('Selected reference prediction method: "{}"'.format(auto_reference_method))
					else:
						print(Fore.RED + 'Invalid reference prediction method: "{}". Please only choose between "slope" and "distance".'.format(auto_reference_method) + Style.RESET_ALL)
						missing_input_status = 1
				else:
					auto_reference_method = 'slope'
					print(Fore.GREEN + 'Reference prediction method is not set. "Slope" method will be used as default.' + Style.RESET_ALL)		
		else:
			print(Fore.RED + 'Unable to create automatically generated reference file because "outdir" has not been assigned yet.' + Style.RESET_ALL)

print('')

if args.config:
	read_config = args.config
	if read_config == 1:
		component_order = 'BBS-CoreSeq-LFS-Barcode-RFS'
	if read_config == 2:
		component_order = 'LFS-Barcode-RFS-CoreSeq-BBS'
	print('Reads configuration: Type {} ({}).'.format(read_config, component_order))
else:
	if read_config != 'nan':
		if read_config == 1:
			component_order = 'BBS-CoreSeq-LFS-Barcode-RFS'
		if read_config == 2:
			component_order = 'LFS-Barcode-RFS-CoreSeq-BBS'
		print('Reads configuration: Type {} ({}).'.format(read_config, component_order))
	else:
		print(Fore.RED + 'Please provide the value for variable "config" either through the argument list or command line flag argument!' + Style.RESET_ALL)
		missing_input_status = 1

print('')

if args.bbs:
	bbs_sequence = args.bbs
	if read_config == 1:
		print('Left BackBone Sequence (BBS) is: {}.'.format(bbs_sequence))
	elif read_config == 2:
		print('Right BackBone Sequence (BBS) is: {}.'.format(bbs_sequence))
	else:
		print('BackBone Sequence (BBS) is: {}.'.format(bbs_sequence))
else:
	if bbs_sequence != 'nan':
		if read_config == 1:
			print('Left BackBone Sequence (BBS) is: {}.'.format(bbs_sequence))
		elif read_config == 2:
			print('Right BackBone Sequence (BBS) is: {}.'.format(bbs_sequence))
		else:
			print('BackBone Sequence (BBS) is: {}.'.format(bbs_sequence))
	else:
		print(Fore.RED + 'Please provide the value for variable "bbs" either through the argument list or command line flag argument!' + Style.RESET_ALL)
		missing_input_status = 1

print('')

if args.lfs:
	lfs_sequence = args.lfs
	print('Left Flanking Sequence (LFS) is: {}.'.format(lfs_sequence))
else:
	if lfs_sequence != 'nan':
		print('Left Flanking Sequence (LFS) is: {}.'.format(lfs_sequence))
	else:
		print(Fore.RED + 'Please provide the value for variable "lfs" either through the argument list or command line flag argument!' + Style.RESET_ALL)
		missing_input_status = 1

print('')

if args.rfs:
	rfs_sequence = args.rfs
	print('Right Flanking Sequence (RFS) is: {}.'.format(rfs_sequence))
else:
	if rfs_sequence != 'nan':
		print('Right Flanking Sequence (RFS) is: {}.'.format(rfs_sequence))
	else:
		print(Fore.RED + 'Please provide the value for variable "rfs" either through the argument list or command line flag argument!' + Style.RESET_ALL)
		missing_input_status = 1

print('')

if args.slen:
	search_length = args.slen
	print('SLICER will search all reads for {}-base matches toward identifier sequences (LFS, RFS, BackBone).'.format(search_length))
else:
	if search_length != 'nan':
		print('SLICER will search all reads for {}-base matches toward identifier sequences (LFS, RFS, BackBone).'.format(search_length))
	else:
		print(Fore.RED + 'Please provide the value for variable "slen" either through the argument list or command line flag argument!' + Style.RESET_ALL)
		missing_input_status = 1

print('')

if args.mapq:
	mapq_score = args.mapq
	print('User-determined minimum MAPQ score: {}.'.format(mapq_score))
else:
	if mapq_score != 'nan':
		print('User-determined minimum MAPQ score: {}.'.format(mapq_score))
	else:
		print(Fore.GREEN + 'The value for variable "mapq" is not provided by user. Setting it to 5.' + Style.RESET_ALL)
		mapq_score = 5

if mapq_score > 40:
	print(Fore.GREEN + 'The value for variable "mapq" exceeds maximum possible MAPQ score. Setting it to 40.' + Style.RESET_ALL)
	mapq_score = 40

print('')

if missing_input_status == 0:
	pass
if missing_input_status == 1:    
	print('\nPlease refer to the manual below and provide missing variable(s) or incorrect path(s):\n')
	parser.print_help()
	print('\n')
	exit()

###############################################################################

if not os.path.isdir(outdir):
	os.mkdir(outdir)

if len(bbs_sequence) < search_length:
	print(Fore.RED + 'Search length (slen) value exceeds the size of BBS (bbs) sequence. Please provide smaller slen value or longer bbs sequence.' + Style.RESET_ALL)
	exit()

if len(lfs_sequence) < search_length:
	print(Fore.RED + 'Search length (slen) value exceeds the size of LFS (lfs) sequence. Please provide smaller slen value or longer lfs sequence.' + Style.RESET_ALL)
	exit()

if len(rfs_sequence) < search_length:
	print(Fore.RED + 'Search length (slen) value exceeds the size of RFS (rfs) sequence. Please provide smaller slen value or longer rfs sequence.' + Style.RESET_ALL)
	exit()
	
# Marker 1 and 2 are used for parsing barcode position in both construct configurations, which are always between LFS end and RFS start
marker_1_seq = lfs_sequence.upper()[-search_length:]
marker_2_seq = rfs_sequence.upper()[:search_length]
print('{} ({}-bases end of the Left Flanking Sequence) will be used to parse for the barcode starting position.'.format(marker_1_seq, search_length))
print('{} ({}-bases start of the Right Flanking Sequence) will be used to parse for the barcode ending position.'.format(marker_2_seq, search_length))

# Marker 3 and 4 are used for parsing core sequence position in both construct configurations, but the flanking elements are different between these configurations
if read_config == 1:
	# In configuration 1, core sequence is located between the left BBS end and the LFS start
	marker_3_seq = bbs_sequence.upper()[-search_length:]
	marker_4_seq = lfs_sequence.upper()[:search_length]
	print('{} ({}-bases end of the Left BackBone Sequence) will be used to parse for the core sequence starting position.'.format(marker_3_seq, search_length))
	print('{} ({}-bases start of the Left Flanking Sequence) will be used to parse for the core sequence ending position.'.format(marker_4_seq, search_length))

if read_config == 2:
	# In configuration 2, core sequence is located between the RFS end and the right BBS start
	marker_3_seq = rfs_sequence.upper()[-search_length:]
	marker_4_seq = bbs_sequence.upper()[:search_length]
	print('{} ({}-bases end of the Right Flanking Sequence) will be used to parse for the core sequence starting position.'.format(marker_3_seq, search_length))
	print('{} ({}-bases start of the Right BackBone Sequence) will be used to parse for the core sequence ending position.'.format(marker_4_seq, search_length))

print('\n')

# Change default font size for plotted graphs
plt.rcParams.update({'font.size': 8})

# Move to the output directory
os.chdir(outdir)

###############################################################################

print('||| Step 1 - Convert uBAM to FASTQ for Manipulation |||\n')
start = time.time()

for file_name in listdir():
	if file_name.endswith('.fxi'):
		os.remove(file_name)

# Filename of the output of step 1
step_01_output_filename = '01_{}.fq.gz'.format(setname)

if input.split('.')[-1].lower() == 'bam' or input.split('.')[-1].lower() == 'ubam':
	if not os.path.isfile(step_01_output_filename):

		# bamToFastq terminal command
		print('Converting input .bam or .ubam to .fq.gz format...')
		bamtofastq_command = 'bamToFastq -i {} -fq /dev/stdout | gzip > {}'.format(input, step_01_output_filename)
		print('\nExecuting: {}...\n'.format(bamtofastq_command))
		os.system(bamtofastq_command)

elif input.split('.')[-2].lower() == 'fq' and input.split('.')[-1].lower() == 'gz':
	step_01_output_filename = input.split('/')[-1]
	print('File already in .fq.gz format, skipping conversion step.')

elif input.split('.')[-2].lower() == 'fastq' and input.split('.')[-1].lower() == 'gz':
	step_01_output_filename = input.split('/')[-1]
	print('File already in .fq.gz format, skipping conversion step.')

else:
	print('Can only process input under extension .bam, .ubam, .fq.gz, or .fastq.gz (case insensitive)')
	exit()

end = time.time()
print('--- Elapsed time: {} seconds. ---\n\n'.format(end - start))

###############################################################################

print('||| Step 2 - Filter and Trim FASTQ File to Appropriate Sequences |||\n')
start = time.time()

# filename of the output of step 2
step_02_read_stats_filename 			= '02_raw_read_stats_{}.tsv'.format(setname)
step_02_discarded_reads_venn_filename 	= '02_discarded_reads_venn_{}.png'.format(setname)
step_02_barcode_stats_filename 			= '02_barcode_stats_{}.tsv'.format(setname)
step_02_auto_reference_filename 		= '02_automatic_reference_{}.fa'.format(setname)
step_02_output_filename 				= '02_filtered_{}.fq'.format(setname)
barcode_coreseq_list_filename 			= 'tmp_barcode_coreseq_list_{}.txt'.format(setname)

os.system('rm -r -f {}.fxi'.format(step_01_output_filename))

# Create the main index for the FASTQ file
start_index = time.time()
raw_fq = pyfastx.Fastq(step_01_output_filename)
end_index = time.time()
# print('Main process FASTQ object successfully generated in {} seconds'.format(end_index - start_index))

raw_fq_len = len(raw_fq)
raw_fq_chunks = np.array_split(range(raw_fq_len), math.ceil(cpu_count / 2))

# open the output fq file for saving
raw_read_stats = open(step_02_read_stats_filename, 'w')
filtered_fq_output = open(step_02_output_filename, 'w')
if reference_exist == False:
	barcode_coreseq_list_output = open(barcode_coreseq_list_filename, 'w')


fastq_filtering_manager = multiprocess.Manager()
fastq_filtering_queue = fastq_filtering_manager.Queue()
fastq_filtering_pool = multiprocess.Pool(processes = math.ceil(cpu_count / 2) + 1)

#put listener to work first
fastq_filtering_watcher = fastq_filtering_pool.apply_async(fastq_filtering_listener, (fastq_filtering_queue,))

#fire off workers
fastq_filtering_job_list = []
for core_counter in range(math.ceil(cpu_count / 2)):
	fastq_filtering_job = fastq_filtering_pool.apply_async(fastq_filtering_function, ([list(raw_fq_chunks[core_counter]), core_counter], fastq_filtering_queue))
	fastq_filtering_job_list.append(fastq_filtering_job)

# collect results from the workers through the pool result queue
for fastq_filtering_job in fastq_filtering_job_list: 
	fastq_filtering_job.get()

#now we are done, kill the listener
fastq_filtering_queue.put('kill')
fastq_filtering_pool.close()
fastq_filtering_pool.join()

raw_read_stats.close()
filtered_fq_output.close()
if reference_exist == False:
	barcode_coreseq_list_output.close()

barcode_sequence_result, barcode_length_mode = fastq_filtering_watcher.get()


# gz the fq file for storage
fastq_zip_command = 'pigz --best -k -f -p {} {}'.format(cpu_count, step_02_output_filename)
print('\nExecuting: {}...\n'.format(fastq_zip_command))
os.system(fastq_zip_command)

unique_barcode_sequence_list, unique_barcode_count_list = np.unique(barcode_sequence_result, return_counts=True)

barcode_df = pd.DataFrame(list(zip(unique_barcode_sequence_list, unique_barcode_count_list))) # All unprocessed barcodes, parsed from all valid reads based on the location of LFS end and RFS start; and the number of instances found
barcode_df.columns = ['Barcode', 'Barcode Count']
barcode_df.sort_values(by = ['Barcode Count'], ascending=False, inplace=True) # Sort barcodes based on the number of instances found. Largest on top.
barcode_df['Barcode Length'] = barcode_df.apply(lambda x: len(x.Barcode), axis=1) # Append the length (number of bases) to each barcode on the list
barcode_df.to_csv(step_02_barcode_stats_filename, header=True, index=None, sep='\t', mode='w') # Save the data as a table

end = time.time()
print('--- Elapsed time: {} seconds. ---\n'.format(end - start))


if reference_exist == False:

	print('Mapping reference file was not provided. A reference file will be predicted by parsing the filtered reads.')
	print('Parsing barcodes in sequenced reads...')

	start = time.time()

	true_barcode_list = []

	selected_barcode_list = []
	collapsed_barcode_list = []
	collapsed_barcode_sequence_list = []

	selected_reference_list = []
	collapsed_reference_list = []
	collapsed_reference_sequence_list = []


	if auto_reference_method == 'slope':

		print('Slope method is going to be used for automatic reference generation.')
		print('Trying to automatically identify barcode sequences...')
		print('Assumption: barcode actual length = parsed barcodes mode length = {}.'.format(barcode_length_mode))

		filtered_barcode_df = barcode_df.loc[barcode_df['Barcode Length'] == barcode_length_mode]
		filtered_barcode_list = filtered_barcode_df[['Barcode', 'Barcode Count']].values.tolist()
		filtered_barcode_count_list = filtered_barcode_df['Barcode Count'].values.tolist()
		filtered_barcode_sequence_list = filtered_barcode_df['Barcode'].values.tolist()

		nonmode_barcode_count_list = barcode_df.loc[barcode_df['Barcode Length'] != barcode_length_mode]['Barcode Count'].values.tolist()

		barcode_count_slope_list = []

		for sorted_barcode_count_counter in range(len(filtered_barcode_count_list) - 1):
			barcode_count_slope = filtered_barcode_count_list[sorted_barcode_count_counter] / filtered_barcode_count_list[sorted_barcode_count_counter + 1]
			barcode_count_slope_list.append(barcode_count_slope)

		print('Maximum slope value: {}.'.format(round(max(barcode_count_slope_list), 3)))
		print('Maximum slope at index number: {}.'.format(np.argmax(barcode_count_slope_list)))
		print('Slope from the minimum mode-lengthed barcode count to the maximum non-mode lengthed barcode count: {}.'.format(
			round((min(filtered_barcode_count_list) / max(nonmode_barcode_count_list)), 3)))

		if (min(filtered_barcode_count_list) / max(nonmode_barcode_count_list)) > max(barcode_count_slope_list):
			number_of_reference = len(barcode_count_slope_list) + 1
		else:
			number_of_reference = np.argmax(barcode_count_slope_list) + 1

		print('Predicted number of barcodes based on slope method: {}.'.format(number_of_reference))
		
		collapsed_barcode_sequence_list = filtered_barcode_sequence_list[:number_of_reference]
		collapsed_barcode_list = filtered_barcode_list[:number_of_reference]


	if auto_reference_method == 'distance':

		print('Distance method is going to be used for automatic reference generation.')
		print('Trying to automatically identify barcode sequences...')
		print('Assumption: barcode actual length = parsed barcodes mode length = {}.'.format(barcode_length_mode))

		for barcode_data_row in list(zip(unique_barcode_sequence_list, unique_barcode_count_list)):
			if len(barcode_data_row[0]) == barcode_length_mode: # Parsed barcodes with length equal to mode...
				true_barcode_list.append(barcode_data_row) # ...are considered to have the right barcode length, and (potentially) have the actual barcode sequences

		print('Number of potential barcodes (with mode length): {}.'.format(len(true_barcode_list)))


		print('Calculating pairwise hamming distance between barcodes and grouping together barcodes with potential SNP variations...')

		for true_barcode_a in true_barcode_list: # Among those with the right barcode length...
			similar_barcode_list = []
			
			for true_barcode_b in true_barcode_list: # ...iterate to check amongst themselves...
				hamming_distance = hamming(list(true_barcode_a[0]), list(true_barcode_b[0]))
				if hamming_distance <= 0.25: # ...for sequence similarity by means of pairwise hamming distance
					similar_barcode_list.append(true_barcode_b) # If hamming distance is small enough between two potential barcode sequences, one of the pair is assumed to be SNP variants of the other 

			unzipped_similar_barcode_list = list(map(list, zip(*similar_barcode_list))) # Unzip to separate barcode sequence and count lists
			maxcount_similar_barcode_index = np.argmax(unzipped_similar_barcode_list[1]) # Find the index of maximum value in the count list
			selected_barcode = [unzipped_similar_barcode_list[0][maxcount_similar_barcode_index], sum(unzipped_similar_barcode_list[1])] # Group together all barcode counts under the majority sequence
			selected_barcode_list.append(selected_barcode) # For each barcode pair sharing similarity, only the majority sequence are selected and all minority sequences will be disregarded


		print('Grouping together and removing redundant entries with identical barcode...')

		for selected_barcode_a in selected_barcode_list: # Among those selected majority sequences...
			redundant_barcode_count_list = []
			for selected_barcode_b in selected_barcode_list: # ...iterate to check amongst themselves...
				if selected_barcode_a[0] == selected_barcode_b[0]: # ...for redundant entries caused by reciprocal search result (from e.g. A finds B similar and B finds A similar)
					redundant_barcode_count_list.append(selected_barcode_b[1]) # For entries with identical barcode sequence, record every total count value
			
			# Since redundant entries should have the same majority sequence, only one will be needed and the rest can be discarded
			# However, redundant entries may have different total count (e.g. when A finds B and C similar, but both B and C only find A similar and not each other --> total count of A will be the highest)
			
			redundant_barcode_max = max(redundant_barcode_count_list) # Thus, take the highest total count, since that would cover the barcode counts from all possible matches
			collapsed_barcode = [selected_barcode_a[0], redundant_barcode_max] # Therefore, right here, one element contains one barcode sequence + the highest total count

			# At this point, list_b no longer have redundant entries
			# However, there are still redundant entries from list_a, and this needs to be filtered out

			if collapsed_barcode[0] not in collapsed_barcode_sequence_list: # Therefore, only add the barcode to the list if it is not already there from any previous entry in list_a
				collapsed_barcode_sequence_list.append(collapsed_barcode[0]) # Checklist, just to prevent duplicate entries of the same barcode
				collapsed_barcode_list.append(collapsed_barcode) # Result: List of barcodes grouped and collapsed based on barcode sequence similarities.

		print('Number of potential barcodes (after merging of similar barcodes based on Hamming distance): {}.'.format(len(collapsed_barcode_list)))

	end = time.time()
	print('--- Elapsed time: {} seconds. ---\n'.format(end - start))


	start = time.time()
	print('Generating references...')

	# Execute the function: generate reference sequences for the automatically generated reference file
	pool_reference_generator_function = multiprocess.Pool(processes = cpu_count)
	reference_data_array = pool_reference_generator_function.map(reference_generator_function, collapsed_barcode_list)
	pool_reference_generator_function.close()
	pool_reference_generator_function.join()

	print('Number of potential references (with unique barcode + core sequence combination): {}.'.format(len(reference_data_array)))


	if auto_reference_method == 'slope':

		collapsed_reference_list = reference_data_array


	if auto_reference_method == 'distance':

		print('Calculating pairwise levenshtein distance between core sequences and calling for potential indel or SNP variations...')

		# Before trying to identify variant barcodes based on their coreseq sequence, we should first group together and collapse indel or SNP close variants of the coreseqs we have
		for reference_data_array_row_a in reference_data_array: # Among those reference sequences with consensus coreseq sequences...
			similar_coreseq_list = []
			for reference_data_array_row_b in reference_data_array: # ...iterate to check amongst themselves...
				levenshtein_distance = lev.distance(reference_data_array_row_a[2], reference_data_array_row_b[2]) / np.mean([len(reference_data_array_row_a[2]), len(reference_data_array_row_b[2])])
				if levenshtein_distance <= 0.012: # ...for sequence similarity by means of pairwise levenshtein distance
					similar_coreseq_list.append(reference_data_array_row_b) # If levenshtein distance is small enough between two potential coreseq sequences, one of the pair is assumed to be indel or SNP variants of the other


			unzipped_similar_coreseq_list = list(map(list, zip(*similar_coreseq_list))) # Unzip to separate reference data and the coreseq count lists
			maxcount_similar_coreseq_index = np.argmax(unzipped_similar_coreseq_list[4]) # Find the index of maximum value in the coreseq count list
			selected_reference = [unzipped_similar_coreseq_list[0][maxcount_similar_coreseq_index],
								unzipped_similar_coreseq_list[1][maxcount_similar_coreseq_index],
								unzipped_similar_coreseq_list[2][maxcount_similar_coreseq_index], 
								sum(unzipped_similar_coreseq_list[3]), # Group together all barcode counts under the majority sequence
								sum(unzipped_similar_coreseq_list[4])] # Group together all coreseq counts under the majority sequence
			selected_reference_list.append(selected_reference) # For each reference pair sharing similar coreseq sequence, only the majority sequence are selected and all minority sequences will be disregarded


		print('Grouping together and removing redundant references with identical core sequences...')

		for selected_reference_a in selected_reference_list: # Among those selected majority sequences...
			redundant_barcode_count_list = []
			redundant_coreseq_count_list = []
			for selected_reference_b in selected_reference_list: # ...iterate to check amongst themselves...
				if selected_reference_a[2] == selected_reference_b[2]: # ...for redundant entries caused by reciprocal search result (from e.g. A finds B similar and B finds A similar)
					redundant_barcode_count_list.append(selected_reference_b[3]) # For entries with identical reference sequence, record every total count value
					redundant_coreseq_count_list.append(selected_reference_b[4]) # For entries with identical reference sequence, record every total count value
			
			# Since redundant entries should have the same majority set of sequences, only one will be needed and the rest can be discarded
			# However, redundant entries may have different total count (e.g. when A finds B and C similar, but both B and C only find A similar and not each other --> total count of A will be the highest)
			
			redundant_barcode_max = max(redundant_barcode_count_list) # Thus, take the highest total count, since that would cover the reference counts from all possible matches
			redundant_coreseq_max = max(redundant_coreseq_count_list) # Thus, take the highest total count, since that would cover the reference counts from all possible matches
			collapsed_reference = [selected_reference_a[0], selected_reference_a[1], selected_reference_a[2], redundant_barcode_max, redundant_coreseq_max] # Therefore, right here, one element contains one set of reference sequences + the highest total count

			# At this point, list_b no longer have redundant entries
			# However, there are still redundant entries from list_a, and this needs to be filtered out

			if collapsed_reference[2] not in collapsed_reference_sequence_list: # Therefore, only add the reference to the list if the coreseq sequence is not already there from any previous entry in list_a
				collapsed_reference_sequence_list.append(collapsed_reference[2]) # Checklist, just to prevent duplicate entries of the same reference
				collapsed_reference_list.append(collapsed_reference) # Result: List of references grouped and collapsed based on coreseq sequence similarities.

		print('Number of potential barcodes (after merging of references with similar core sequences based on Levenshtein distance): {}.'.format(len(collapsed_reference_list)))
		print('Predicted number of barcodes based on distance method: {}.'.format(len(collapsed_reference_list)))


	collapsed_reference_list_transposed = list(map(list, zip(*collapsed_reference_list)))
	
	collapsed_reference_dict = {'Full Sequence': collapsed_reference_list_transposed[0], 
								'Barcode Sequence': collapsed_reference_list_transposed[1],
								'Core Sequence': collapsed_reference_list_transposed[2], 
								'Barcode Count': collapsed_reference_list_transposed[3],
								'Coreseq Count': collapsed_reference_list_transposed[4]}

	collapsed_reference_df = pd.DataFrame(data = collapsed_reference_dict)

	collapsed_reference_df.sort_values(['Barcode Count', 'Barcode Sequence'], ascending = [False, True], inplace = True)
	collapsed_reference_list_updated = collapsed_reference_df.values.tolist()

	with open(step_02_auto_reference_filename, 'w') as auto_reference_output:
		for reference_number, final_reference_reference_data in enumerate(collapsed_reference_list_updated):
			final_reference_full_sequence = final_reference_reference_data[0]
			final_reference_barcode_sequence = final_reference_reference_data[1]
			final_reference_coreseq_sequence = final_reference_reference_data[2]
			final_reference_barcode_count = final_reference_reference_data[3]
			final_reference_coreseq_count = final_reference_reference_data[4]
			# print(final_reference_barcode_sequence, final_reference_barcode_count, final_reference_coreseq_count)

			auto_reference_output.write('>Ref{}_{}\n'.format(reference_number + 1, final_reference_barcode_sequence))
			auto_reference_output.write('{}\n'.format(final_reference_full_sequence))

	end = time.time()
	print('--- Elapsed time: {} seconds. ---\n\n'.format(end - start))

os.remove(barcode_coreseq_list_filename)

###############################################################################

print('||| Step 3 - Align, Sort, Index and MAPQ5 Filter Aligned Reads |||\n')
start = time.time()

# filename of the output of step 2
step_03_output_filename = '03_aligned_filtered_{}.sam'.format(setname)

# map to reference using minimap2
minimap_alignment_command = 'minimap2 -ax splice:hq -t {} -2 --secondary=no -L {} {} -o {}_inclSuppl.sam'.format(
	cpu_count, minimap_reference_path, step_02_output_filename, step_03_output_filename[:-4])
print('Executing: {}...\n'.format(minimap_alignment_command))
os.system(minimap_alignment_command)

# filter out supplemental/chimeric reads from alignments
samtools_filter_command = 'samtools view -@ {} -F 2048 -bo {} {}_inclSuppl.sam'.format(cpu_count, step_03_output_filename, step_03_output_filename[:-4])
print('\nExecuting: {}...\n'.format(samtools_filter_command))
os.system(samtools_filter_command)

# convert, sort and index bam
samtools_view_command = 'samtools view -Sb -@ {} {}_inclSuppl.sam > {}_inclSuppl.bam'.format(cpu_count, step_03_output_filename[:-4], step_03_output_filename[:-4])
print('\nExecuting: {}...\n'.format(samtools_view_command))
os.system(samtools_view_command)

samtools_sort_command = 'samtools sort -@ {} {} -o {}.sorted.bam'.format(cpu_count, step_03_output_filename, step_03_output_filename[:-4])
print('\nExecuting: {}...\n'.format(samtools_sort_command))
os.system(samtools_sort_command)

samtools_index_command = 'samtools index -@ {} {}.sorted.bam'.format(cpu_count, step_03_output_filename[:-4])
print('\nExecuting: {}\n...'.format(samtools_index_command))
os.system(samtools_index_command)

# clean up unnecessary sam files
os.system('rm {}'.format(step_03_output_filename))
os.system('rm {}_inclSuppl.sam'.format(step_03_output_filename[:-4]))

#########################################

# generate quality plots before mapq filtering
bam = pysam.AlignmentFile('{}.sorted.bam'.format(step_03_output_filename[:-4]), 'rb')

quals = [read.mapping_quality for read in bam.fetch()]

sns.set_style('dark')
sns.set_context('paper')

plt.figure(figsize=(10,6))

ax = sns.kdeplot(quals, fill = True)
ax.set_ylabel('Density')
ax.set_xlabel('Mapping Quality Score')

plt.savefig('{}.sorted.png'.format(step_03_output_filename[:-4]))
plt.clf()

qualimap_command = 'qualimap bamqc --bam {}.sorted.bam -c -nr 100000 -nt {} --outdir {}/qualimapQC_before_mapq_filter --outformat HTML --java-mem-size=12G'.format(
	step_03_output_filename[:-4], cpu_count, outdir)
print('\nExecuting: {}...\n'.format(qualimap_command))
os.system(qualimap_command)

#########################################

# filter reads based on MAP5 of 5 alignment which should cleanup the data as wanted
samtools_view_command = 'samtools view -@ {} -bSq {} {}.sorted.bam > {}.sorted.MAPQ5.bam'.format(cpu_count, mapq_score, step_03_output_filename[:-4], step_03_output_filename[:-4])
print('\nExecuting: {}...\n'.format(samtools_view_command))
os.system(samtools_view_command)

# index the filtered aligned bam file for viewing in IGV
samtools_index_command = 'samtools index -@ {} {}.sorted.MAPQ5.bam'.format(cpu_count, step_03_output_filename[:-4])
print('\nExecuting: {}...\n'.format(samtools_index_command))
os.system(samtools_index_command)

# generate quality plots before mapq filtering
bam = pysam.AlignmentFile('{}.sorted.MAPQ5.bam'.format(step_03_output_filename[:-4]), 'rb')

quals = [read.mapping_quality for read in bam.fetch()]

sns.set_style('dark')
sns.set_context('paper')

plt.figure(figsize=(10,6))

ax = sns.kdeplot(quals, fill = True)
ax.set_ylabel('Density')
ax.set_xlabel('Mapping Quality Score')

plt.savefig('{}.sorted.MAPQ5.png'.format(step_03_output_filename[:-4]))
plt.clf()

qualimap_command = 'qualimap bamqc --bam {}.sorted.MAPQ5.bam -c -nr 100000 -nt {} --outdir {}/qualimapQC_after_mapq5_filter --outformat HTML --java-mem-size=12G'.format(
	step_03_output_filename[:-4], cpu_count, outdir)
print('\nExecuting: {}\n'.format(qualimap_command))
os.system(qualimap_command)

end = time.time()
print('--- Elapsed time: {} seconds. ---\n\n'.format(end - start))

###############################################################################

print('||| Step 4 - Split by "References" |||\n')
start = time.time()

step_04_output_directory = '{}/04_split_by_ref/'.format(outdir)

# Make output directory
os.system('rm -rf {}'.format(step_04_output_directory))
os.mkdir(step_04_output_directory)

# copy bam file to split directory
os.system('cp {}.sorted.MAPQ5.bam {}/04_split_by_ref'.format(step_03_output_filename[:-4], outdir))

# move to the split directory
os.chdir(step_04_output_directory)

# split bam file by reference
bam_split_command = 'bamtools split -in {}.sorted.MAPQ5.bam -reference'.format(step_03_output_filename[:-4])
print('Executing: {}...\n'.format(bam_split_command))
os.system(bam_split_command)

# delete the original bam file before proceeding
# copy bam file to split directory
os.system('rm {}.sorted.MAPQ5.bam'.format(step_03_output_filename[:-4]))

print('')

# figure out what the names of the bam files are and convert each to a fq file
for filename in listdir():
	if filename.endswith('.bam'):
		bamtofastq_command = 'bamToFastq -i {} -fq /dev/stdout | gzip > {}.fq.gz'.format(filename, filename[:-4])
		print('Executing: {}...'.format(bamtofastq_command))
		os.system(bamtofastq_command)

print('')

end = time.time()
print('--- Elapsed time: {} seconds. ---\n\n'.format(end - start))

###############################################################################

print('||| Step 5 - Extract Barcodes and Counts for Each Reference and Merge Tables |||\n')
start = time.time()

# move to the split directory
os.chdir(step_04_output_directory)

# get listing of all files and sort alphabetically such that file pairs are together
dir_files = []

for file_name in listdir():
	if file_name.endswith('fq.gz'):
		dir_files.append(file_name)
	if file_name.endswith('.fxi'):
		os.remove(file_name)

# sort the fq filenames, coz.... why not?
dir_files.sort()

# make a global extracted barcode dataframe
global_barcodes_dataframe = []

step_05_global_barcode_list_filename = '{}/05_{}_quantified_barcodes_all_references.tsv'.format(outdir, setname)
step_05_global_collapsed_barcode_list_filename = '{}/05_{}_quantified_barcodes_all_references_combined.tsv'.format(outdir, setname)
step_05_global_collapsed_barcode_chart_filename = '{}/05_{}_quantified_barcodes_all_references_combined.png'.format(outdir, setname)

global_barcode_output = open(step_05_global_barcode_list_filename, 'w')
global_collapsed_barcode_output = open(step_05_global_collapsed_barcode_list_filename, 'w')
		
dir_files_len = len(dir_files)
dir_files_chunks = np.array_split(range(dir_files_len), cpu_count)


# Run the parallelized barcode extracting function for running on Linux
barcode_extracting_manager = multiprocess.Manager()
barcode_extracting_queue = barcode_extracting_manager.Queue()    
barcode_extracting_pool = multiprocess.Pool(processes = cpu_count)

#put listener to work first
barcode_extracting_watcher = barcode_extracting_pool.apply_async(barcode_extracting_listener, (barcode_extracting_queue,))

#fire off workers
barcode_extracting_job_list = []
for core_counter in range(cpu_count):
	barcode_extracting_job = barcode_extracting_pool.apply_async(barcode_extracting_function, (list(dir_files_chunks[core_counter]), barcode_extracting_queue))
	barcode_extracting_job_list.append(barcode_extracting_job)

# collect results from the workers through the pool result queue
for barcode_extracting_job in barcode_extracting_job_list: 
	barcode_extracting_job.get()

#now we are done, kill the listener
barcode_extracting_queue.put('kill')
barcode_extracting_pool.close()
barcode_extracting_pool.join()

global_barcode_output.close()
global_collapsed_barcode_output.close()


os.system('sort -nr -k 3,3 -o {} {}'.format(step_05_global_barcode_list_filename, step_05_global_barcode_list_filename))
os.system('sort -nr -k 2,2 -o {} {}'.format(step_05_global_collapsed_barcode_list_filename, step_05_global_collapsed_barcode_list_filename))

if sys.platform == 'linux' or sys.platform == 'linux2':
	os.system("sed -i '1i Reference\tBarcodes\tBarcode Counts' {}".format(step_05_global_barcode_list_filename))
	os.system("sed -i '1i Reference\tBarcode Counts' {}".format(step_05_global_collapsed_barcode_list_filename))

if sys.platform == 'darwin':
	os.system("gsed -i '1i Reference\tBarcodes\tBarcode Counts' {}".format(step_05_global_barcode_list_filename))
	os.system("gsed -i '1i Reference\tBarcode Counts' {}".format(step_05_global_collapsed_barcode_list_filename))


global_collapsed_barcode_df = pd.read_csv(step_05_global_collapsed_barcode_list_filename, delimiter = '\t')

sns.set_style('white')
sns.set_context(None)
plt.rcParams.update({'font.size': 2.8})
bar_chart_colors = sns.color_palette('deep')
ax = global_collapsed_barcode_df.plot.barh(x = 'Reference', y = 'Barcode Counts', color = bar_chart_colors, legend = None)
plt.gca().invert_yaxis()
plt.ylabel('[reference]_[barcode]', fontsize = 5)
ax.tick_params(axis = 'y', direction = 'out', length = 2, colors = 'black', labelsize = 2.8)
# ax.yaxis.set_tick_params(labelsize = 2.8)
ax.bar_label(ax.containers[0])
ax.spines[['right', 'top', 'bottom']].set_visible(False)
ax.get_xaxis().set_visible(False)
plt.savefig(step_05_global_collapsed_barcode_chart_filename, dpi = 600, bbox_inches = 'tight')
plt.clf()

end = time.time()
print('Elapsed time: {} seconds'.format(end - start))

print('\n\n||| {} - SLICER - Program Ended |||\n\n'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())))