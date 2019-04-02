#!/usr/bin/env python3

import os
import subprocess
import argparse
from shutil import copyfile
from sample_tsv_builder import SarekSampleFile

if __name__ == "__main__":

	parser = argparse.ArgumentParser()

	# input options
	input_group = parser.add_argument_group("Input Options")
	input_group_me = input_group.add_mutually_exclusive_group(required=True)
	input_group_me.add_argument("-d", "--samples_directory", help="Directory containing ngs 
raw samples")

	input_group.add_argument("-s", "--ngs_script", help="Script that is applied to each 
sample",
		required=True)

	# output group
	output_group = parser.add_argument_group("Output Options")
	output_group.add_argument("-o", "--output_directory", default="ngs-workspace", 
		help="Base directory where further working-directories for ngs-analysis runs will 
be set up")

	args = parser.parse_args()
	
	
	if not os.path.exists(args.ngs_script):
		print("[!] Script not found!")
		exit(1)

	with open(os.path.abspath(args.ngs_script), "r") as file:
		script = file.read()

	if not os.path.exists(args.samples_directory):
		print("[!] Invalid Samples Directory")
		exit(1)

	samples_dir = os.path.abspath(args.samples_directory)
	# transform to absolute paths
	samples = list(map(lambda x: os.path.abspath(os.path.join(args.samples_directory, x)),
			 os.listdir(samples_dir)))
	# only keep directories
	samples = list(filter(lambda x: os.path.isdir(x), samples))

	if args.output_directory:
		base_dir = os.path.abspath(args.output_directory)
	else:
		base_dir = os.path.abspath("ngs-workspace")

	# create workdir if not existent
	if not os.path.exists(base_dir):
		os.mkdir(base_dir)

	os.chdir(base_dir)
	for ngs_sample_dir in samples:
		print("Collecting samples in {}".format(ngs_sample_dir))
		run_dir = os.path.split(ngs_sample_dir)[1]

		# tsv files
		samples_tsv = SarekSampleFile(ngs_sample_dir)
		
		# only for valid ngs_sample_dir paths, a SarekSampleFile object is returned
		if samples_tsv and samples_tsv.is_valid:	
			
			# create a workdir 
			if not os.path.exists(run_dir):
				os.mkdir(run_dir)

			# change to workdir
			os.chdir(run_dir)

			# write to tsv
			sample_tsv = os.path.join(base_dir, run_dir, "samples.tsv")
			samples_tsv.write_tsv(sample_tsv)

			# place a modified version of the  script in workdir
			script_dest = os.path.join(base_dir, run_dir, 
"sarek-germline-analysis.sh")
			placeholder = "#DUMMY-STRING#" 
			with open(script_dest, "w") as file:
				# wite location information into the script
				content = script.replace(placeholder,
					os.path.join(base_dir, run_dir))
				file.write(content)

			# run ngs-bash-script in background
			cmdline = "qsub -q long {}".format(script_dest)
			print("submitting script via {}".format(cmdline))
			os.system(cmdline)

			# change back to base directory
			os.chdir(base_dir)

