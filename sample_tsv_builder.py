#!/usr/bin/env python3

"""
Tool to build samples.tsv files for NGS Sampes from QBIC.

A samples.tsv file has the following contents:

<patient-id>\t<sex=XX/XY>\t<status=disease/control>\t<sample>\t<lane>[_<flow_cell>]\t<read-pair1>\t<read-pair2>

A Samples.tsv file is build either from a <sample>.metadata file or read from a complete sample 
directory.

Currently assumes FastQ files and pair-ended reads.
"""

import os
import argparse
import json
import re

# STATIC

# file ending of json description file
METAINFO_SUFFIX = ".metadata"


class SarekSampleFile:

    class Patient:
        def __init__(self, id, sex, diseased):
            self.id = id
            self.sex = sex
            self.diseased = "0" if diseased == "no" else "1"

    class ReadInfo:
        def __init__(self, r1, r2, lane, sample_id):
            self.sample_id = sample_id
            self.r1 = r1
            self.r2 = r2
            self.lane = lane

    def __init__(self, directory):
        self.is_valid = False
	# directoy is the absolute path to the sample-directory
        self.directory = directory
        self.patient = None
        self.reads = None
        self.metainfo_file = None
        self.read_info = []

        if os.path.exists(self.directory) and self._contains_ngs_samples(self.directory):

            self.is_valid = True

            readfiles = list(map(lambda file: os.path.join(self.directory, file), 
filter(self._is_readfile, os.listdir(self.directory))))
            self.reads = readfiles

            mb_metainfo_files = list(filter(lambda file: file.endswith(METAINFO_SUFFIX), 
os.listdir(self.directory)))
            if len(mb_metainfo_files) == 1 and mb_metainfo_files[0]:
                self.metainfo_file = os.path.join(self.directory, mb_metainfo_files[0])

            self._collect_pinfo()
            self._collect_rinfo()

        else:
            print("[!] Invalid Path or no Readfiles! ") 

    def _is_readfile(self, file):
        """
        return True if the given file is a readfile
        """
        return file.endswith("fasta") or \
               file.endswith("fasta.gz") or \
               file.endswith(".fastq") or \
               file.endswith(".fastq.gz")

    def _contains_ngs_samples(self, directory):
        """
        Return True, if the given diretory contains readfiles
        """
        if os.path.exists(directory):
            return any([self._is_readfile(f) for f in os.listdir(directory)])
        return False

    def _collect_pinfo(self):
        """
        Extract Information about Sample and Patient from dict with metainfo
        """
        if self.metainfo_file:

            with open(self.metainfo_file, "r") as file:
                metainfo = json.loads(file.read())
                readfiles = list(
                    map(lambda file: os.path.join(self.directory, file), 
                    filter(self._is_readfile, metainfo["files"])))
                self.reads = readfiles

            diseased = metainfo["sample1"]["tumor"]

            id_str = metainfo["sample1"]["id_genetics"]
            if re.search(r"_\d{2}", id_str):
                patient_id, _ = id_str.split("_")

        else:
            reads = list(filter(lambda f: self._is_readfile(f), os.listdir(self.directory)))
            patient_id = reads[0].split("_")[0]

            # alway set to diseased for now
            diseased = "1"

        # TODO: Get information about sex
        sex = "XX"

        self.patient = self.Patient(patient_id, sex, diseased)

    def _collect_rinfo(self):
        """
        Extract the read information from a list of all fastq files.
        """

        if len(self.reads) % 2 != 0:
            print("[!] Not all reads could be paired")

        r1_filter = lambda r: r.find("R1") != -1
        r2_filter = lambda r: r.find("R2") != -1
        left_reads = list(filter(r1_filter, self.reads))
        right_reads = list(filter(r2_filter, self.reads))

        # regex for sample id
        sample_id_re = re.compile(self.patient.id + r"_\d{2}_")
        # regex for lane information
        lane_info_re = re.compile(r"_L\d{3}_")

        # TODO: VERIFY!
        # flow cell identifier to solve ambiguities with lane numbering
        flow_cell_re = re.compile(r"_[A-Z]{3}\d{3}_")

        # regex for read-pair information
        readnum_re = re.compile(r"_R[12]{1}_")
        # anonymous function for stripping the read filename of its "read-code" (R1 or R2)
        replace_readnum = lambda read: re.sub(readnum_re, "", read)

        while left_reads:

            r1 = left_reads.pop()
            r1_stripped = replace_readnum(r1)

            for i in range(len(right_reads)):
                read = right_reads[i]
                if r1_stripped == replace_readnum(read):
                    r2 = right_reads.pop(i)
                    break

            if not r2:
                print("[!] No read-pair for {}".format(r1))
                exit(1)

            # extract lane information
            lane_re_match = re.search(lane_info_re, r1)
            lane = lane_re_match.group(0).replace("_", "")

            # extract flow_cell information, if exists
            flow_cell_match = re.search(flow_cell_re, r1)
            if flow_cell_match:
                flow_cell = flow_cell_match.group(0).replace("_", "")
                lane += "_" + flow_cell

            # extract sample id
            sample_id = self.patient.id

            # gather info into ReadInfo object
            read_pair = self.ReadInfo(r1, r2, lane, sample_id)

            # collect
            self.read_info.append(read_pair)

    def to_tsv(self, sample_id_filter=""):
        """
        Print the object to a tsv string
        """

        line_template = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
        output = ""

        for read_pair in self.read_info:
            if re.match(sample_id_filter, read_pair.sample_id):
                output += line_template.format(
                    self.patient.id,
                    self.patient.sex,
                    self.patient.diseased,
                    # this makes more sense
                    read_pair.sample_id,
                    read_pair.lane,
                    read_pair.r1,
                    read_pair.r2)

        return(output)

    def write_tsv(self, outfile="samples.tsv"):
        """
        Write the tsv to a file
        """
        if self.is_valid:
               with open(outfile, "w") as file:
                    file.write(self.to_tsv())
        else:
               print("[!] tsv for invalid sample directory can't be written!") 
    
    def __str__(self):
        if self.is_valid:
            return(self.to_tsv())
        else:
            return("Invalid Sample directory!")

if __name__ == "__main__":

    parser = argparse.ArgumentParser("Sarek Samples Script")

    input_group = parser.add_argument_group("Input", "Input files or directory")
    input_group_me = input_group.add_mutually_exclusive_group(required=True)
    input_group_me.add_argument("-d", "--directory", type=str,
        help="Absolute or relative path to the read directory.")

    options_group = parser.add_argument_group("Filters", "Filter options for reads")
    options_group.add_argument("-s", "--sample_id", type=str,
        help="Regex or string to describe the sample IDs that should be included in the .tsv")

    args = parser.parse_args()
    if args.directory:
        sample_tsv = SarekSampleFile(args.directory)
        
        if args.sample_id:
            print(sample_tsv.to_tsv(sample_id_filter=args.sample_id))
        else:
            print(sample_tsv.to_tsv())


