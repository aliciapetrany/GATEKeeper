import numpy as np
import pandas as pd
import os
from scipy.sparse.csgraph import minimum_spanning_tree
import warnings
from multiprocessing import Process, Manager
import glob
import sys

# TODO
# arg checking
class GATEKeeper:

    #PARAMS
    seq_limit = float("inf")  ##Only read in x number of sequences from the top of the file
    root_pos = 0  # The position of the "root" sequence in the fasta file. Only necessary in test mode
    trial_name = "1"  # Tag added to end of files to indicate different trials
    nperms = 1000  # number of permuations to run the mst resampling for
    nthreads = 8 #Number of threads to run parallelized GATEKeeper calls over
    vcf_level = 2  # 2=allvall, 1 = root v all, 0 = no vcf output, note 2 does not produce standard vcf
    test_mode = False #If true, run time validation
    time_metadata_path = "" #path to ncbi file containing time info, does not need to be filtered
    gatekeeper_args = [70,15,200,25,200, False, False, False] #arguments to pass through to GATEKeeper, don't touch it <3
    bin_path = "bin/" #path to bin within GATEKeeper directory
    verbosity = 2 #2 = all warnings and prints, #1 = warnings only, #0 = no printing
    # initializes GATEKeeper object
    # params:
        # fasta_file_path = path to a fasta file containing target sequences, assuming root node is at root_pos
    def __init__(self, fasta_file_path):
        self.fasta_file_path = fasta_file_path

    # This sucks make it faster, also remove i in final version
    def load_sequences(self, file):
        i = 0
        with open(file, "r") as file:
            sequence_dict = {}
            firstkey = True
            seq = ""
            key = ""
            for line in file:
                line = line.rstrip()
                if not line.startswith(">"):
                    seq = seq + line
                    i = i + 1
                else:
                    if firstkey:
                        key = line
                        firstkey = False
                    else:
                        sequence_dict[key] = seq
                        key = line
                        seq = ""
                if i >= self.seq_limit:
                    break

            # gets last record
            sequence_dict[key] = seq
            return (sequence_dict)

    # main control function for python component of GATEKeeper
    def run(self):
        rerun = False
        self.checkargs()
        self.seq_dict = self.load_sequences(self.fasta_file_path)
        self._mkdirs()
        # prevents redundant Gatekeeper calls
        if not os.path.exists("output_files/mutation_matrix_" + self.trial_name + ".csv"):
            # multiprocessed gatekeeper calls
            self._run_GATEKeeper()
        else:
            if self.verbosity > 0:
                with warnings.catch_warnings():
                    warnings.simplefilter("always")
                    warnings.warn(
                        "Existing mutation matrix found. Skipping GATEKeeper calls. If This wasn't intended, please change trial name.")
            rerun = True
            self.adj_mat = np.loadtxt("output_files/mutation_matrix_" + self.trial_name + ".csv",
                                      delimiter=",")

        if self.vcf_level == 1 or self.vcf_level == 2:
            if rerun and os.path.exists("output_files/" + self.trial_name + ".vcf"):
                self.vcf = np.genfromtxt("output_files/" + self.trial_name + ".vcf",
                                        delimiter="\t",
                                        dtype=str)
            else:
                self._output_vcf()

        self._run_mst()

        if self.test_mode:
            self._time_series_validation()

        self._write_outputs()
        if self.verbosity >= 2:
            print("GATEKeeper ran successfully. All output files are in the 'output_files/' directory")

    #makes temp file and output file directories
    def _mkdirs(self):
        if not os.path.exists("temp_files"):
            os.makedirs("temp_files")
        if not os.path.exists("output_files"):
            os.makedirs("output_files")
    def checkargs(self):
        if not os.path.exists(self.fasta_file_path):
            sys.exit("Provided fasta file does not exist.")
        if type(self.trial_name) is not str:
            sys.exit("Invalid trial name")
        if self.root_pos < 0 or type(self.root_pos) is not int:
            sys.exit("Invalid root position specified")
        if self.nperms < 1 or type(self.nperms) is not int:
            sys.exit("Invalid number of permutations")
        if type(self.vcf_level) is not int:
            sys.exit("VCF level must be an integer")
        if self.vcf_level < 0 or self.vcf_level > 2:
            sys.exit("Invalid VCF output level. vcf_level must fall between 0-2.")
        if self.test_mode and not os.path.exists(self.time_metadata_path):
            sys.exit("Path to time-containing metadata does not exist")
        if self.seq_limit < 2 or type(self.seq_limit) is not int:
            sys.exit("Invalid seq limit")
        if self.nthreads < 1 or type(self.nthreads) is not int:
            sys.exit("Invalid number of threads")
        if not os.path.exists(self.bin_path):
            sys.exit("Invalid bin path. Please provide path to the bin folder within the GATEKeeper directory")

    #top function facilitating GATEKeeper calls. Manages all threads.
    def _run_GATEKeeper(self):
        jobs = []
        manager = Manager()
        return_dict = manager.dict()

        if self.nthreads > len(self.seq_dict.keys()):
            self.nthreads = len(self.seq_dict.keys())

        for proc in range(self.nthreads):
            # allocate even spread of work across all threads
            fasta_ids = np.array(list(self.seq_dict.keys()))
            fasta_ids = fasta_ids[np.arange(0 + proc, len(self.seq_dict.keys()), self.nthreads)]
            p = Process(target=self._single_process_GATEKeeper_calls,
                        args=[return_dict, fasta_ids, list(self.seq_dict.keys())])
            jobs.append(p)
            p.start()

        for proc in jobs:
            proc.join()

        self._organize_GATEKeeper_output(return_dict)

    #Facilitates behavior within threads
    #params:
        # return dict - a manager dictionary that is merged together after all threads finish.
        #              Each entry has an integer key representing row # in final adjacency matrix
        # fasta_ids - a list of fasta ids assigned to thread
        # keys - all keys from the sequence dict
    def _single_process_GATEKeeper_calls(self, return_dict, fasta_ids, keys):
        for fasta_id in fasta_ids:

            i = keys.index(fasta_id)
            adj = np.zeros(len(keys))

            for z in range(i + 1, len(keys)):
                self._write_fastas_for_GATEkeeper(i, z)
                command = self._gen_gatekeeper_command_string(i, z)
                os.system(command)
                adj[z] = self._get_n_mutations("temp_files/output_" + str(i) + "_" + str(z) + ".vcf")
                self._rm_temps(i, z)
            return_dict[i] = adj

    #Removes temp files after each Gatekeeper call
    #params:
        #i/z, sequence positions in sequence dict
    def _rm_temps(self, i, z):
        rm = glob.glob("temp_files/*_" + str(i) + "_" + str(z) + ".*")
        if (self.vcf_level == 1):
            if i == self.root_pos or z == self.root_pos:
                for f in rm:
                    if ".vcf" not in f:
                        os.remove(f)
            else:
                for f in rm:
                    os.remove(f)

        elif (self.vcf_level == 2):
            for f in rm:
                if ".vcf" not in f:
                    os.remove(f)
        else:
            for f in rm:
                os.remove(f)

    #Writes fastas for GATEkeeper calls
    #params:
        #i/j, sequence positions in sequence dict
    def _write_fastas_for_GATEkeeper(self, i, j):
        key1 = list(self.seq_dict.keys())[i]
        key2 = list(self.seq_dict.keys())[j]
        with open("temp_files/fa1_" + str(i) + "_" + str(j) + ".fa", "w") as file:
            file.write(key1)
            file.write("\n")
            file.write(self.seq_dict[key1])
            file.close()
        with open("temp_files/fa2_" + str(i) + "_" + str(j) + ".fa", "w") as file:
            file.write(key2)
            file.write("\n")
            file.write(self.seq_dict[key2])
            file.close()

    # retreives number of mutations in vcf file after GATEKeeper call, called in make_GATEKeeper_calls
    # params:
        # path - path to vcf file
    def _get_n_mutations(self, path):
        nmuts = 0
        with open(path, "r") as file:
            for line in file:
                if not line.startswith("#"):
                    nmuts = nmuts + 1
        return nmuts

    # Generates the string use to call GATKeeper, passes GATKeeper specific params
    #params:
        #i/j, sequence positions in sequence dict
    def _gen_gatekeeper_command_string(self, i, j):
        strl = [self.bin_path + "GATEkeeper",
                "-r", "temp_files/fa1_" + str(i) + "_" + str(j) + ".fa",
                "-q", "temp_files/fa2_" + str(i) + "_" + str(j) + ".fa",
                "-o", "temp_files/output_" + str(i) + "_" + str(j),
                "-idy", str(self.gatekeeper_args[0]),
                "-slen", str(self.gatekeeper_args[1]),
                "-alen", str(self.gatekeeper_args[2]),
                "-ind", str(self.gatekeeper_args[3]),
                "-clr", str(self.gatekeeper_args[4])]
        if (self.gatekeeper_args[5]):
            strl.append("-u")
        if (self.gatekeeper_args[6]):
            strl.append("-sen")
        if (self.gatekeeper_args[7]):
            strl.append("-dp")
        return " ".join(strl)

    # Compiles dictionary containing thread outputs into an adjacency matrix
    # params:
        # d - Manager dictionary with pooled results from each thread.
    def _organize_GATEKeeper_output(self, d):
        d = dict(d)
        d = sorted(d.items())
        d = [value for key, value in d]
        self.adj_mat = np.vstack(d)
        self.adj_mat = self.adj_mat + np.transpose(self.adj_mat)

    #Manually reads in individual vcfs and compiles them
    def _output_vcf(self):

        irange = [self.root_pos]
        if (self.vcf_level == 2):
            irange = range(len(self.seq_dict))

        vcf_df = pd.DataFrame({"CHROM": ["0"],
                               "POS": ["0"],
                               "ID": [0],
                               "REF": ["0"],
                               "ALT": ["0"],
                               "QUAL": ["0"],
                               "FILTER": ["*"],
                               "INFO": ["0"]})

        for i in irange:
            for j in range(i + 1, len(self.seq_dict.keys())):
                with open("temp_files/output_" + str(i) + "_" + str(j) + ".vcf", "r") as file:
                    for line in file:
                        if not line.startswith("#"):
                            line.rstrip()
                            line = line.split('\t')
                            temp = pd.DataFrame({"CHROM": [line[0] + ";" + list(self.seq_dict.keys())[j]],
                                                 "POS": [line[1]],
                                                 "ID": [line[2]],
                                                 "REF": [line[3]],
                                                 "ALT": [line[4]],
                                                 "QUAL": [line[5]],
                                                 "FILTER": [line[6]],
                                                 "INFO": [line[7]]})
                            vcf_df = pd.concat([vcf_df, temp])
                os.remove("temp_files/output_" + str(i) + "_" + str(j) + ".vcf")

        vcf_df = vcf_df.iloc[1:, ]
        with open("output_files/" + self.trial_name + ".vcf", 'w') as file:
            file.write("##fileformat=VCFv4.1\n")
            file.write("##reference=temp_files/fa1_0_9.fa\n")
            file.write("##source=GATEkeeper 1.0.0\n")
            file.write(
                "##INFO=<ID=TYPE,Number=1,Type=String,Description='The type of allele, either SUBSTITUTE, INSERT, or DELETE.'>\n")
            file.write("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for index, row in vcf_df.iterrows():
                outstr = '\t'.join(str(value) for value in row)
                file.write(outstr)
            file.close()
        self.vcf_df = vcf_df

    #top function for mst analysis
    def _run_mst(self):
        self.similarity_matrix = self._resample_mst()
        self.mst = minimum_spanning_tree(self.nperms * 2 - self.similarity_matrix).toarray().astype(int)
        self.mst = np.where(self.mst != 0, self.similarity_matrix, 0) / self.nperms

    #performs mst resampling
    def _resample_mst(self):
        similarity_matrix = np.zeros((self.adj_mat.shape[0], self.adj_mat.shape[1]))

        for i in range(self.nperms):
            shuffle_i = np.random.choice(np.arange(0, self.adj_mat.shape[0]),
                                         size=self.adj_mat.shape[0],
                                         replace=False)
            shuffled = minimum_spanning_tree(self.adj_mat[shuffle_i, :][:, shuffle_i])
            shuffled = shuffled.toarray().astype(int)
            x = np.argsort(shuffle_i)
            unshuffled = shuffled[x, :][:, x]
            add_array = np.where(unshuffled != 0, 1, 0)
            similarity_matrix = similarity_matrix + add_array

        return (similarity_matrix + np.transpose(similarity_matrix))

    #performs a BFS and compares time differences between child and parent nodes
    def _time_series_validation(self):
        variant_labels = np.array([s[1:] for s in list(self.seq_dict.keys())])
        self._read_in_time_metadata(variant_labels)

        time_series_df = pd.DataFrame({"Root": ["filler"],
                                       "Interactor": ["filler"],
                                       "Time_diff": [0]})

        new_roots = np.array([])
        roots = [variant_labels[self.root_pos]]
        already_traversed = np.array([])
        flag = True

        while 7 == 7:
            for root in roots:
                if root not in already_traversed:
                    root_pos = int(np.where(root == variant_labels)[0])  # get interactors with root
                    root_ints = self.mst[root_pos, :] + self.mst[:, root_pos]  # collect all interactions
                    int_poses = np.where(root_ints != 0)[0]
                    if len(int_poses) == 0:
                        continue
                    new_roots = np.concatenate([new_roots, variant_labels[int_poses]])  # save for later
                    for int_pos in int_poses:
                        if variant_labels[int_pos] in already_traversed:
                            continue
                        t_int = np.array(self.time_dict[variant_labels[int_pos]]).astype(int)
                        t_root = np.array(self.time_dict[root]).astype(int)

                        # NOTE: TIME DIFFERENCE IS CALCULATED SO THAT LATER DATES ARE POSTIVE
                        time_diff = t_int[0:2] - t_root[0:2]
                        diff_in_months = time_diff[0] * 12 + time_diff[1]
                        temp = pd.DataFrame({"Root": [root],
                                             "Interactor": [variant_labels[int_pos]],
                                             "Time_diff": [diff_in_months]})
                        time_series_df = pd.concat((time_series_df, temp))
            already_traversed = np.concatenate([already_traversed, roots])
            already_traversed = np.unique(already_traversed)
            if len(already_traversed) == len(variant_labels):  # exit condition
                break
            roots = new_roots

        self.time_series_df = time_series_df.iloc[1:]

        n_success = len(self.time_series_df[self.time_series_df["Time_diff"] > 0])
        n_zero = len(self.time_series_df[self.time_series_df["Time_diff"] == 0])
        n_fail = len(self.time_series_df[self.time_series_df["Time_diff"] < 0])
        n_total = len(self.time_series_df)

        if self.verbosity >= 2:
            print(f"Test Mode Results\nSuccesses: {n_success}\nFailures: {n_fail}\nZeros: {n_zero}\nSuccess Rate: {(n_success + n_zero)/n_total}\nFailure Rate: {n_fail/n_total}")

    #reads in the metadata file from the ncbi
    def _read_in_time_metadata(self, variant_labels):
        variant_labels = np.array([s[1:] for s in list(self.seq_dict.keys())])
        dates = pd.read_csv(self.time_metadata_path)[["Accession", "Collection_Date"]]
        dates = dates[dates["Accession"].isin(variant_labels)]
        self.time_dict = {}
        for i, date in enumerate(dates["Collection_Date"]):
            if type(date) == str:
                self.time_dict[dates["Accession"].iloc[i]] = date.split("-")

    #writes outputs to output file
    def format_cytoscape_output(self):
        variant_labels = np.array([s[1:] for s in list(self.seq_dict.keys())])
        hits = np.nonzero(self.mst)
        cyto_df = pd.DataFrame({"Interactor 1": ["0"],
                                "Interactor 2": ["0"],
                                "Confidence": [0]})
        for row, col, in zip(hits[0], hits[1]):
            temp = pd.DataFrame({"Interactor 1": [variant_labels[row]],
                                 "Interactor 2": [variant_labels[col]],
                                 "Confidence": [self.mst[row, col]]})
            cyto_df = pd.concat([cyto_df, temp])
        cyto_df = cyto_df.iloc[1:]
        self.interaction_df = cyto_df
    def _write_outputs(self):
        self.format_cytoscape_output()
        np.savetxt(fname = "output_files/mutation_matrix_" + self.trial_name + ".csv",
                   X = self.adj_mat,
                   delimiter = ",")
        np.savetxt(fname = "output_files/minimum_spanning_tree_" + self.trial_name + ".csv",
                   X = self.mst,
                   delimiter = ",")
        self.interaction_df.to_csv("output_files/cytoscape_output_" + self.trial_name + ".csv")
        if self.test_mode:
            self.time_series_df.to_csv("output_files/time_valdiation_results.csv",
                                       index = False)
