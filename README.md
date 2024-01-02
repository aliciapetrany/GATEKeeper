# THIS PART OF THE README IS FOR THE GIRLIES AND WILL BE DELETED LATER <3

This is a temporary repository for us to work on GATEKeeper. The final respository will be on Yong's github. 

### Usage
This pipeline can be used from both the command line and within a python script. For both cases, you need to download it, cd into the GATEKeeper_complete directory, and run "make"

### Command line usage
It's kinda janky, so for it to work, you MUST be in the Gatekeeper_complete directory. Thats something we probably need to fix.  
For an extended list of params, you can check run_GATEKeeper.py, I'm too lazy to type them all out so heres the short version of the usage:

```
chmod 744 run_GATEKeeper.py
./run_GATEKeeper.py -f <path to fasta file> 
```
If you wanna run it with the time validation:
```
./run_GATEKeeper.py -f <path to fasta file> -m <path to metadata file> -x
```
The pipline also defaults with no VCF output, so if you want a vcf file run the following, where n is either 1 or 2. 1 produces a vcf for the reference against all others, while 2 produces all against all. 2 can get very unruly with a lot of sequnces. Also note that for 1, it assumes that the first sequence is the root, unless you specify otherwise with -r. 
```
./run_GATEKeeper.py -f <path to fasta file> -v <n>
```
All outputs are written to the output_file folder :)

### Python Usage
I tried to make it easy.
```
from GATEKeeperUtils import GATEKeeper
g = GATEKeeper(<path to fasta file>)
g.run()
```
You can access the outputs in the same output_folder as above, or you can access them from within the object:
```
g.mst #the final weighted minimum spanning tree
g.adj_mat #the mutation matrix
g.time_series_df #time series valdiation info
g.vcf_df #vcf info
```
To set parameters, just set the variable within the object before running. Check the GATEKeeperUtils file for all parameters im lazy and left comments there
```
from GATEKeeperUtils import GATEKeeper
g = GATEKeeper(<path to fasta file>)
g.trial_name = "yee_haw"
g.time_metadata_path = "../yee_haw.csv"
g.test_mode = True
g.run()
```
# OK I'M DRAFTING THE REAL README NOW
# GATEKeeper
Gatekeeper defines phylogenetic relationships by performing pairwise alignments across all input sequences, then contructing a minimum spanning tree based on relative sequence similarities. GATEKeeper is currently only available for linux based operating systems. 
*** INSERT WORKFLOW IMAGE HERE ***

## Installation
To install GATEKeeper, please run the following series of commands within a linux terminal:
```
git clone https://github.com/aliciapetrany/GATEKeeper.git
cd GATEKeeper
make
chmod +x run_GATEKeeper.py
chmod +x bin/GATEkeeper
```
## Requirements
The following are required for the python portion of the pipeline to function:
```
numpy>=1.26.2
pandas>=2.1.4
scipy>=1.11.4
```
The following are required for the c++ portion of the pipeline to function:
```
```
## Usages
GATEKeeper is available for use under python and command line implementations. The full pipeline can be run from both python and the command line, however, single rapid pairwise alignments can only be run from the command line at this time.
### Python Usage
All of GATEKeeper's python functionality is encompassed within the GATEKeeper class. Upon initialization, a path to a valid fasta file must be provided, then GATEKeeper must be run manually:
```
from GATEKeeperUtils import GATEKeeper
g = GATEKeeper("../Selected_1000_Genomes.fasta")
g.run()
```
GATEKeeper has optional attributes that can be adjusted with the following syntax prior to running it:
```
from GATEKeeperUtils import GATEKeeper
g = GATEKeeper("../Selected_1000_Genomes.fasta")
g.trial_name = "test_trial"
g.vcf_level = 2
g.run()
```
The exhaustive list of GATEKeeper attributes is as follows:  
  - `trial_name` [string] The suffix to be added to each output file.
  - `vcf_level` [int] The level of detail in vcf outputs. 0 = No VCF output, 1 = Root vs. all VCF output, and 2 = all vs. all VCF output. 
  - `bin_path` [string] Path to the 'bin' directory within the GATEKeeper directory. Must change if working outside of the GATEKeeper directory
  - `nperms` [int] The number of permutations run when resampling the minimum spanning tree
  - `nthreads` [int] The number of threads to run GATEKeeper calls across
  - `verbosity` [int] Level of verbosity in the output. 0 = No outputs, 1 = warnings only, 2 = warnings and messages
  - `seq_limit` [int] The number of sequences to cut the alignments off at, if neceessary. If specified, GATKeeper will only read in 'seq_limit" number of sequences, then truncate the remaining sequences in the fasta file. 
  - `root_pos` [int] Position of the root sequence in the fasta file, defaults to the first sequence. Only necessary if 'vcf_level' = 1, or if running GATEKeeper in test mode 
  - `test_mode` [boolean] If true, runs time series validation of the final minimum spanning tree. For this mode to work properly, a 'root_pos' and 'time_metadata_path' must be properly defined
  - `time_metadata_path` [string] Path to NCBI file containing information for each sequence. For more information about this file, please refer to the "Notes" section of this readme file.

## Notes
