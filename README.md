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

```
### Python Usage
All of GATEKeeper's python functionality is encompassed within the GATEKeeper class. Upon initialization, a path to a valid fasta file must be provided, then GATEKeeper must be run manually:
```
from GATEKeeperUtils import GATEKeeper
```
