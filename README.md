# scTIGER2.0 (Single-cell Temporal Inference of Gene Regulatory Networks 2.0)
Single-cell Temporal Inference of Gene Regulatory 2.0 (scTIGER2.0) Networks is an extension of scTIGER, a computational method of predicting gene regulatory networks (GRNs) that uses paired datasets of case versus control experiments. scTIGER2.0 can accept paired datasets and single case or control experiments to predict GRNs. scTIGER2.0 also selects a downsample of 200 cells from input dataset unless an override argument is provided -- a feature that is not present in the original program. After building a gene co-differential expression network, scTIGER2.0 combines cell-based pseudotiming, an attention-based convolutional neural network, and permutation-based significance testing to infer GRNs) among gene modules.

Corresponding paper for scTIGER: Dautle M, Zhang S, Chen Y. scTIGER: A Deep-Learning Method for Inferring Gene Regulatory Networks from Case versus Control scRNA-seq Datasets. International Journal of Molecular Sciences. 2023; 24(17):13339. https://doi.org/10.3390/ijms241713339

## Prerequisites 
### General
- Linux
- Python
- Pytorch
- Optional: CUDA

### Required Python Packages
- numpy
- pandas
- matplotlib
- scanpy
- scipy
- networkx
- bambi
- arviz
- leidenalg
- numba

### Installation
We recommend installing the packages using a conda environment. Information on downloading Anaconda can be found [here]([url](https://www.anaconda.com/download)). You can use the following steps to install the necessary python packages into a new environment. 
```
conda create -n scTIGER2.0 python=3.11
conda activate scTIGER2.0
git clone https://github.com/chenyongrowan/scTIGER2.0
cd scTIGER2.0
chmod +x run_scTIGER.py
unzip Data/ProstateCancer/Patient4_Benign_endothelial.zip
conda install pytorch==2.1.2 torchvision torchaudio cpuonly==2.0 -c pytorch
```

CUDA capable installation:
```
conda create -n scTIGER2.0 python=3.11
conda activate scTIGER2.0
git clone https://github.com/chenyongrowan/scTIGER2.0
cd scTIGER2.0
chmod +x run_scTIGER.py
unzip Data/ProstateCancer/Patient4_Benign_endothelial.zip
conda install pytorch==2.1.2 torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
conda install -c conda-forge bambi=0.15.0
conda install -c conda-forge scanpy anndata
```

### Input
Required data: One or two scRNA-seq dataset(s) in the format of a non-normalized CSV file with genes in the first column of the file and cells in the remaining columns. 
Required flags: 
- -goi/--geneOfInterest: One or more genes of interest. Separate multiple genes with a "+" (ex. Arc+Bdnf)
- -exp/--experimental: Path to the csv file containing the case cells. Provide file with cells as columns and genes as rows. The gene names should be the first column in the file. The file must contain at least 10 cells

Optional Flags:
- -p/--permutations: Number of permutations to run. Default 100
- -top/--numTopGenes: Number of top correlated genes selected. Default 50
- -zero/--zeroThresh: Threshold for number of 0's tolerated for a gene. Default 0.30
- -t/--timesteps: Maximum number of steps allowed between interactions. 0 implies a direct, causal interaction. 1 implies one interaction between the source and target gene. Default is 0.
- -s/--start: Starting point for scTIGER. Default 1 (Run scTIGER and and generate GRN files). 2 uses existing scTIGER output to generate GRN visualization files if you'd like to change the alpha level.
- --cuda: CUDA use on when flag included. Leave flag out if using CPU based discovery
- -o/--output: Output directory name. Default is 'scTIGER_Output'
- -a/--alpha: Alpha value for determining significant gene interactions by scTIGER discovery. Default 0.05
- -od/--override_downsample: Overrides selection of 200 cells. Default is false.
- -ctrl/--control: Path to the csv file containing control cells. Provide file with cells as columns and genes as rows. The gene names should be the first column in the file. The file must contain at least 10 cells. This flag activates scTIGER instead of scTIGER2.0 to use a co-differential expression network. 


Notes: 
If you choose to input two datasets (case and control), they should either be the same cell type and two different experimental conditions OR the same experimental condition and two different cell types. 
Sample datasets are provided in the Data folder. There are multiple sample datasets under the Data folder. 
1. The ProstateCancer folder contains datasets for one patient. The files were processed to contain only one cell type. They are also separated into benign and tumor cells.
2. The RemoteMemoryFormation folder contains preprocessed datasets. Datasets contain only neurons. Only fear conditioned (FC) and controls were selected for the Chen2 dataset.
3. The K562 folder contains a filtered dataset from the K562 cell line

### Running
The main folder of this repository contains three main files:
1. run_scTIGER.py (script to run scTIGER and scTIGER2.0)
2. scTIGER.py (definitions and functions used in run_scTIGER.py)
3. 10x_preprocess.py (process 10x files into gene expression matrix for scTIGER)

#### Preprocessing
We included a file (10x_preprocess.py) to convert 10x sequencing files to the gene expression matrix scTIGER2.0 takes in. It has a required flag (-d/--directoryPath) that takes in the path to the directory with the following 10x sequencing output files:
- features.tsv.gz
- barcodes.tsv.gz
- matrix.mtx.gz

The escript will output the gene expression matrix to input into scTIGER2.0. The command to run this script should be in the following format:
```
./10x_preprocess.py -d ./Path_to_dir_with_10x_files
```

#### scTIGER2.0 
scTIGER2.0 is set up as a single-line command using the flags defined in the input section of this page. If you are using one gene expression matrix, make sure to put the file path after the -exp flag and not the -ctrl flag. The command with required flags should be in the following format:
```
./run_scTIGER.py -goi Gene1+Gene2+Gene3 -exp ./path_to_file
```
The command with optional flags added (in any order) to adjust default parameters should be in the following format:
```
./run_scTIGER.py -goi Gene1+Gene2+Gene3 -ctrl ./path_to_file -exp ./path_to_file -p 50 -top 90 -zero 0.20 -t 1 --cuda -o NameOfDirectory_Output -s 2
```

#### Example 
To run scTIGER2.0 on the sample K562 dataset included in the Data folder, use the following command:
```
./run_scTIGER.py -goi JUNB+STAT5B+ATF5+HES1+MXD1 -exp ./Data/K562/K562.csv -p 100 -top 50 -zero 0.00 -o SampleResult_K562
```

We have also included an Example folder which allows users to run the example dataset provided to make sure their download of scTIGER2.0 is functioning. To use it, simply download the scTIGER2.0 package, make the Example folder your working directory, and run ./Example/runExample.py in your terminal. It will output the SampleResult_K562 directory. The example does not run scTIGER2.0 with CUDA. (Be sure to unzip the data in the Example folder first).

This directory contains sample output for scTIGER2.0. The main folder contains the raw counts of gene interactions detected using each entered gene of interest, a .txt file of command details with parameter values, as well as 3 directories. The sig_gene_networks directory contains filtered gene interaction files without background noise, keeping only interactions with enough counts to be deemed significant. The Graphs directory contains histograms displaying the number of interactions detected with a percentage recovery as the function of the percent recovery for each gene. The GRN_Visualization directory contains .graphml files for each gene of interest to be opened by a visualization software, such as Cytoscape.

### Visualization
scTIGER2.0 outputs .graphml files for each individual gene of interest as well as a merged overall interaction map to be uploaded to your visualization software (we used Cytoscape). We included a style file to use in Cytoscape that automatically adds directionality arrows and displays a T at the end of the arrow for downregulated interactions and an arrow for upregulated interactions. 
![alt text](https://www.dropbox.com/scl/fi/hco7ekmk9mgr12p34khhx/AML_GRN.png?rlkey=0kkeddl5m34gswefl6skdxq23&st=p0of37ya&dl=0)
