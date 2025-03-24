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
- leidenalg??
- numba??

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
conda install -c conda-forge bambi=0.15.0
conda install -c conda-forge scanpy anndata
```
