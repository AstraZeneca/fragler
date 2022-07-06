![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red)

# Fragler

## Overview

This tool implements an alogrithm to recycle fragments for modular DNA assembley. For details, see Öling et al ("FRAGLER: A Fragment Recycler Application Enabling Rapid and Scalable Modular DNA Assembly").

## Dependencies

The tool uses [conda](https://docs.conda.io/en/latest/) to manage its dependencies. Follow installation instructions in either the official conda repository or any other conda flavour. We recommend [mamba](https://github.com/mamba-org/mamba) as a drop-in replacement for conda.

The tool is also dependent on acces to the GeneArt API (https://www.thermofisher.com/se/en/home/life-science/cloning/gene-synthesis/geneart-gene-synthesis.html) and requires an access key. To obtain an access key, see instructions http://assets.thermofisher.com/TFS-Assets/BID/Reference-Materials/geneart-gene-synthesis-api-presentation.pdf  and contact geneartapi@thermofisher.com 

## Installation instructions

1. Clone the repo:
```bash
git@github.com:DS-QuBi/golden-gate-cloning.git
```
2. Install dependencies
```bash
cd golden-gate-cloning
conda env create -f environment.yaml -n golden-gate-cloning
```
	
## Usage
1. Activate conda environment

```bash
conda activate golden-gate-cloning
```

2. Start the server
```bash
bash start_server.sh
```
3. Access the app through a browser using the url reported in the previous step
