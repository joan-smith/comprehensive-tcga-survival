# Overview
Code for analyses described in "Genome-wide identification and analysis of prognostic features in human cancers."

Main entrypoint for the computationally intensive analysis is in comprehensive_tcga_survival/main.py.
Running the code in main.py requires substantial computational power, at least 30 cores, and several hundred GB of RAM. Individual sections can be run piecemeal.


Analyses used to produce the figures are in figs_and_figdata.

# Installation

Note there are a substantial number of dependencies. The author recommends anaconda for managing python environments.

`python setup.py install`
 (or with `-e .` for local development).

# Tests
`pytest`

# Reproduction

1. Download data. `comprehensive-tcga-survival/download_data.py` is the download script. The default directory for downloads and analysis is ~/comprehensive-tcga-survival-analysis. To change this to another directory, replace line 9 with your chosen directory. Note this script downloads approximately 5GiB.

2. Run the analysis. `comprehensive-tcga-survival/mian.py` is the main analysis script. If you modified the download directory above, also modify the analysis directory on line 30. You may also want to set parallelism, on line 36. This script performs a substantial amount of analysis. Expect to need several hundred GiB of RAM, and if, running the whole script at once, many hours.

