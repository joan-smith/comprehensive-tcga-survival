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

