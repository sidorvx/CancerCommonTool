# An algorithm for ranking the potential efficacy of tyrosine kinase inhibitors for different cancer subtypes
## Description
Cancer is a serious problem of our time. It is the second most common cause of death after cardiac diseases. Drug development is a very long and laborious process.
This tool allows you to assess the effectiveness of a drug for a particular type of cancer. The algorithm is a generalized version of another one developed by Niek A. Peters, et al. (https://doi.org/10.3389/fonc.2022.969855).
## Installation
You need to download all the files of this project into one directory. It is recommended to use the library versions according to the requirements.txt file. Also in the directory you need to add a file with gene expressions of your patients and name it transposed_expressions.tsv. The list of genes associated with a certain type of cancer should be represented as a python list and inserted into line 62 of main.py 
## Usage
After performing all the steps in Installation. It is necessary to run main.py. After that you will get the final information.
## Authors
+ Sidorov Daniil. Moscow Institute of Physics and Technology. https://github.com/sidorvx
+ Andrey Kravets. BostonGene Technologies. https://github.com/laiwas
