This is a repository to hold additional code for the manuscript "Aggregation of symbionts on hosts depends on interaction type and host traits" ( https://doi-org.pitt.idm.oclc.org/10.1111/ecog.06858). 

The main analysis, results, and data are all stored in a dryad repository (https://doi.org/10.5281/zenodo.8329060) as desired by Ecography, however this repository holds additional code used for cleaning the dataset and calculating other important metrics used in the analysis. 

In this repository you will find:
Aggregation_Analysis_ExtraGraphs.R - used to generate extra graphs for visualizing the code
Aggregation_Analysis_PrGy_Trinidad.R - Initial exploration of the code before formal analysis was moved into the Rmarkdown document
FeasSet.yml - python environment to run the feasible set. Note this is in python 2 and requires specific dependencies found in this environment
AggGraphsforPub.R - code used to generate graphs that are used in the manuscript. 
Cleaning_and_Seperating_for_aggregation_Fall2019 - This is original cleaning script for preparing the original host and parasite counts into their community data to be used in the feasible set. 
Gyro_Aggregation.py - python script for running feasible set modeling to generate the median variance for each community.
Feas#Extract_f.py - python script for generating and extracting the central feasible set for comparison with observed distributions. 

If you have any issues or questions please don't hesitate to reach out!
