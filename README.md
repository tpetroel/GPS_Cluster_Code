# GPS_Cluster_Code
## R code for identifying clusters of GPS locations to identify potential carnivore predation sites.

Cluster Identification Code, written by Albon Guillemot, Nathan Svoboda, and Tyler Petroelje
Last edited 15 August, 2014

This code was written to identify potential predation sites ("clusters") by using GPS locations from GPS collared carnivores.

Download both the .R file and example .dbf into a single folder. 
Load the R file into Program R or R Studio
Change the working directory in the .R file to the folder that contains the example .dbf

Make sure all needed packages are installed prior to runnnig this code.
Run all lines, final output is a shapefile of cluster locations.

The example data were collected with a Lotek 7000MU GPS collar. Use of GPS data collected by other types of collars will need to adjust column names according to specific collar company formats.
