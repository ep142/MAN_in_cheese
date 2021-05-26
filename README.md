# Microbial Association Networks in Cheese.  

This workflow is meant to illustrate a potential approach for the inference and analysis of microbial association networks from cheese microbiota data extracted from [DairyFMBN](https://data.mendeley.com/datasets/3cwf729p34/4), a repository of data on the composition of bacterial communities in dairy products.  
The workflow makes use of the [NetCoMi package](https://github.com/stefpeschel/NetCoMi), using [phyloseq](https://joey711.github.io/phyloseq/) objects extracted from DairyFMBN (only 5 are used here, but the workflow has been tested with >10 even with a computer with 8 GB RAM). The workflow has also been tested on phyloseq objects obtained from raw sequence data (see for example [SRP212264](https://www.ncbi.nlm.nih.gov/sra/?term=SRP212264)) analyzed using the [DADA2 pipeline](https://benjjneb.github.io/dada2/tutorial_1_8.html) (with SILVA as a taxonomic database).  
To test the workflow open the .Rproj file in the workflow folder, and/or create your own project in RStudio using the same folder, then open the MAN_NetCoMiFMBN_genus_example.Rmd file.  
The workflow makes use of a number of user defined functions, which can be found in the subfolder source in folder workflow. Note that a lookup table is used to (optionally) convert the taxon names of the genus _Lactobacillus_ to the new classification described by [Zheng et al. (2020)](https://www.microbiologyresearch.org/content/journal/ijsem/10.1099/ijsem.0.004107).  

If you want to try it on your own data:  

1. put the phyloseq objects in the input_data folder (you should use .rds files)

2. remove all content from the netcompare_output folder (here it is were all output files will be saved)

3. set the options in chunks 1 and 2

The workflow has been tested on two computers with 8 GB RAM running MacOS 10.14.6 with R 4.0.5 and 4.1.  

Comments and suggestions are welcome.  