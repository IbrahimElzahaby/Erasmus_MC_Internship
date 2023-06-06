# Plasmidome investigations to address the main transfer route of anti-microbial resistance (AMR) genes in bacteria

This project aims to analyze large-scale genomic datasets (441,120 RefSeq plasmids) to annotate and discover the presence and absence of AMR genes. Regarding the large-scale datasets that we are dealing with, we firstly started trial out the pipeline workflow on a small dataset contains 261 fasta sequences. More details about the trial data can be found [here](https://github.com/IbrahimElzahaby/Erasmus_MC_Internship/tree/1f1734ae406ae491c1f0a07d1fbf759891fded06/dummy_data).

In this [directory](https://github.com/IbrahimElzahaby/Erasmus_MC_Internship/tree/1f1734ae406ae491c1f0a07d1fbf759891fded06/original_data) that contains the original data of our project, we [fetched](https://github.com/IbrahimElzahaby/Erasmus_MC_Internship/blob/1f1734ae406ae491c1f0a07d1fbf759891fded06/original_data/retrieve_fasta.py) 66,147 fasta records and this [script](https://github.com/IbrahimElzahaby/Erasmus_MC_Internship/blob/1f1734ae406ae491c1f0a07d1fbf759891fded06/original_data/progress_fasta.py) used to continue fetching the remaining sequences without fetching errors with a total fasta records of 441,120 files.


## The pipeline workflow underlying according to the following structure:

### Identify Inc types

[plasmidfinder.py](https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/README.md) is utilized to identify plasmid incompatibility types.


### Whole genome annotation

For plasmidome annotation we used Prokka tool. more information can be found through [Torsten Seemann](https://doi.org/10.1093/bioinformatics/btu153) paper


### Pan-Plasmidome investigation

For Pan-Plasmidome investigation, we utilized Panaroo pipeline. More information can be found [here](https://doi.org/10.1186/s13059-020-02090-4)


### Mass Screening of contigs for virulence genes

[ABRicate](https://github.com/tseemann/abricate) package used for mass screening of the virulence genes


### AMR gene discovery

Resistance Gene Identifier (RGI) was utilized to predict the antibiotic resistomes from our RefSeq plasmids that retrieved from NCBI. More information can be found on [CARD database](https://doi.org/10.1093/nar/gkz935)



Note: This an ongoing project. The changes will be updated accordingly.
