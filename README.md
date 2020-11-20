######## Supplementary for Melanoma network manuscript

These are collection of R codes used in the manuscript. 
The collection of data are too large to be part of GitHub, and can be downloaded by the following dropbox link: 

https://www.dropbox.com/s/tixh6czclvqu3py/processed_data.zip?dl=0

Once downloaded, the data have to be unzipped under a folder named: "processed_data". 

The folder "codes" include the scripts used in the data analysis. For each code, "root.dir" (i.e. the root directory variable) has to be set properly, where this root directory holds the "processed_data" and "codes" folder. The codes includes: 

summarize_modules.R - reproduces the module ranking

MEGENA_pSKCM.R - co-expression network codes

network_validation_code_via_RRHO.R - uses RRHO algorithm to confirm enrichments of differentially expressed genes by siRNAs in the respective network neighborhoods,

Run_Primary2GTEX_MDC.R - performs module differential connectivity analysis between primary melanoma and GTEx skin data to fish out modules similarly present in normal skin transcriptome

run_CIT.R - performed causality inference test (CIT) as applied in the manuscript to fish out key methylation changes that mediated by PTPN6 and PTPRCAP. 
