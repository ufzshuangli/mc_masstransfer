# mc_masstransfer
R script for 'Stabilizing microbial communities by looped mass transfer'
The folder code_lis_210406 include following materials:
1. R-scripts:
   code_lis: main analysis and visualization in the manuscript (relatively independent analysis are listed as followed)
   code_betapart: partitioning of beta-diversity
   code_growth rate_lis: calculation of netgrowth rates
   code_correlation_lis: analysis of SC vs. SC and SC vs. Para. correlation (SC:subcommunity, Para.:parameter)
   code_stability_lis: stability properties calculation
   code_sequencing_lis: analysis based on 16S rRNA amplicon sequencing data
  (each R-script can be run independently with following datasets)
 
 2.datasets:
   parameters: biotic and abiotic parameters collected during reactors running
   seq table: taxonomic tables from 16S rRNA amplicon sequencing data
   RA.txt: relative cell abundance table
   Aabundance.csv: absolute cell abundance table
   all.csv: relative cell abundance combined with all parameters
   stabilitydata.csv: relative cell abundance combined with Canberra distance
   stabilitytbl.csv: stability properties result
   
