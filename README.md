# Microbiome-Analyses
This repository contains code for statistical analyses of microbiome data.

## Background and Context
The code contained in "Longitudinal-Microbiome-Analysis.R" was used to carry out longitudinal statistical analyses with limma in R to identify differences in gut microbial communities between mice with a mutation in Mecp2 and wild-type controls. Mutations in Mecp2 cause Rett syndrome, a neurodevelopmental disorder that also features gastrointestinal and metabolic dysfunction. In this experiment, fecal samples were collected from mutant and wild-type mice on a weekly basis from 5-19 weeks of age and then sequenced via 16S sequencing to characterize the gut microbiome longitudinally across disease course. The findings are published [in Communications Biology.](https://rdcu.be/cDkCI) 

## Aproach
Analyses are based on Amplicon Sequence Variants (ASVs) obtained using phyloseq. All analyses were stratified by sex. Limma voom was used to normalize the data, which was then fit to a linear model with lmFit() followed by emperical Bayes (eBayes()) statistics for estimating differential abundance of ASVs. 

The first analysis in the script is a cross-sectional analysis with stratification by age. The second analysis is a longitudinal analysis which uses the duplicateCorrelation() approach to account for longitudinal data with repeated measures. 

### Additional Information
Microbiome data from this experiment was also integrated with metabolomics data. This analysis can be found in the [Metabolomics Data Integration Repository.](https://github.com/karineier/Metabolomics-Data-Integration) Please feel free to contact me if you have any questions and if you are interested in implementing these approaches in your work!
