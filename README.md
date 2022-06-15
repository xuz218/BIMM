# Uncovering the dynamic effects of DEX  treatment on lung cancer by integrating  bioinformatic inference and multiscale  modeling of scRNA-seq and proteomics data

Dexamethasone (DEX) has shown anti-cancer efficacy and anti-estrogenic activity in human non-small cell lung cancer  (NSCLC)； hence, understanding the underlying mechanisms that affect the implementation and effectiveness of lung cancer therapeutics is vital. In this study, we combine the power of  Bioinformatics and Systems Biology to uncover functional and signaling pathways of  drug treatment using bioinformatics inference and multiscale modeling of both scRNA-seq data and  proteomics data. Here, we propose a comprehensive multiscale model of tumor regulation centered on both TGF-β-induced and ERBB-amplified signaling pathways to characterize the dynamics effects of DEX therapy on lung cancer cells. We also provided predictions of  different doses to illustrate the trend and therapeutic potential of DEX treatment, and therefore can be further applied to other computational studies in tumorigenesis and oncotherapy.

## Structure of the repository

### DEX model
For Windows users: matlab -nodisplay -nosplash -nodesktop -r “run_model ; exit”<br />
For Mac users: /Applications/MATLAB_R2022a.app/bin/matlab -nodisplay -nosplash -nodesktop -r “run_model ; exit”<br />
Lines above run the model and plot the graphs with the parameter values loaded. Notice the MATLAB version should be modified to the version you have. Please make sure MATLAB R2017a or any later version is downloaded.<br />

```
┌──ga_cancer.m
├──direct.m
├──model.m
├──fitness.m
├──run_model.m
└──level_doses.m
```

- `ga_cancer.m`: contains the script for running Genetic Algorithm for optimizing parameters.
- `direct.m`: contains the script for running Direct Optimization Algorithm for optimizing parameters.
- `model.m`: contains the script for ordinary differential equations of tumorigenesis regulatory network model
- `fitness.m`: contains the script for objective functions of parameter optimizations intumorigenesis regulatory network model
- `run_model.m`: contains the script for running model with parameters loaded and plotting simulated protein profiles and expression level of signature genes
- `level_doses.m`: contains the script for plotting simulated levels of signature genes under one, two, and three doses of DEX treatment respectively


### Bioinformatics_anlaysis
```
┌──DEG_anlaysis.R
├──TF_analysis.R
├──enrichment_analysis.R
├──clinical_analysis.R
```
- `DEG_anlaysis.R`: contains the script for identifying differential expressed genes after DEX treatment.
- `TF_analysis.R`: contains the script for identifying the upstream transcriptional factors of DEG list.
- `enrichment_analysis.R`: contains the script for functional enrichment analysis 
- `clinical_analysis.R`: contains the script for analyzing the clinical significance

## Software Requirement

MATLAB; R

## Contact

If you have any further questions or suggestions, please contact [chenm@wfu.edu](mailto:chenm@wfu.edu) or [qsong@wakehealth.edu](mailto:qsong@wakehealth.edu)
