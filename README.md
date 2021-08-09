# Computational Analysis of Oesophageal Cancer Data to Predict Personalised Drug Regimens


## Databases
This project used:
- Oesophageal cancer mutation data (genome screens) from v92 of the COSMIC (Catalogue of Somatic 
Mutations in Cancer) database https://cancer.sanger.ac.uk/cosmic
- Oncogene and tumour suppressor data from the Cancer Gene Census https://cancer.sanger.ac.uk/census
- Synthetic lethal target information from the SLORTH database http://slorth.biochem.sussex.ac.uk/welcome/index
- Drug information from the DrugBank database https://go.drugbank.com/
- Drug information from canSAR https://cansarblack.icr.ac.uk/
- Clinical trial information from https://clinicaltrials.gov/
- Approved drug information from https://www.cancer.gov/about-cancer/treatment/drugs/esophageal

database csv files and xml files used in this project can be found at 
https://drive.google.com/drive/folders/1hvZ4EtuhUdydkiHwRYL7EIwIk8WpWk3L?usp=sharing

## Usage

1. Run **01_cosmic_data.py** to remove duplicate genes from the data and to get information about the histologies
2. Run **02_histologies.py** to separate data into csv files based on histologies (adenocarcinoma, squamous cell carcinoma 
and (NS) non-specified in this case, to get counts of the number of samples in each histology and create a pie chart.
3. Run **03_adenocarcinoma_high_frequency.py** to remove silent coding mutations, calculate genes that are mutated in
\>= 5% of patient samples and create a csv. Create another csv containing the percentage of samples the gene is mutated 
in and create a histogram of this data.
4. Run **04_squamous_cell_high_frequency.py** to remove silent coding mutations, calculate genes that are mutated in
\>= 5% of patient samples and create a csv. Create another csv containing the percentage of samples the gene is mutated 
in and create a histogram of this data.
5. run **05_adenocarcinoma_oncogenes.py** to compare the high frequency genes found in file 03 to cancer gene census 
data and find which genes are oncogenes. Create a csv of the data.
6. run **06_squamous_cell_oncogenes.py** to compare the high frequency genes found in file 04 to cancer gene census 
data and find which genes are oncogenes. Create a csv of the data.
7. run **07_adenocarcinoma_tumour_suppressors.py** to compare the high frequency genes found in file 03 to cancer gene 
census data and find which genes are tumour suppressors. Create a csv of the data.
8. run **08_squamous_cell_tumour_suppressors.py** to compare the high frequency genes found in file 04 to cancer gene 
census data and find which genes are tumour suppressors. Create a csv of the data.
9. run **09a_adenocarcinoma_synthetic_lethal_targets.py** to compare tumour suppressor data from file 07 to SSL data from 
the SLORTH database to find synthetic lethal partners and then count how many partners there are.
10. run **09b_adenocarcinoma_synthetic_lethal_druggable_targets.py** to compare the targets found in file 09a to the canSAR
database data and find druggable synthetic lethal targets.
11. run **10a_squamous_cell_synthetic_lethal_targets.py** to compare tumour suppressor data from file 08 to SSL data from 
the SLORTH database to find synthetic lethal partners and then count how many partners there are.
12. run **10b_squamous_cell_synthetic_lethal_druggable_targets.py** to compare the targets found in file 10a to the canSAR
database data and find druggable synthetic lethal targets.
13. run **11_drugbank.py** to parse the DrugBank xml file and create a csv of drug data
14. run **12_create_drug_databases.py** to create databases of drug/target information by comparing the list of approved
drugs to drugbank to get target information and also manipulated the clinical trial data to create a new database of 
drugs/targets. Finally, created a new database of drugs/targets from the canSAR cpat results.
15. run **13_adenocarcinoma_personalised_regimes.py** to compare all the databases made in file 12 to the high frequency
genes from part 3 and create a new database of predicted personalised regimes for adenocarcinoma.
16. run **14_squamous_cell_personalised_regimes.py** to compare all the databases made in file 12 to the high frequency
genes from part 4 and create a new database of predicted personalised regimes for squamous cell carcinoma.
17. run **15_adenocarcinoma_oncogene_druggability.py** to create a stacked bar chart based on counts of druggability of 
oncogene targets in adenocarcinoma.
18. run **16_squamous_cell_oncogene_druggability.py** to create a stacked bar chart based on counts of druggability of 
oncogene targets in squamous cell carcinoma.
19. run **17_adenocarcinoma_gene_B_SSL_druggability.py** to create a stacked bar chart based on counts of druggability 
of synthetic lethal partners of tumour suppressors in adenocarcinoma.
20. run **18_squamous_cell_gene_B_SSL_druggability.py** to create a stacked bar chart based on counts of druggability 
of synthetic lethal partners of tumour suppressors in adenocarcinoma.
21. run **19_full_patient_druggability.py** to create stacked bar charts based on the patient regimes made in files 13
and 14.
22. run **20_traditional_vs_personalised.py** to get counts of the number of patients treatable with each method of 
treatment e.g. approved, clinical trial and repositioned drugs
23. run **21_violin_plots_mutation_distribution.py** to calculate the average number of mutations in each histology
and create violin plots to visualise the data.

