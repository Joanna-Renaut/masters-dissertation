import pandas

"""
Get the drug targets for approved drugs and clinical trial drugs by comparing to drugbank.
Clean the cansar drug files to get a list of drug targets and their drugs
"""

# # # 1. compare the approved drug list to the drugbank database to look for drug targets

drugbank = pandas.read_csv("../00-database-csv/drugbank_full.csv", header=0, sep=",")
approved_drugs = pandas.read_csv(
    "../making-drug-databases/approved-drugs/12_approved_treatments_for_oesophageal_cancer.csv",
    header=0,
    sep=",",
)

gene_name = "Gene Name"
approved_drug = "Drugs_upper"
drugbank_drug = "Drug Name"

approved_targets = drugbank.loc[
    drugbank[drugbank_drug].isin(approved_drugs[approved_drug])
][[gene_name, drugbank_drug]]
# print(approved_targets)
approved_targets["Drug Status"] = "Approved"
# print(approved_targets)
approved_targets.to_csv("../00-database-csv/approved_drugs_full.csv", index=False)


# # # 2. prepare clinical trial drugs database

oc_clinical_trials = pandas.read_csv(
    "../00-database-csv/clinical_trials_.gov_OC_trials.csv", header=0, sep=","
)

# # # split the drugs column and explode each value into a new row
ct_drugs = oc_clinical_trials["Interventions"].str.split("|").explode("Interventions")
# print(ct_drugs)
ct_drugs.to_csv("making-drug-databases-csv/temporary_ct_drugs.csv", index=False)

# # # # drop anything that isn't a drug
clinical_trial_drugs = pandas.read_csv(
    "making-drug-databases-csv/temporary_ct_drugs.csv", header=0, sep=","
)
drop_non_drugs = clinical_trial_drugs[
    clinical_trial_drugs["Interventions"].astype(str).str.startswith("Drug: ")
]
# print(drop_non_drugs)

# # # # strip the Drug: from the start of the column
replace = {x.replace("Drug: ", "") for x in drop_non_drugs["Interventions"]}
remove_drug_prefix = drop_non_drugs["Interventions"].str.lstrip("Drug: ")
print(remove_drug_prefix)
remove_drug_prefix.to_csv(
    "making-drug-databases-csv/clinical_trial_drugs_only.csv", index=False
)


# # # 3. compare clinical trial drugs to drugbank to see if any different drug targets

clinical_trial_drugs = pandas.read_csv(
    "../making-drug-databases/clinical-trial-drugs/12_clinical_trial_drugs_drugs_"
    "only.csv",
    header=0,
    sep=",",
)

ct_drug = "Interventions"

ct_targets = drugbank.loc[drugbank[drugbank_drug].isin(clinical_trial_drugs[ct_drug])][
    [gene_name, drugbank_drug]
]
# print(ct_targets)
ct_targets["Drug Status"] = "In clinical trial"
# print(ct_targets)
ct_targets.to_csv("../00-database-csv/clinical_trial_drugs_full.csv", index=False)


# # # 4. get cansar repositionable drugs for SSL targets adenocarcinoma

cansar_SSL = pandas.read_csv(
    "../making-drug-databases/cansar-drugs/12_adenocarcinoma_SSL_cpat_results.csv",
    header=0,
    sep=",",
)

# # drop na values
cansar_SSL = cansar_SSL.dropna(subset=["nearest_drug_target"])
# print(cansar_SSL['nearest_drug_target'])

# # strip end
cansar_SSL = cansar_SSL["nearest_drug_target"].str.rstrip("|Whole:::Protein").to_frame()

# # split off the start
cansar_SSL = cansar_SSL["nearest_drug_target"].str.split("|", n=2).str[2].to_frame()
# print(cansar_SSL['nearest_drug_target'])

# # split off the gene name into a new column
cansar_SSL = cansar_SSL.nearest_drug_target.str.split("|", expand=True, n=1)

# # name the columns
cansar_SSL.rename(columns={0: "Gene Name", 1: "Drugs"}, inplace=True)
# print(cansar_SSL)

# # set the index
cansar_SSL.set_index("Gene Name")
# print(cansar_SSL)

# # create a dataframe from a dictionary, turn first column into list and then split second column and turn into a list
cansar_result = pandas.DataFrame(
    {
        "Gene Name": cansar_SSL["Gene Name"].values.tolist(),
        "Drugs": cansar_SSL["Drugs"].str.split(",").values.tolist(),
    }
)

# # explode the data
cansar_SSL = cansar_result.explode("Drugs")
# print(cansar_SSL)
cansar_SSL["Drug Status"] = "Repositionable cancer drug - synthetic lethality"
print(cansar_SSL)
cansar_SSL.to_csv("../00-database-csv/adenocarcinoma_SSL_cansar_full.csv", index=False)

# # # 5. get cansar respositionable drugs for oncogenes squamous cell

cansar_onco = pandas.read_csv(
    "../making-drug-databases/cansar-drugs/12_squamous_cell_oncogenes_cpat_results.csv",
    header=0,
    sep=",",
)

# # drop na values
cansar_onco = cansar_onco.dropna(subset=["nearest_drug_target"])
print(cansar_onco["nearest_drug_target"])

# # strip end
cansar_onco = (
    cansar_onco["nearest_drug_target"].str.rstrip("|Whole:::Protein").to_frame()
)

# # split off the start
cansar_onco = cansar_onco["nearest_drug_target"].str.split("|", n=2).str[2].to_frame()
print(cansar_onco["nearest_drug_target"])

# # split off the gene name into a new column
cansar_onco = cansar_onco.nearest_drug_target.str.split("|", expand=True, n=1)

# # name the columns
cansar_onco.rename(columns={0: "Gene Name", 1: "Drugs"}, inplace=True)
print(cansar_onco)

# # set the index
cansar_onco.set_index("Gene Name")
print(cansar_onco)

# # create a dataframe from a dictionary, turn first column into list and then split second column and turn into a list
cansar_result = pandas.DataFrame(
    {
        "Gene Name": cansar_onco["Gene Name"].values.tolist(),
        "Drugs": cansar_onco["Drugs"].str.split(",").values.tolist(),
    }
)

# # explode the data
cansar_onco = cansar_result.explode("Drugs")
# print(cansar_onco)
cansar_onco["Drug Status"] = "Repositionable cancer drug - GoF oncogene"
print(cansar_onco)
cansar_onco.to_csv(
    "../00-database-csv/squamous_cell_oncogenes_cansar_full.csv", index=False
)


# # # 6. get cansar repositionable drugs for SSL targets squamous cell

cansar_SSL = pandas.read_csv(
    "../making-drug-databases/cansar-drugs/12_squamous_cell_SSL_cpat_results.csv",
    header=0,
    sep=",",
)

# # drop na values
cansar_SSL = cansar_SSL.dropna(subset=["nearest_drug_target"])
# print(cansar_SSL['nearest_drug_target'])

# # strip end
cansar_SSL = cansar_SSL["nearest_drug_target"].str.rstrip("|Whole:::Protein").to_frame()

# # split off the start
cansar_SSL = cansar_SSL["nearest_drug_target"].str.split("|", n=2).str[2].to_frame()
# print(cansar_SSL['nearest_drug_target'])

# # split off the gene name into a new column
cansar_SSL = cansar_SSL.nearest_drug_target.str.split("|", expand=True, n=1)

# # name the columns
cansar_SSL.rename(columns={0: "Gene Name", 1: "Drugs"}, inplace=True)
# print(cansar_SSL)

# # set the index
cansar_SSL.set_index("Gene Name")
# print(cansar_SSL)

# # create a dataframe from a dictionary, turn first column into list and then split second column and turn into a list
cansar_result = pandas.DataFrame(
    {
        "Gene Name": cansar_SSL["Gene Name"].values.tolist(),
        "Drugs": cansar_SSL["Drugs"].str.split(",").values.tolist(),
    }
)

# # explode the data
cansar_SSL = cansar_result.explode("Drugs")
# print(cansar_SSL)
cansar_SSL["Drug Status"] = "Repositionable cancer drug - synthetic lethality"
print(cansar_SSL)
cansar_SSL.to_csv("../00-database-csv/squamous_cell_SSL_cansar_full.csv", index=False)
