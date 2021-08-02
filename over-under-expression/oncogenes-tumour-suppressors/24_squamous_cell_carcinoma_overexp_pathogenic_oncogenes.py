import pandas

"""
Check which of the high frequency genes are over expressed.
Check which of the high frequency over expressed genes are oncogenes.
Check which are pathogenic.
"""

ts_onc_data = pandas.read_csv(
    "../../00-database-csv/cancer_gene_census.csv", header=0, sep=","
)

expression_data = pandas.read_csv(
    "../../00-database-csv/gene_expression_duplicates_removed.csv",
    header=0,
    sep=",",
    na_values="0",
    low_memory=False,
)

high_frequency = pandas.read_csv(
    "../../squamous-cell-carcinoma/squamous-cell-carcinoma-csv/04_squamous-cell-carcinoma_"
    "high_frequency_mutations(above_five_percent).csv",
    header=0,
    sep=",",
)

full_data = pandas.read_csv(
    "../../squamous-cell-carcinoma/squamous-cell-carcinoma-csv/04_filtered_cosmic_squamous_cell.csv",
    header=0,
    sep=",",
)


# # # 1. Rename columns
ts_onc_data.rename(columns={"Gene Symbol": "Gene Name"}, inplace=True)
full_data.rename(columns={"GENE_NAME": "Gene Name"}, inplace=True)


# # # 2. set variables
gene_name = "Gene Name"
over_expressed = "Over"
role_in_cancer = "Role in Cancer"
fathmm = " FATHMM_PREDICTION"


# # # 3. check for adenocarcinoma gain of function
over_exp = expression_data.loc[
    expression_data[gene_name].isin(high_frequency[gene_name])
][[gene_name, over_expressed]]
print(over_exp)
over_exp.dropna(subset=["Over"], how="any", inplace=True)
print(over_exp)
count = over_exp.count()
print(count)


# # # 4. check for adenocarcinoma oncogenes
oncogenes = ts_onc_data.loc[ts_onc_data[gene_name].isin(over_exp[gene_name])][
    [gene_name, role_in_cancer]
]
# print(oac_onc)
# print(type(oac_onc))
oncogenes.dropna(subset=["Role in Cancer"], how="any", inplace=True)
# print(oac_onc)
filtered_onc = oncogenes[oncogenes["Role in Cancer"].str.contains("oncogene")]
print(filtered_onc)


# # # 5. check if oncogene is pathogenic
pathogenic = full_data.loc[full_data[gene_name].isin(filtered_onc[gene_name])][
    [gene_name, fathmm]
]
# print(pathogenic)
pathogenic.dropna(subset=[" FATHMM_PREDICTION"], how="any", inplace=True)
# print(pathogenic)
filtered_pathogenic = pathogenic[
    pathogenic[" FATHMM_PREDICTION"].str.contains("PATHOGENIC")
]
# print(filtered_pathogenic)
filtered_pathogenic.drop_duplicates(subset=["Gene Name"], inplace=True)
# print(filtered_pathogenic)


# # # 6. merge fathmm and oncogene
onc_merge = pandas.merge(
    left=filtered_pathogenic, right=filtered_onc, how="outer", on="Gene Name"
)
# print(oac_onc_merge)
onc_merge.dropna(subset=[" FATHMM_PREDICTION"], how="any", inplace=True)
onc_merge.rename(columns={" FATHMM_PREDICTION": "Fathmm Prediction"}, inplace=True)
print(onc_merge)
onc_merge.to_csv(
    "../over-under-expression-csv/24_squamous_cell_over_exp_oncogenes.csv", index=False
)
