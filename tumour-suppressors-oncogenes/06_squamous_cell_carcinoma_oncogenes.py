import pandas

"""
Check which of the high frequency genes are gain of function.
Check which of the high frequency gain of function genes are oncogenes.
"""

ts_onc_data = pandas.read_csv(
    "../00-database-csv/cancer_gene_census.csv", header=0, sep=","
)

oscc_high_frequency = pandas.read_csv(
    "../squamous-cell-carcinoma/squamous-cell-carcinoma-csv/04_squamous-cell-carcinoma_high_frequency_mutations(above_five_percent).csv",
    header=0,
    sep=",",
)


# # # 1. Rename columns
ts_onc_data.rename(columns={"Gene Symbol": "Gene Name"}, inplace=True)

# # # 2. set variables
gene_name = "Gene Name"
role_in_cancer = "Role in Cancer"

# # # 3. check for squamous cell oncogenes
oscc_onc = ts_onc_data.loc[ts_onc_data[gene_name].isin(oscc_high_frequency[gene_name])][
    [gene_name, role_in_cancer]
]
# print(oscc_onc)
# print(type(oscc_onc))
oscc_onc.dropna(subset=["Role in Cancer"], how="any", inplace=True)
# print(oscc_onc)
filtered_onc = oscc_onc[oscc_onc["Role in Cancer"].str.contains("oncogene")]
print(filtered_onc)

filtered_onc.to_csv(
    "../tumour-suppressors-oncogenes/tumour-suppressors-oncogenes-csv/"
    "06_squamous_cell_carcinoma_oncogenes.csv",
    index=False,
)
