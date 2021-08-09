import pandas

"""
Check which of the high frequency genes are loss of function.
Check which of the high frequency loss of function genes are tumour suppressors.
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

# # # 3. check for squamous cell tumour suppressors
oscc_ts = ts_onc_data.loc[ts_onc_data[gene_name].isin(oscc_high_frequency[gene_name])][
    [gene_name, role_in_cancer]
]
oscc_ts.dropna(subset=["Role in Cancer"], how="any", inplace=True)
filtered_ts = oscc_ts[oscc_ts["Role in Cancer"].str.contains("TSG")]
print(filtered_ts)

filtered_ts.to_csv(
    "../tumour-suppressors-oncogenes/tumour-suppressors-oncogenes-csv/"
    "08_squamous cell_tumour_suppressors.csv",
    index=False,
)
