import pandas

"""
Check which of the high frequency genes are gain of function.
Check which of the high frequency gain of function genes are oncogenes.
"""

ts_onc_data = pandas.read_csv('../00-database-csv/cancer_gene_census.csv', header=0, sep=',')

expression_data = pandas.read_csv('../00-database-csv/gene_expression_duplicates_removed.csv', header=0, sep=',',
                                  na_values='0', low_memory=False)

oscc_high_frequency = pandas.read_csv(
    '../squamous-cell-carcinoma/squamous-cell-carcinoma-csv/04_squamous-cell-carcinoma_high_frequency_mutations(above_five_percent).csv', header=0,
    sep=',')

oscc_full_data = pandas.read_csv(
    '../squamous-cell-carcinoma/squamous-cell-carcinoma-csv/04_filtered_cosmic_squamous_cell.csv', header=0, sep=',')


# # # 1. Rename columns
ts_onc_data.rename(columns={'Gene Symbol':'Gene Name'}, inplace=True)
oscc_full_data.rename(columns={'GENE_NAME':'Gene Name'}, inplace=True)

# # # 2. set variables

gene_name = 'Gene Name'
gain_of_function = 'Gain'
role_in_cancer = 'Role in Cancer'
fathmm = ' FATHMM_PREDICTION'

# # # 2. check for squamous cell gain of function

oscc_gain = expression_data.loc[expression_data[gene_name].isin(oscc_high_frequency[gene_name])][[gene_name,
                                                                                                  gain_of_function]]
# print(oscc_gain)
oscc_gain.dropna(subset=['Gain'], how='any', inplace=True)
# print(oscc_gain)
count = oscc_gain.count()
print(count)

# # # 2. check for squamous cell oncogenes

oscc_onc = ts_onc_data.loc[ts_onc_data[gene_name].isin(oscc_gain[gene_name])][[gene_name, role_in_cancer]]
# print(oscc_onc)
# print(type(oscc_onc))
oscc_onc.dropna(subset=['Role in Cancer'], how='any', inplace=True)
# print(oscc_onc)
filtered_onc = oscc_onc[oscc_onc['Role in Cancer'].str.contains('oncogene')]
print(filtered_onc)

# # # 3. check if oncogene is pathogenic

pathogenic = oscc_full_data.loc[oscc_full_data[gene_name].isin(filtered_onc[gene_name])][[gene_name, fathmm]]
# print(pathogenic)
pathogenic.dropna(subset=[' FATHMM_PREDICTION'], how='any', inplace=True)
# print(pathogenic)
filtered_pathogenic = pathogenic[pathogenic[' FATHMM_PREDICTION'].str.contains('PATHOGENIC')]
# print(filtered_pathogenic)
filtered_pathogenic.drop_duplicates(subset=['Gene Name'], inplace=True)
print(filtered_pathogenic)

# # # 4. merge fathmm and oncogene

oscc_onc_merge = pandas.merge(left=filtered_pathogenic, right=filtered_onc, how='outer', on='Gene Name')
print(oscc_onc_merge)
oscc_onc_merge.dropna(subset=[' FATHMM_PREDICTION'], how='any', inplace=True)
oscc_onc_merge.rename(columns={' FATHMM_PREDICTION': 'Fathmm Prediction'}, inplace=True)
print(oscc_onc_merge)
oscc_onc_merge.to_csv('../tumour-suppressors-oncogenes/tumour-suppressors-oncogenes-csv/'
                     '06_squamous_cell_carcinoma_GOF_pathogenic_oncogenes.csv', index=False)





