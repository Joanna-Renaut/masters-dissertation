import pandas

"""
Check which of the high frequency genes are loss of function.
Check which of the high frequency loss of function genes are tumour suppressors.
"""

ts_onc_data = pandas.read_csv('../00-database-csv/cancer_gene_census.csv', header=0, sep=',')

expression_data = pandas.read_csv('../00-database-csv/gene_expression_duplicates_removed.csv', header=0, sep=',',
                                  na_values='0', low_memory=False)

oac_high_frequency = pandas.read_csv(
    '../adenocarcinoma/adenocarcinoma-csv/03_adenocarcinoma_high_frequency_mutations(above_five_percent).csv', header=0, sep=',')

oac_full_data = pandas.read_csv('../adenocarcinoma/adenocarcinoma-csv/03_filtered_cosmic_adenocarcinoma.csv', header=0,
                                sep=',')


# # # 1. Rename columns
ts_onc_data.rename(columns={'Gene Symbol':'Gene Name'}, inplace=True)
oac_full_data.rename(columns={'GENE_NAME':'Gene Name'}, inplace=True)

# # # 2. set variables

gene_name = 'Gene Name'
loss_of_function = 'Loss'
role_in_cancer = 'Role in Cancer'
fathmm = ' FATHMM_PREDICTION'

# # # 2. check for adenocarcinoma loss of function

oac_loss = expression_data.loc[expression_data[gene_name].isin(oac_high_frequency[gene_name])][[gene_name,
                                                                                                loss_of_function]]
# print(oac_loss)
oac_loss.dropna(subset=['Loss'], how='any', inplace=True)
# print(oac_loss)
count = oac_loss.count()
print(count)

# # # 2. check for adenocarcinoma tumour suppressors

oac_ts = ts_onc_data.loc[ts_onc_data[gene_name].isin(oac_loss[gene_name])][[gene_name, role_in_cancer]]
# print(oac_ts)
# print(type(oac_ts))
oac_ts.dropna(subset=['Role in Cancer'], how='any', inplace=True)
# print(oac_ts)
filtered_onc = oac_ts[oac_ts['Role in Cancer'].str.contains('TSG')]
print(filtered_onc)

# # # 3. check if tumour suppressor is pathogenic

pathogenic = oac_full_data.loc[oac_full_data[gene_name].isin(filtered_onc[gene_name])][[gene_name, fathmm]]
# print(pathogenic)
pathogenic.dropna(subset=[' FATHMM_PREDICTION'], how='any', inplace=True)
# print(pathogenic)
filtered_pathogenic = pathogenic[pathogenic[' FATHMM_PREDICTION'].str.contains('PATHOGENIC')]
# print(filtered_pathogenic)
filtered_pathogenic.drop_duplicates(subset=['Gene Name'], inplace=True)
print(filtered_pathogenic)

# # # 4. merge fathmm and tumour suppressors

oac_ts_merge = pandas.merge(left=filtered_pathogenic, right=filtered_onc, how='outer', on='Gene Name')
# print(oac_ts_merge)
oac_ts_merge.dropna(subset=[' FATHMM_PREDICTION'], how='any', inplace=True)
oac_ts_merge.rename(columns={' FATHMM_PREDICTION': 'Fathmm Prediction'}, inplace=True)
print(oac_ts_merge)
oac_ts_merge.to_csv('../tumour-suppressors-oncogenes/tumour-suppressors-oncogenes-csv/'
                     '07_adenocarcinoma_pathogenic_LOF_tumour_suppressors.csv', index=False)





