import pandas
import plotly.express as px
import plotly.io as pio
pio.renderers.default = "browser"

"""
Create violin plots for average number of mutations per sample across dataset filtered for 5%+ mutated.
Work out mutation value counts for a student t test.
Check MUC16 mutation rates.
"""


# # # 1. create a violin plot for the average number of mutations per sample in adenocarcinoma (at least 5% mutated)
five_percent = pandas.read_csv('../adenocarcinoma/adenocarcinoma-csv/03_adenocarcinoma_high_frequency_mutations'
                               '(above_five_percent).csv', header=0, sep=',')

adenocarcinoma = pandas.read_csv('../adenocarcinoma/adenocarcinoma-csv/03_filtered_cosmic_adenocarcinoma.csv',
                                 header=0, sep=',')

patient = ' ID_SAMPLE'
gene_name = 'GENE_NAME'

five_percent_patients = adenocarcinoma.loc[adenocarcinoma[gene_name].isin(five_percent['Gene Name'])][[patient,
                                                                                                       gene_name]]
# print(five_percent_patients)

full_patients = adenocarcinoma.loc[adenocarcinoma[patient].isin(five_percent_patients[patient])][[patient]]
# print(full_patients)

oac_count = full_patients.value_counts().to_frame('OAC Counts').reset_index()
print(oac_count)

total_count = len(oac_count[' ID_SAMPLE'].value_counts())  #298
print(total_count)
mutations = oac_count['OAC Counts'].sum()
print(mutations)  # 27899
average = mutations / total_count
print(average)  # 93.6


fig = px.violin(oac_count,
                y="OAC Counts",
                box=True,
                color_discrete_sequence=px.colors.qualitative.G10)
fig.update_yaxes(title="Number of mutations")
fig.update_layout(xaxis_title_font_size=25, yaxis_title_font_size=25)
fig.update_layout(showlegend=False)
fig.update_xaxes(tickangle=45, tickfont=dict(size=20))
fig.update_yaxes(tickfont=dict(size=20))
fig.update_yaxes(range=[-100, 1700])
fig.show()
fig.write_image('violin-plots/adenocarcinoma_high_frequency_average_mutations_violin.png', width=1900, height=1080,
                scale=1)


# # # 2. Work out distribution of squamous cell mutations in samples that have at least 5%+ mutated genes
five_percent = pandas.read_csv('../squamous-cell-carcinoma/squamous-cell-carcinoma-csv/04_squamous-cell-carcinoma_high_'
                               'frequency_mutations(above_five_percent).csv', header=0, sep=',')

squamous = pandas.read_csv('../squamous-cell-carcinoma/squamous-cell-carcinoma-csv/04_filtered_cosmic_squamous_'
                           'cell.csv', header=0, sep=',')

patient = ' ID_SAMPLE'
gene_name = 'GENE_NAME'

five_percent_patients = squamous.loc[squamous[gene_name].isin(five_percent['Gene Name'])][[patient, gene_name]]
# print(five_percent_patients)

full_patients = squamous.loc[squamous[patient].isin(five_percent_patients[patient])][[patient]]
# print(full_patients)

oscc_count = full_patients.value_counts().to_frame('OSCC Counts').reset_index()
print(oscc_count)

total_count = len(oscc_count[' ID_SAMPLE'].value_counts())  #744
print(total_count)
mutations = oscc_count['OSCC Counts'].sum()
print(mutations)  # 61855
average = mutations / total_count
print(average)  # 83.13

fig = px.violin(oscc_count,
                y="OSCC Counts",
                box=True,
                color_discrete_sequence=px.colors.qualitative.Set1)
fig.update_yaxes(title="Number of mutations")
fig.update_layout(xaxis_title_font_size=25, yaxis_title_font_size=25)
fig.update_layout(showlegend=False)
fig.update_xaxes(tickangle=45, tickfont=dict(size=20))
fig.update_yaxes(tickfont=dict(size=20))
fig.update_yaxes(range=[-100, 1700])
fig.show()
fig.write_image('violin-plots/squamous_high_frequency_average_mutations_violin.png', width=1900, height=1080,
                scale=1)

# # # 3. stick the frames together for the t test data
oscc_count.drop(columns=[' ID_SAMPLE'], inplace=True)
oac_count.drop(columns=[' ID_SAMPLE'], inplace=True)
total_distribution = pandas.concat([oscc_count['OSCC Counts'], oac_count['OAC Counts']], axis=1)
# print(total_distribution)
total_distribution.to_csv('../00-database-csv/total_mutation_distribution_data.csv', index=False)


# # # 4. check average number of mutations in samples that contain a MUC16 mutation
muc_check_full_patients = adenocarcinoma.loc[adenocarcinoma[patient].isin(five_percent_patients[patient])][[patient,
                                                                                                            gene_name]]
# print(muc_check_full_patients)
muc_count = muc_check_full_patients.groupby(' ID_SAMPLE')['GENE_NAME'].apply(lambda x: x[x == 'MUC16'])
# print(muc_count)
mutational_load_filter = muc_check_full_patients.groupby(' ID_SAMPLE')['GENE_NAME'].apply(lambda x: x[
    x == 'MUC16']).to_frame().reset_index()
mutational_load_filter.drop(columns='level_1', inplace=True)
# print(mutational_load_filter)

compare = muc_check_full_patients.loc[muc_check_full_patients[' ID_SAMPLE'].isin(mutational_load_filter[' ID_SAMPLE'])]
# print(compare)

count = compare[' ID_SAMPLE'].value_counts()
print(count)

total = count.sum()
print(total)

average = total / len(count)
print(average)


# # # 5. check average number of mutations in samples that contain do not have MUC16 mutation
muc_check_full_patients = adenocarcinoma.loc[adenocarcinoma[patient].isin(five_percent_patients[patient])][[patient,
                                                                                                            gene_name]]
# print(muc_check_full_patients)

no_muc_count = muc_check_full_patients.groupby(' ID_SAMPLE')['GENE_NAME'].apply(lambda x: x[x != 'MUC16'])
# print(no_muc_count)
mutational_load_filter = muc_check_full_patients.groupby(' ID_SAMPLE')['GENE_NAME'].apply(lambda x: x[
    x == 'MUC16']).to_frame().reset_index()
mutational_load_filter.drop(columns='level_1', inplace=True)
# print(mutational_load_filter)

compare = muc_check_full_patients.loc[~muc_check_full_patients[' ID_SAMPLE'].isin(mutational_load_filter[' ID_SAMPLE'])]
# print(compare)

count = compare[' ID_SAMPLE'].value_counts()
print(count)

total = count.sum()
print(total)

average = total / len(count)
print(average)



