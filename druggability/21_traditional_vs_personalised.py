import pandas

"""
Counts of the treatable patients for each histology
"""


adenocarcinoma = pandas.read_csv('../druggability/druggability-csv/unstacking_adenocarcinoma_full.csv',
                           header=0, sep=',', na_values='NaN')
squamous = pandas.read_csv('../druggability/druggability-csv/unstacking_squamous_full.csv', header=0, sep=',')


# # # 1. Get total counts of treatable patients for each histology
total_count = adenocarcinoma.value_counts('Patient ID')
print(total_count)  # 219

total_count = squamous.value_counts('Patient ID')
print(total_count)  # 448


# # # 2. Get counts of treatable patients using only approved treatments
adenocarcinoma.drop(columns=['Repositionable cancer drug - synthetic lethality', 'In clinical trial'], axis=1,
                    inplace=True)
adenocarcinoma.dropna(subset=['Approved'], inplace=True)
print(adenocarcinoma)  # 62

squamous.drop(columns=[
    'Repositionable cancer drug - synthetic lethality', 'In clinical trial', 'Repositionable cancer drug - GoF oncogene'
], axis=1, inplace=True)
squamous.dropna(subset=['Approved'], inplace=True)
print(squamous)  #97


# # # 3. Get counts of treatable patients using SSL
adenocarcinoma.drop(columns=['Approved', 'In clinical trial'], axis=1, inplace=True)
adenocarcinoma.dropna(subset=['Repositionable cancer drug - synthetic lethality'], inplace=True)
print(adenocarcinoma)  # 138


