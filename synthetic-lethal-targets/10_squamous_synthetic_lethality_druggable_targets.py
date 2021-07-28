import pandas

"""
Compare druggable synthetic lethal targets with the squamous cell SSL targets and create a csv of gene A with their SLL
partners and druggability
"""


squamous_SSL_drugs = pandas.read_csv('../00-database-csv/squamous_cell_SSL_cansar_full.csv', header=0, sep=',')
squamous_SSL_targets = pandas.read_csv('../synthetic-lethal-targets/synthetic-lethal-targets-csv/10_squamous_cell_SSL_'
                                       'targets.csv', header=0, sep=',')

compare = squamous_SSL_targets.loc[squamous_SSL_targets['Gene Name B'].isin(squamous_SSL_drugs['Gene Name'])]
# print(compare)


compare_2 = squamous_SSL_drugs.loc[squamous_SSL_drugs['Gene Name'].isin(compare['Gene Name B'])]
# print(compare_2)

compare_2.rename(columns={'Gene Name': 'Gene Name B'}, inplace=True)

merge = pandas.merge(left=compare, right=compare_2, how='outer', on='Gene Name B')
print(merge)
merge.to_csv('../00-database-csv/squamous_SSL_druggable_targets.csv', index=False)