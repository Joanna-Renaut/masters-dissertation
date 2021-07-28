import pandas

"""
Compare druggable synthetic lethal targets with the adenocarcinoma SSL targets and create a csv of gene A with their SLL
partners and druggability
"""


adenocarcinoma_SSL_drugs = pandas.read_csv('../00-database-csv/adenocarcinoma_SSL_cansar_full.csv', header=0, sep=',')
adenocarcinoma_SSL_targets = pandas.read_csv('../synthetic-lethal-targets/synthetic-lethal-targets-csv/09_adenocarcinoma_SSL'
                                       '_targets.csv', header=0, sep=',')

compare = adenocarcinoma_SSL_targets.loc[adenocarcinoma_SSL_targets['Gene Name B'].isin(adenocarcinoma_SSL_drugs['Gene Name'])]
# print(compare)


compare_2 = adenocarcinoma_SSL_drugs.loc[adenocarcinoma_SSL_drugs['Gene Name'].isin(compare['Gene Name B'])]
# print(compare_2)

compare_2.rename(columns={'Gene Name': 'Gene Name B'}, inplace=True)

merge = pandas.merge(left=compare, right=compare_2, how='outer', on='Gene Name B')
print(merge)
merge.to_csv('../00-database-csv/adenocarcinoma_SSL_druggable_targets.csv', index=False)
