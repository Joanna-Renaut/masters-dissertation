import pandas

"""
Compare pathogenic under expressed tumour suppressors from adenocarcinoma to the SLORTH database to find synthetic 
lethal targets.

Count how many tumour suppressors had SSL targets and how many SSL partners there are.

Compare druggable synthetic lethal targets with the adenocarcinoma SSL targets and create a csv of gene A with their SSL
partners and druggability
"""

oac_data = pandas.read_csv(
    "../over-under-expression-csv/24_adenocarcinoma_under_exp_"
    "tumour_suppressors.csv",
    header=0,
    sep=",",
)

slorth_data = pandas.read_csv(
    "../../00-database-csv/slorth_data.csv", header=0, sep=","
)


# # # 1. Check for SSL partners
gene = slorth_data["Gene Name A"]
ssl_target = slorth_data["Gene Name B"]
SSL_targets = slorth_data.loc[gene.isin(oac_data["Gene Name"])]
# print(SSL_targets)
SSL_targets = SSL_targets[SSL_targets["Synthetic Lethality Score"] >= 0.75]
# print(SSL_targets)
SSL_targets = SSL_targets[["Gene Name A", "Gene Name B", "Synthetic Lethality Score"]]
print(SSL_targets)

SSL_targets.to_csv(
    "../over-under-expression-csv/27_adenocarcinoma_under_exp_synthetic_lethal_targets.csv",
    index=False,
)


# # # 2. Count how many tumour suppressors have SSL partners
count_genes = SSL_targets["Gene Name A"].value_counts()
print(count_genes)
genes = count_genes.count()
print(genes)


# # # 3. Count how many SSL partners there are
count_targets = SSL_targets["Gene Name B"].value_counts()
print(count_targets)
target_genes = count_targets.count()
print(target_genes)


# # # 4. compare drugs to targets to get mutated target, SSL partner and drug that targets partner
adenocarcinoma_SSL_drugs = pandas.read_csv(
    "../../00-database-csv/adenocarcinoma_SSL_cansar_full.csv", header=0, sep=","
)

adenocarcinoma_SSL_targets = pandas.read_csv(
    "../over-under-expression-csv/26_adenocarcinoma_under_exp_synthetic_lethal"
    "_targets.csv",
    header=0,
    sep=",",
)

compare = adenocarcinoma_SSL_targets.loc[
    adenocarcinoma_SSL_targets["Gene Name B"].isin(
        adenocarcinoma_SSL_drugs["Gene Name"]
    )
]
# print(compare)
compare_2 = adenocarcinoma_SSL_drugs.loc[
    adenocarcinoma_SSL_drugs["Gene Name"].isin(compare["Gene Name B"])
]
# print(compare_2)
compare_2.rename(columns={"Gene Name": "Gene Name B"}, inplace=True)
merge = pandas.merge(left=compare, right=compare_2, how="outer", on="Gene Name B")
print(merge)
merge.to_csv(
    "../../00-database-csv/adenocarcinoma_under_exp_SSL_druggable_targets.csv",
    index=False,
)
