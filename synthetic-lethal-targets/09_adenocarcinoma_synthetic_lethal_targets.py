import pandas

"""
Compare pathogenic tumour suppressors from adenocarcinoma to the SLORTH database to find synthetic lethal targets.
Count how many tumour suppressors had SSL targets and how many SSL partners there are.
"""

oac_data = pandas.read_csv(
    "../tumour-suppressors-oncogenes/tumour-suppressors-oncogenes-csv/"
    "07_adenocarcinoma_pathogenic_LOF_tumour_suppressors.csv",
    header=0,
    sep=",",
)

slorth_data = pandas.read_csv("../00-database-csv/slorth_data.csv", header=0, sep=",")


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
    "../synthetic-lethal-targets/synthetic-lethal-targets-csv/09_adenocarcinoma_SSL_targets.csv",
    index=False,
)

# # 2. Count how many tumour suppressors have SSL partners
count_genes = SSL_targets["Gene Name A"].value_counts()
print(count_genes)
genes = count_genes.count()
print(genes)

# # # 3. Count how many SSL partners there are
count_targets = SSL_targets["Gene Name B"].value_counts()
print(count_targets)
target_genes = count_targets.count()
print(target_genes)
