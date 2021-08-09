import pandas

"""
Read in the starting cosmic data and remove duplicate genes
"""


# # # 1. read original data csv
cosmic_data = "../00-database-csv/cosmic_data.csv"
headers = pandas.read_csv(cosmic_data, index_col=None, nrows=0).columns.tolist()
# print(headers)


data = pandas.read_csv(cosmic_data, header=0, sep=",")


# # # 2. get gene Name values
print(data["GENE_NAME"].head(50))


# # # 3. remove duplicates
filtered_data = data[~data["GENE_NAME"].str.contains("(?:ENST\d*)$")]
print(filtered_data)
filtered_data.to_csv(
    "../00-database-csv/cosmic_data_duplicates_removed.csv", index=False
)


# # # 4. get gene Name values
print(filtered_data["GENE_NAME"].head(50))


# # # 5. get values of histology column
print(filtered_data[" HISTOLOGY_SUBTYPE_1"].value_counts())


# # # 6. get type of mutation
print(filtered_data[" MUTATION_DESCRIPTION"].value_counts())
