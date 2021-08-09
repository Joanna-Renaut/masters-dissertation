import pandas

"""
Counts of the treatable patients for each histology.
Run these one at a time.
"""


adenocarcinoma = pandas.read_csv(
    "../druggability/druggability-csv/unstacking_adenocarcinoma_full.csv",
    header=0,
    sep=",",
    na_values="NaN",
)
squamous = pandas.read_csv(
    "../druggability/druggability-csv/unstacking_squamous_full.csv", header=0, sep=","
)


# # # 1. Get total counts of treatable patients for each histology
total_count = adenocarcinoma.value_counts("Patient ID")
print(total_count)  # 236

total_count = squamous.value_counts("Patient ID")
print(total_count)  # 537

# # # 2. Get counts of treatable patients using only approved treatments
adenocarcinoma.drop(
    columns=["Repositionable - SSL target", "In clinical trial"],
    axis=1,
    inplace=True,
)
adenocarcinoma.dropna(subset=["Approved"], inplace=True)
print(adenocarcinoma)  # 62

# # # 3. Get counts of treatable patients using SSL
adenocarcinoma.drop(
    columns=["Approved", "In clinical trial", "Repositionable - GoF oncogene target"],
    axis=1,
    inplace=True,
)
adenocarcinoma.dropna(subset=["Repositionable - SSL target"], inplace=True)
print(adenocarcinoma)  # 139

# # # 4. Get counts of treatable patients using repositioned oncogene
adenocarcinoma.drop(
    columns=["Approved", "In clinical trial", "Repositionable - SSL target"],
    axis=1,
    inplace=True,
)
adenocarcinoma.dropna(subset=["Repositionable - GoF oncogene target"], inplace=True)
print(adenocarcinoma)  # 16

# # # 5. Get counts of treatable patients using Clinical Trial Drugs
adenocarcinoma.drop(
    columns=[
        "Approved",
        "Repositionable - SSL target",
        "Repositionable - GoF oncogene target",
    ],
    axis=1,
    inplace=True,
)
adenocarcinoma.dropna(subset=["In clinical trial"], inplace=True)
print(adenocarcinoma)  # 200

# # # 6. Get counts of treatable patients using only approved treatments
squamous.drop(
    columns=["Repositionable - SSL target", "In clinical trial"],
    axis=1,
    inplace=True,
)
squamous.dropna(subset=["Approved"], inplace=True)
print(squamous)  # 97

# # # 7. Get counts of treatable patients using SSL
squamous.drop(
    columns=["Approved", "In clinical trial", "Repositionable - GoF oncogene target"],
    axis=1,
    inplace=True,
)
squamous.dropna(subset=["Repositionable - SSL target"], inplace=True)
print(squamous)  # 389

# # # 8. Get counts of treatable patients using repositioned oncogene
squamous.drop(
    columns=["Approved", "In clinical trial", "Repositionable - SSL target"],
    axis=1,
    inplace=True,
)
squamous.dropna(subset=["Repositionable - GoF oncogene target"], inplace=True)
print(squamous.count())  # 21

# # # 9. Get counts of treatable patients using Clinical Trial Drugs
squamous.drop(
    columns=[
        "Approved",
        "Repositionable - SSL target",
        "Repositionable - GoF oncogene target",
    ],
    axis=1,
    inplace=True,
)
squamous.dropna(subset=["In clinical trial"], inplace=True)
print(squamous)  # 396
