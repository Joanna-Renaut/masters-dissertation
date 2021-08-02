import pandas
import plotly.express as px
import plotly.io as pio

pio.renderers.default = "browser"


"""
Remove the silent coding mutations from the squamous cell data.
Filter the data down to genes that are mutated in only >=5% of samples.
Create a histogram of the high frequency mutations
Search for frequency of individual genes that are significant in data
"""

oscc_data = pandas.read_csv(
    "squamous-cell-carcinoma-csv/02_cosmic_squamous_cell_only.csv",
    header=0,
    sep=",",
    index_col=None,
    low_memory=False,
)

headers = oscc_data.columns.tolist()
# print(headers)
# print(oscc_data.head())


# # # 1. Remove silent coding OSCC data and create new csv
filtered_data = oscc_data[
    ~oscc_data[" MUTATION_AA"].str.contains("(?:p.\?)$")
    & ~oscc_data[" MUTATION_DESCRIPTION"].str.contains(
        "(?:Substitution - coding silent)$"
    )
]
# print(filtered_data.value_counts().head())
filtered_data.to_csv(
    "../squamous-cell-carcinoma/squamous-cell-carcinoma-csv/04_filtered_cosmic_squamous_cell.csv",
    index=False,
)


# # # 2. Genes mutated in at least 5% of samples (high frequency genes)
sample_count = len(filtered_data[" ID_SAMPLE"].value_counts())  # samples = 788
print(sample_count)

five_percent = sample_count * 0.05  # 39.25
print(five_percent)

gene_count = filtered_data["GENE_NAME"].value_counts()
print(gene_count)
five_percent_plus = gene_count[gene_count >= five_percent].reset_index()
five_percent_plus.columns = ["GENE_NAME", "Count"]
# print(five_percent_plus)
five_percent_plus = five_percent_plus.rename(columns={"GENE_NAME": "Gene Name"})
# print(five_percent_plus)
five_percent_plus["Percentage"] = (
    (five_percent_plus["Count"] / sample_count) * 100
).to_frame()
print(five_percent_plus)
five_percent_plus.to_csv(
    "../squamous-cell-carcinoma/squamous-cell-carcinoma-csv/"
    "04_squamous-cell-carcinoma_high_frequency_mutations(above_five_percent).csv",
    index=False,
)


# # # 3. Histogram of high frequency genes
fig = px.histogram(
    five_percent_plus,
    x=["Gene Name"],
    y="Percentage",
    title=None,
    color_discrete_sequence=px.colors.qualitative.G10,
)
fig.update_layout(
    yaxis_title_text="Percentage of samples gene is mutated in", xaxis_title_text="Gene"
)
fig.update_layout(legend_title=None, legend_font_size=25)
fig.update_layout(xaxis_title_font_size=25, yaxis_title_font_size=25)
fig.update_layout(showlegend=False)
fig.update_xaxes(tickangle=45, tickfont=dict(size=20))
fig.update_yaxes(tickfont=dict(size=20))
fig.update_yaxes(range=[0, 100])
fig.show()

fig.write_image(
    "../figures/histograms/high_frequency_squamous_cell.png",
    width=1920,
    height=1080,
    scale=1,
)


# # # 4. check percentage of samples that contain a specific mutation

five_percent = pandas.read_csv(
    "../squamous-cell-carcinoma/squamous-cell-carcinoma-csv/04_squamous-cell-carcinoma_high_"
    "frequency_mutations(above_five_percent).csv",
    header=0,
    sep=",",
)

filtered_squamous = pandas.read_csv(
    "../squamous-cell-carcinoma/squamous-cell-carcinoma-csv/04_filtered_cosmic_s"
    "quamous_cell.csv",
    header=0,
    sep=",",
)

patient = " ID_SAMPLE"
gene_name = "GENE_NAME"

five_percent_patients = filtered_squamous.loc[
    filtered_squamous[gene_name].isin(five_percent["Gene Name"])
]
# print(five_percent_patients)

check_full_patients = filtered_squamous.loc[
    filtered_squamous[patient].isin(five_percent_patients[patient])
][[patient, gene_name]]
# print(check_full_patients)

# # # FBXW7
mutational_load_filter = (
    check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
    .apply(lambda x: x[x == "FBXW7"])
    .to_frame()
    .reset_index()
)
mutational_load_filter.drop(columns="level_1", inplace=True)
# print(mutational_load_filter)
compare = check_full_patients.loc[
    check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
]
# print(compare)
count = compare[" ID_SAMPLE"].value_counts()
print(count)
samples = count.count()
print(samples)
percentage = (samples / 744) * 100
print(percentage)


# # # NOTCH1
mutational_load_filter = (
    check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
    .apply(lambda x: x[x == "NOTCH1"])
    .to_frame()
    .reset_index()
)
mutational_load_filter.drop(columns="level_1", inplace=True)
# print(mutational_load_filter)
compare = check_full_patients.loc[
    check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
]
# print(compare)
count = compare[" ID_SAMPLE"].value_counts()
print(count)
samples = count.count()
print(samples)
percentage = (samples / 744) * 100
print(percentage)


# # # KMT2D
mutational_load_filter = (
    check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
    .apply(lambda x: x[x == "KMT2D"])
    .to_frame()
    .reset_index()
)
mutational_load_filter.drop(columns="level_1", inplace=True)
# print(mutational_load_filter)
compare = check_full_patients.loc[
    check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
]
# print(compare)
count = compare[" ID_SAMPLE"].value_counts()
# print(count)
samples = count.count()
print(samples)
percentage = (samples / 744) * 100
print(percentage)


# # # PIK3CA
mutational_load_filter = (
    check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
    .apply(lambda x: x[x == "PIK3CA"])
    .to_frame()
    .reset_index()
)
mutational_load_filter.drop(columns="level_1", inplace=True)
# print(mutational_load_filter)
compare = check_full_patients.loc[
    check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
]
# print(compare)
count = compare[" ID_SAMPLE"].value_counts()
# print(count)
samples = count.count()
print(samples)
percentage = (samples / 744) * 100
print(percentage)


# # # NFE2L2
mutational_load_filter = (
    check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
    .apply(lambda x: x[x == "NFE2L2"])
    .to_frame()
    .reset_index()
)
mutational_load_filter.drop(columns="level_1", inplace=True)
# print(mutational_load_filter)
compare = check_full_patients.loc[
    check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
]
# print(compare)
count = compare[" ID_SAMPLE"].value_counts()
# print(count)
samples = count.count()
print(samples)
percentage = (samples / 744) * 100
print(percentage)


# # # FAT1
mutational_load_filter = (
    check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
    .apply(lambda x: x[x == "FAT1"])
    .to_frame()
    .reset_index()
)
mutational_load_filter.drop(columns="level_1", inplace=True)
# print(mutational_load_filter)
compare = check_full_patients.loc[
    check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
]
# print(compare)
count = compare[" ID_SAMPLE"].value_counts()
# print(count)
samples = count.count()
print(samples)
percentage = (samples / 744) * 100
print(percentage)


# # # LRP1B
mutational_load_filter = (
    check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
    .apply(lambda x: x[x == "LRP1B"])
    .to_frame()
    .reset_index()
)
mutational_load_filter.drop(columns="level_1", inplace=True)
# print(mutational_load_filter)
compare = check_full_patients.loc[
    check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
]
# print(compare)
count = compare[" ID_SAMPLE"].value_counts()
# print(count)
samples = count.count()
print(samples)
percentage = (samples / 744) * 100
print(percentage)


# # # EP300
mutational_load_filter = (
    check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
    .apply(lambda x: x[x == "EP300"])
    .to_frame()
    .reset_index()
)
mutational_load_filter.drop(columns="level_1", inplace=True)
# print(mutational_load_filter)
compare = check_full_patients.loc[
    check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
]
# print(compare)
count = compare[" ID_SAMPLE"].value_counts()
# print(count)
samples = count.count()
print(samples)
percentage = (samples / 744) * 100
print(percentage)


# # # CDKN2A
mutational_load_filter = (
    check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
    .apply(lambda x: x[x == "CDKN2A"])
    .to_frame()
    .reset_index()
)
mutational_load_filter.drop(columns="level_1", inplace=True)
# print(mutational_load_filter)
compare = check_full_patients.loc[
    check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
]
# print(compare)
count = compare[" ID_SAMPLE"].value_counts()
# print(count)
samples = count.count()
print(samples)
percentage = (samples / 744) * 100
print(percentage)


# # # FBXW7
mutational_load_filter = (
    check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
    .apply(lambda x: x[x == "FBXW7"])
    .to_frame()
    .reset_index()
)
mutational_load_filter.drop(columns="level_1", inplace=True)
# print(mutational_load_filter)
compare = check_full_patients.loc[
    check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
]
# print(compare)
count = compare[" ID_SAMPLE"].value_counts()
# print(count)
samples = count.count()
print(samples)
percentage = (samples / 744) * 100
print(percentage)


# # # MUC16
mutational_load_filter = (
    check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
    .apply(lambda x: x[x == "MUC16"])
    .to_frame()
    .reset_index()
)
mutational_load_filter.drop(columns="level_1", inplace=True)
# print(mutational_load_filter)
compare = check_full_patients.loc[
    check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
]
# print(compare)
count = compare[" ID_SAMPLE"].value_counts()
# print(count)
samples = count.count()
print(samples)
percentage = (samples / 744) * 100
print(percentage)
