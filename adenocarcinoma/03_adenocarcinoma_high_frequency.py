import pandas
import plotly.express as px
import plotly.io as pio

pio.renderers.default = "browser"

"""
Remove the silent coding mutations from the adenocarcinoma data.
Filter the data down to genes that are mutated in only >=5% of samples.
Create a histogram of the high frequency mutations
Search for frequency of individual genes that are significant in data
"""

adenocarcinoma = pandas.read_csv(
    "adenocarcinoma-csv/02_cosmic_adenocarcinoma_only.csv",
    header=0,
    sep=",",
    index_col=None,
    low_memory=False,
)

headers = adenocarcinoma.columns.tolist()
# print(headers)
# print(oac_data.head())


# # # 1. Remove silent coding adenocarcinoma data and create new csv
filtered_data = adenocarcinoma[
    ~adenocarcinoma[" MUTATION_AA"].str.contains("(?:p.\?)$")
    & ~adenocarcinoma[" MUTATION_DESCRIPTION"].str.contains(
        "(?:Substitution - coding silent)$"
    )
]
# print(filtered_data.value_counts().head())
filtered_data.to_csv(
    "../adenocarcinoma/adenocarcinoma-csv/03_filtered_cosmic_adenocarcinoma.csv",
    index=False,
)


# # # 2. Genes mutated in at least 5% of samples (high frequency genes)
sample_count = len(
    filtered_data[" ID_SAMPLE"].value_counts()
)  # samples = 323 # 10% = 32 genes
# print(sample_count)

five_percent = sample_count * 0.05  # 16.15
# print(five_percent)

gene_count = filtered_data["GENE_NAME"].value_counts()
# print(gene_count)
five_percent_plus = gene_count[gene_count >= five_percent].reset_index()
five_percent_plus.columns = ["GENE_NAME", "Count"]
# print(five_percent_plus)
five_percent_plus = five_percent_plus.rename(columns={"GENE_NAME": "Gene Name"})
# print(five_percent_plus)
five_percent_plus["Percentage"] = (
    (five_percent_plus["Count"] / sample_count) * 100
).to_frame()
# print(five_percent_plus)
five_percent_plus.to_csv(
    "../adenocarcinoma/adenocarcinoma-csv/03_adenocarcinoma_high_frequency_mutations(above_five_percent).csv",
    index=False,
)

# # # 3. count percentage in >5%
fivePercentData = pandas.read_csv(
    "../adenocarcinoma/adenocarcinoma-csv/03_adenocarcinoma_high_frequency_mutations(above_five_percent).csv",
    header=0,
    sep=",",
)

fivePercent = filtered_data.loc[
    filtered_data["GENE_NAME"].isin(fivePercentData["Gene Name"])
]

sample_count = len(fivePercent[" ID_SAMPLE"].value_counts())
# print(sample_count)
gene_count = fivePercent["GENE_NAME"].value_counts().to_frame("counts").reset_index()
# print(gene_count)
data = gene_count
data["Percentage"] = ((gene_count["counts"] / sample_count) * 100).to_frame()
# print(data)
data.to_csv(
    "../adenocarcinoma/adenocarcinoma-csv/03_adenocarcinoma_5%_percentages.csv",
    index=False,
)


# # # 4. Histogram of high frequency genes
fig = px.histogram(
    data,
    x=["index"],
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
    "../figures/histograms/high_frequency_adenocarcinoma.png",
    width=1920,
    height=1080,
    scale=1,
)


# # # 5. check percentage of samples that contain a specific mutation

five_percent = pandas.read_csv(
    "../adenocarcinoma/adenocarcinoma-csv/03_adenocarcinoma_high_frequency_mutations"
    "(above_five_percent).csv",
    header=0,
    sep=",",
)

filtered_adenocarcinoma = pandas.read_csv(
    "../adenocarcinoma/adenocarcinoma-csv/03_filtered_cosmic_adenocarcinoma.csv",
    header=0,
    sep=",",
)

patient = " ID_SAMPLE"
gene_name = "GENE_NAME"

five_percent_patients = filtered_adenocarcinoma.loc[
    filtered_adenocarcinoma[gene_name].isin(five_percent["Gene Name"])
]
# print(five_percent_patients)

check_full_patients = filtered_adenocarcinoma.loc[
    filtered_adenocarcinoma[patient].isin(five_percent_patients[patient])
][[patient, gene_name]]
# print(check_full_patients)
# print(check_full_patients[' ID_SAMPLE'].value_counts())

# # # # ARID1A
# mutational_load_filter = (
#     check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
#     .apply(lambda x: x[x == "ARID1A"])
#     .to_frame()
#     .reset_index()
# )
# mutational_load_filter.drop(columns="level_1", inplace=True)
# # print(mutational_load_filter)
# compare = check_full_patients.loc[
#     check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
# ]
# # print(compare)
# count = compare[" ID_SAMPLE"].value_counts()
# # print(count)
# samples = count.count()
# print(samples)
# percentage = (samples / 298) * 100
# print(percentage)
#
#
# # # # ERBB4
# mutational_load_filter = (
#     check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
#     .apply(lambda x: x[x == "ERBB4"])
#     .to_frame()
#     .reset_index()
# )
# mutational_load_filter.drop(columns="level_1", inplace=True)
# # print(mutational_load_filter)
# compare = check_full_patients.loc[
#     check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
# ]
# # print(compare)
# count = compare[" ID_SAMPLE"].value_counts()
# # print(count)
# samples = count.count()
# print(samples)
# percentage = (samples / 298) * 100
# print(percentage)
#
# # # # APC
# mutational_load_filter = (
#     check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
#     .apply(lambda x: x[x == "APC"])
#     .to_frame()
#     .reset_index()
# )
# mutational_load_filter.drop(columns="level_1", inplace=True)
# # print(mutational_load_filter)
# compare = check_full_patients.loc[
#     check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
# ]
# # print(compare)
# count = compare[" ID_SAMPLE"].value_counts()
# # print(count)
# samples = count.count()
# print(samples)
# percentage = (samples / 298) * 100
# print(percentage)
#
#
# # # # SMAD4
# mutational_load_filter = (
#     check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
#     .apply(lambda x: x[x == "SMAD4"])
#     .to_frame()
#     .reset_index()
# )
# mutational_load_filter.drop(columns="level_1", inplace=True)
# # print(mutational_load_filter)
# compare = check_full_patients.loc[
#     check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
# ]
# # print(compare)
# count = compare[" ID_SAMPLE"].value_counts()
# print(count)
# samples = count.count()
# print(samples)
# percentage = (samples / 298) * 100
# print(percentage)
#
#
# # # # LRP1B
# mutational_load_filter = (
#     check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
#     .apply(lambda x: x[x == "LRP1B"])
#     .to_frame()
#     .reset_index()
# )
# mutational_load_filter.drop(columns="level_1", inplace=True)
# # print(mutational_load_filter)
# compare = check_full_patients.loc[
#     check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
# ]
# # print(compare)
# count = compare[" ID_SAMPLE"].value_counts()
# # print(count)
# samples = count.count()
# print(samples)
# percentage = (samples / 298) * 100
# print(percentage)
#
#
# # # # FAT4
# mutational_load_filter = (
#     check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
#     .apply(lambda x: x[x == "FAT4"])
#     .to_frame()
#     .reset_index()
# )
# mutational_load_filter.drop(columns="level_1", inplace=True)
# # print(mutational_load_filter)
# compare = check_full_patients.loc[
#     check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
# ]
# # print(compare)
# count = compare[" ID_SAMPLE"].value_counts()
# # print(count)
# samples = count.count()
# print(samples)
# percentage = (samples / 298) * 100
# print(percentage)
#
#
# # # # PTPRD
# mutational_load_filter = (
#     check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
#     .apply(lambda x: x[x == "PTPRD"])
#     .to_frame()
#     .reset_index()
# )
# mutational_load_filter.drop(columns="level_1", inplace=True)
# # print(mutational_load_filter)
# compare = check_full_patients.loc[
#     check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
# ]
# # print(compare)
# count = compare[" ID_SAMPLE"].value_counts()
# # print(count)
# samples = count.count()
# print(samples)
# percentage = (samples / 298) * 100
# print(percentage)
#
#
# # # # CDKN2A
# mutational_load_filter = (
#     check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
#     .apply(lambda x: x[x == "CDKN2A"])
#     .to_frame()
#     .reset_index()
# )
# mutational_load_filter.drop(columns="level_1", inplace=True)
# # print(mutational_load_filter)
# compare = check_full_patients.loc[
#     check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
# ]
# # print(compare)
# count = compare[" ID_SAMPLE"].value_counts()
# # print(count)
# samples = count.count()
# print(samples)
# percentage = (samples / 298) * 100
# print(percentage)
#
#
# # # # CTNNA2
# mutational_load_filter = (
#     check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
#     .apply(lambda x: x[x == "CTNNA2"])
#     .to_frame()
#     .reset_index()
# )
# mutational_load_filter.drop(columns="level_1", inplace=True)
# # print(mutational_load_filter)
# compare = check_full_patients.loc[
#     check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
# ]
# # print(compare)
# count = compare[" ID_SAMPLE"].value_counts()
# # print(count)
# samples = count.count()
# print(samples)
# percentage = (samples / 298) * 100
# print(percentage)


# # # # MUC16
# mutational_load_filter = (
#     check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
#     .apply(lambda x: x[x == "MUC16"])
#     .to_frame()
#     .reset_index()
# )
# mutational_load_filter.drop(columns="level_1", inplace=True)
# # print(mutational_load_filter)
# compare = check_full_patients.loc[
#     check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
# ]
# # print(compare)
# count = compare[" ID_SAMPLE"].value_counts()
# # print(count)
# samples = count.count()
# print(samples)
# percentage = (samples / 298) * 100
# print(percentage)


# # # # TTN
# mutational_load_filter = (
#     check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
#     .apply(lambda x: x[x == "TTN"])
#     .to_frame()
#     .reset_index()
# )
# mutational_load_filter.drop(columns="level_1", inplace=True)
# # print(mutational_load_filter)
# compare = check_full_patients.loc[
#     check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
# ]
# # print(compare)
# count = compare[" ID_SAMPLE"].value_counts()
# # print(count)
# samples = count.count()
# print(samples)
# percentage = (samples / 298) * 100
# print(percentage)


# # # # TP53
# mutational_load_filter = (
#     check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
#     .apply(lambda x: x[x == "TP53"])
#     .to_frame()
#     .reset_index()
# )
# mutational_load_filter.drop(columns="level_1", inplace=True)
# print(mutational_load_filter)
# compare = check_full_patients.loc[
#     check_full_patients[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
# ]
# print(compare)
# count = compare[" ID_SAMPLE"].value_counts()
# print(count)
# samples = count.count()
# print(samples)
# percentage = (samples / 298) * 100
# print(percentage)


# # # Mutation rates of TTN

# # # has TTN
# mutational_load_main_data = check_full_patients
# muc_count = check_full_patients.groupby(' ID_SAMPLE')['GENE_NAME'].apply(lambda x: x[x == 'TTN'])
# # print(muc_count)
# mutational_load_filter = check_full_patients.groupby(' ID_SAMPLE')['GENE_NAME'].apply(lambda x: x[x == 'TTN']).to_frame().reset_index()
# mutational_load_filter.drop(columns='level_1', inplace=True)
# # print(mutational_load_filter)
#
# compare = mutational_load_main_data.loc[mutational_load_main_data[' ID_SAMPLE'].isin(mutational_load_filter[' ID_SAMPLE'])]
# # print(compare)
#
# count = compare[' ID_SAMPLE'].value_counts()
# print(count)
#
# total = count.sum()
# print(total)
#
# average = total / len(count)
# print(average)

# # # NO TTN
mutational_load_main_data = check_full_patients
ttn_count = check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"].apply(
    lambda x: x[x != "TTN"]
)
# print(muc_count)
mutational_load_filter = (
    check_full_patients.groupby(" ID_SAMPLE")["GENE_NAME"]
    .apply(lambda x: x[x == "TTN"])
    .to_frame()
    .reset_index()
)
mutational_load_filter.drop(columns="level_1", inplace=True)
# print(mutational_load_filter)
#
compare = mutational_load_main_data.loc[
    ~mutational_load_main_data[" ID_SAMPLE"].isin(mutational_load_filter[" ID_SAMPLE"])
]
# print(compare)

count = compare[" ID_SAMPLE"].value_counts()
print(count)

total = count.sum()
print(total)

average = total / len(count)
print(average)
