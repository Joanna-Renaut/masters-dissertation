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
