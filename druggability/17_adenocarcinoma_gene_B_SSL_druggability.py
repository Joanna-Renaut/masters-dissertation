import pandas
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"

""" 
Creating a stacked bar chart to show druggability of synthetic lethal targets in adenocarcinoma
"""

SSL_targets = pandas.read_csv(
    "../synthetic-lethal-targets/synthetic-lethal-targets-csv/09_adenocarcinoma_SSL_targets.csv",
    header=0,
    sep=",",
)

approved_drugs = pandas.read_csv(
    "../00-database-csv/approved_drugs_full.csv", header=0, sep=","
)
clinical_trial_drugs = pandas.read_csv(
    "../00-database-csv/clinical_trial_drugs_full.csv", header=0, sep=","
)

repositionable_SSL_drugs = pandas.read_csv(
    "../00-database-csv/adenocarcinoma_SSL_druggable_targets.csv", header=0, sep=","
)

# # # 1. merge patients with approved drugs
SSL_targets = SSL_targets.rename(columns={"Gene Name B": "Gene Name"})
approved_drugs = approved_drugs.rename(columns={"Drug Name": "Approved Drugs"})
approved_drugs.drop(columns=["Drug Status"], axis=1, inplace=True)

merge_one = pandas.merge(
    left=SSL_targets, right=approved_drugs, how="left", on="Gene Name"
)
# print(merge_one)

# # # 2. merge first merge with clinical trial drugs
clinical_trial_drugs = clinical_trial_drugs.rename(
    columns={"Drug Name": "Clinical Trial Drugs"}
)
clinical_trial_drugs.drop(columns=["Drug Status"], axis=1, inplace=True)

merge_two = pandas.merge(
    left=merge_one, right=clinical_trial_drugs, how="left", on="Gene Name"
)
# print(merge_two)


# # # 3. merge third merge with repositionable oncogene drugs
repositionable_SSL_drugs = repositionable_SSL_drugs.rename(
    columns={"Drugs": "Repositionable SSL Drugs"}
)
repositionable_SSL_drugs.drop(
    columns=["Drug Status", "Gene Name A"], axis=1, inplace=True
)
repositionable_SSL_drugs.rename(columns={"Gene Name B": "Gene Name"}, inplace=True)

merge_three = pandas.merge(
    left=merge_two, right=repositionable_SSL_drugs, how="left", on="Gene Name"
)
# print(merge_three)

merge_three.drop_duplicates(inplace=True)

# # # get total number of SSL targets for each gene
df = merge_three[['Gene Name A', 'Gene Name']]
df.drop_duplicates(inplace=True)

count = df.groupby(['Gene Name A', 'Gene Name']).size().to_frame('Counts').reset_index()
print(count)

total = (
    count.groupby("Gene Name A")["Counts"]
    .sum('Counts')
)
print(total)

druggability_data = merge_three

# # # 4. Create columns for the counts of each of the drugs and then merge to one file

approved_drugs_count = (
    druggability_data.groupby("Gene Name")["Approved Drugs"]
    .value_counts(normalize=True)
    .to_frame("Normalised Approved Drug Counts")
    .reset_index()
)
# print(approved_drugs_count)
approved_drugs_count["Approved Drug Counts"] = (
    1 / approved_drugs_count["Normalised Approved Drug Counts"]
)
# print(approved_drugs_count)
approved_drugs_count.to_csv(
    "../figures/stacked-bars/stacked-bar-csv/adenocarcinoma_gene_B_SSL_approved_count.csv",
    index=False,
)
approved_drugs_count = pandas.read_csv(
    "../figures/stacked-bars/stacked-bar-csv/adenocarcinoma_gene_B_SSL_approved_count.csv",
    header=0,
    sep=",",
)
approved_drugs_count = approved_drugs_count.drop(
    columns=["Approved Drugs", "Normalised Approved Drug Counts"]
)
# print(approved_drugs_count)
approved_drugs_count = approved_drugs_count.groupby("Gene Name", as_index=False).max()
print(approved_drugs_count)


clinical_trial_drugs_count = (
    druggability_data.groupby("Gene Name")["Clinical Trial Drugs"]
    .value_counts(normalize=True)
    .to_frame("Normalised Clinical Trial Drug Counts")
    .reset_index()
)
clinical_trial_drugs_count["Clinical Trial Drug Counts"] = (
    1 / clinical_trial_drugs_count["Normalised Clinical Trial Drug Counts"]
)
clinical_trial_drugs_count.to_csv(
    "../figures/stacked-bars/stacked-bar-csv/adenocarcinoma_gene_B_SSL_clinical_"
    "trial_count.csv",
    index=False,
)
clinical_trial_drugs_count = pandas.read_csv(
    "../figures/stacked-bars/stacked-bar-csv/adenocarcinoma_gene_B_SSL_clinical_trial_count.csv",
    header=0,
    sep=",",
)
clinical_trial_drugs_count = clinical_trial_drugs_count.drop(
    columns=["Clinical Trial Drugs", "Normalised Clinical Trial Drug Counts"]
)
clinical_trial_drugs_count = clinical_trial_drugs_count.groupby(
    "Gene Name", as_index=False
).max()
print(clinical_trial_drugs_count)


repositionable_SSL_drugs_count = (
    druggability_data.groupby("Gene Name")["Repositionable SSL Drugs"]
    .value_counts(normalize=True)
    .to_frame("Normalised Repositionable SSL Drug Counts")
    .reset_index()
)
repositionable_SSL_drugs_count["Repositionable SSL Drug Counts"] = (
    1 / repositionable_SSL_drugs_count["Normalised Repositionable SSL Drug Counts"]
)
repositionable_SSL_drugs_count.to_csv(
    "../figures/stacked-bars/stacked-bar-csv/adenocarcinoma_repositionable_gene_B_SSL"
    "_count.csv",
    index=False,
)
repositionable_SSL_drugs_count = pandas.read_csv(
    "../figures/stacked-bars/stacked-bar-csv/adenocarcinoma_"
    "repositionable_gene_B_SSL_count.csv",
    header=0,
    sep=",",
)
repositionable_SSL_drugs_count = repositionable_SSL_drugs_count.drop(
    columns=["Repositionable SSL Drugs", "Normalised Repositionable SSL Drug Counts"]
)
repositionable_SSL_drugs_count = repositionable_SSL_drugs_count.groupby(
    "Gene Name", as_index=False
).max()
print(repositionable_SSL_drugs_count)


merge_approved_clinical = pandas.merge(
    left=approved_drugs_count,
    right=clinical_trial_drugs_count,
    how="outer",
    on="Gene Name",
)
# print(merge_approved_clinical)
final_merge = pandas.merge(
    left=merge_approved_clinical,
    right=repositionable_SSL_drugs_count,
    how="outer",
    on="Gene Name",
)
# print(final_merge)


# # # 5. Create a stacked bar chart
bar_data = final_merge

labels = bar_data["Gene Name"]

fig = go.Figure()

fig.add_trace(
    go.Bar(
        name="Targets of a repositionable drug",
        x=labels,
        y=bar_data["Repositionable SSL Drug Counts"],
        marker_color="maroon",
    )
)
fig.add_trace(
    go.Bar(
        name="Targets of an oesophageal cancer clinical trial drug",
        x=labels,
        y=bar_data["Clinical Trial Drug Counts"],
        marker_color="orange",
    )
)
fig.add_trace(
    go.Bar(
        name="Targets of a drug approved for oesophageal cancer",
        x=labels,
        y=bar_data["Approved Drug Counts"],
        marker_color="forestgreen",
    )
)

# Change the bar mode
fig.update_layout(barmode="stack")
fig.update_layout(
    legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1,
    ),
    legend_font_size=25,
)
fig.update_layout(xaxis_title_font_size=25, yaxis_title_font_size=25)
fig.update_xaxes(title_text="Genes", tickfont_size=20, tickangle=45)
fig.update_yaxes(
    title_text="Number of drugs available for the target", tickfont=dict(size=20)
)
# fig.update_yaxes(range=[0, 450])
fig.show()
fig.write_image(
    "../figures/stacked-bars/adenocarcinoma_gene_B_SSL_druggability.png",
    width=1920,
    height=1080,
    scale=1,
)
