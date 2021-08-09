import pandas
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"

""" 
Creating a stacked bar chart to show druggability of gain of function oncogenes in squamous cell carcinoma.
"""

oncogenes = pandas.read_csv(
    "../tumour-suppressors-oncogenes/tumour-suppressors-oncogenes-csv/05_adenocarcinoma_oncogenes.csv",
    header=0,
    sep=",",
)

approved_drugs = pandas.read_csv(
    "../00-database-csv/approved_drugs_full.csv", header=0, sep=","
)
clinical_trial_drugs = pandas.read_csv(
    "../00-database-csv/clinical_trial_drugs_full.csv", header=0, sep=","
)

repositionable_onco_drugs = pandas.read_csv(
    "../00-database-csv/adenocarcinoma_oncogenes_cansar_full.csv", header=0, sep=","
)

# # # 1. merge patients with approved drugs
oncogenes = oncogenes.rename(columns={"GENE_NAME": "Gene Name"})
approved_drugs = approved_drugs.rename(columns={"Drug Name": "Approved Drugs"})
approved_drugs.drop(columns=["Drug Status"], axis=1, inplace=True)

merge_one = pandas.merge(
    left=oncogenes, right=approved_drugs, how="left", on="Gene Name"
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
repositionable_onco_drugs = repositionable_onco_drugs.rename(
    columns={"Drugs": "Repositionable Oncogene Drugs"}
)
repositionable_onco_drugs.drop(columns=["Drug Status"], axis=1, inplace=True)

merge_three = pandas.merge(
    left=merge_two, right=repositionable_onco_drugs, how="left", on="Gene Name"
)

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
    "../figures/stacked-bars/stacked-bar-csv/squamous_cell_carcinoma_approved_count.csv",
    index=False,
)
approved_drugs_count = pandas.read_csv(
    "../figures/stacked-bars/stacked-bar-csv/squamous_cell_carcinoma_approved_count.csv",
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
    "../figures/stacked-bars/stacked-bar-csv/squamous_cell_carcinoma_clinical_"
    "trial_count.csv",
    index=False,
)
clinical_trial_drugs_count = pandas.read_csv(
    "../figures/stacked-bars/stacked-bar-csv/squamous_cell_carcinoma_clinical_trial_count.csv",
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


repositionable_onco_drugs_count = (
    druggability_data.groupby("Gene Name")["Repositionable Oncogene Drugs"]
    .value_counts(normalize=True)
    .to_frame("Normalised Repositionable Oncogene Drug Counts")
    .reset_index()
)
repositionable_onco_drugs_count["Repositionable Oncogene Drug Counts"] = (
    1
    / repositionable_onco_drugs_count["Normalised Repositionable Oncogene Drug Counts"]
)
repositionable_onco_drugs_count.to_csv(
    "../figures/stacked-bars/stacked-bar-csv/squamous_cell_carcinoma_repositionable_"
    "onc_count.csv",
    index=False,
)
repositionable_onco_drugs_count = pandas.read_csv(
    "../figures/stacked-bars/stacked-bar-csv/squamous_cell_carcinoma_repositionable_onc_count.csv",
    header=0,
    sep=",",
)
repositionable_onco_drugs_count = repositionable_onco_drugs_count.drop(
    columns=[
        "Repositionable Oncogene Drugs",
        "Normalised Repositionable Oncogene Drug Counts",
    ]
)
repositionable_onco_drugs_count = repositionable_onco_drugs_count.groupby(
    "Gene Name", as_index=False
).max()
print(repositionable_onco_drugs_count)


merge_approved_clinical = pandas.merge(
    left=approved_drugs_count,
    right=clinical_trial_drugs_count,
    how="outer",
    on="Gene Name",
)
# print(oac_merge_approved_clinical)
approved_clinical_onco = pandas.merge(
    left=merge_approved_clinical,
    right=repositionable_onco_drugs_count,
    how="outer",
    on="Gene Name",
)

final_merge = pandas.merge(
    left=approved_clinical_onco, right=oncogenes, how="outer", on="Gene Name"
)

print(final_merge)


# # # 5. Create a stacked bar chart
bar_data = final_merge

labels = bar_data["Gene Name"]

fig = go.Figure()

fig.add_trace(
    go.Bar(
        name="Targets of a repositionable drug",
        x=labels,
        y=bar_data["Repositionable Oncogene Drug Counts"],
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
fig.update_yaxes(range=[0, 100])
fig.show()
fig.write_image(
    "../figures/stacked-bars/adenocarcinoma_oncogene_druggability.png",
    width=1920,
    height=1080,
    scale=1,
)
