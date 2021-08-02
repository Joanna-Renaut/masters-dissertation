import pandas
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"

"""
Create csvs of druggability of each patient and create stacked bars for druggability of each patient
"""

adenocarcinoma = pandas.read_csv(
    "../personalised-regimes/personalised-regimes-csv/FINAL_adenocarcinoma_complete_"
    "patient_regimes.csv",
    header=0,
    sep=",",
)

squamous = pandas.read_csv(
    "../personalised-regimes/personalised-regimes-csv/FINAL_squamous_cell_complete_patient"
    "_regimes.csv",
    header=0,
    sep=",",
)


# # # 1. Adenocarcinoma
approved_treatment_count = (
    adenocarcinoma.groupby("Patient ID")["Drug Status"]
    .value_counts()
    .to_frame("Drug Count")
)
# print(approved_treatment_count)
approved_treatment_count.to_csv(
    "../druggability/druggability-csv/adenocarcinoma_patient_druggability.csv"
)


bar_data = pandas.read_csv(
    "../druggability/druggability-csv/adenocarcinoma_patient_druggability.csv",
    header=0,
    sep=",",
)
bar_data["Patient ID"] = bar_data["Patient ID"].astype(str)

# print(bar_data)

unstacking_adenocarcinoma = (
    bar_data.reset_index()
    .groupby(["Patient ID", "Drug Status"])["Drug Count"]
    .aggregate("first")
    .unstack()
)
print(unstacking_adenocarcinoma)
unstacking_adenocarcinoma.to_csv(
    "../druggability/druggability-csv/unstacking_adenocarcinoma_full.csv"
)

bar_data = pandas.read_csv(
    "../druggability/druggability-csv/unstacking_adenocarcinoma_full.csv",
    header=0,
    sep=",",
)
bar_data["Patient ID"] = bar_data["Patient ID"].astype(str)

labels = bar_data["Patient ID"]

fig = go.Figure()

fig.add_trace(
    go.Bar(
        name="Targets of a drug approved for another cancer",
        x=labels,
        y=bar_data["Repositionable cancer drug - synthetic lethality"],
        marker_color="maroon",
    )
)
fig.add_trace(
    go.Bar(
        name="Targets of an oesophageal cancer clinical trial drug",
        x=labels,
        y=bar_data["In clinical trial"],
        marker_color="orange",
    )
)
fig.add_trace(
    go.Bar(
        name="Targets of a drug approved for oesophageal cancer",
        x=labels,
        y=bar_data["Approved"],
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
fig.update_xaxes(
    title_text="Patients", tickfont_size=20, tickangle=45, showticklabels=False
)
fig.update_yaxes(
    title_text="Number of each type of drug available to the patient",
    tickfont=dict(size=20),
)
# fig.update_yaxes(range=[0, 160])
fig.show()
fig.write_image(
    "../figures/stacked-bars/adenocarcinoma_full_patient_druggability.png",
    width=1920,
    height=1080,
    scale=1,
)


# # # 2. Squamous cell carcinoma
approved_treatment_count = (
    squamous.groupby("Patient ID")["Drug Status"].value_counts().to_frame("Drug Count")
)
# print(approved_treatment_count)
approved_treatment_count.to_csv(
    "../druggability/druggability-csv/squamous_cell_patient_druggability.csv"
)


bar_data = pandas.read_csv(
    "../druggability/druggability-csv/squamous_cell_patient_druggability.csv",
    header=0,
    sep=",",
)
bar_data["Patient ID"] = bar_data["Patient ID"].astype(str)

# print(bar_data)

unstacking_squamous = (
    bar_data.reset_index()
    .groupby(["Patient ID", "Drug Status"])["Drug Count"]
    .aggregate("first")
    .unstack()
)
print(unstacking_squamous)
unstacking_squamous.to_csv(
    "../druggability/druggability-csv/unstacking_squamous_full.csv"
)

bar_data = pandas.read_csv(
    "../druggability/druggability-csv/unstacking_squamous_full.csv", header=0, sep=","
)
bar_data["Patient ID"] = bar_data["Patient ID"].astype(str)

labels = bar_data["Patient ID"]

fig = go.Figure()

fig.add_trace(
    go.Bar(
        name="Targets of a drug approved for another cancer",
        x=labels,
        y=bar_data["Repositionable cancer drug - synthetic lethality"],
        marker_color="maroon",
    )
)
fig.add_trace(
    go.Bar(
        name="Targets of an oesophageal cancer clinical trial drug",
        x=labels,
        y=bar_data["In clinical trial"],
        marker_color="orange",
    )
)
fig.add_trace(
    go.Bar(
        name="Targets of a drug approved for oesophageal cancer",
        x=labels,
        y=bar_data["Approved"],
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
fig.update_xaxes(
    title_text="Patients", tickfont_size=20, tickangle=45, showticklabels=False
)
fig.update_yaxes(
    title_text="Number of drugs available for the target", tickfont=dict(size=20)
)
# fig.update_yaxes(range=[0, 160])
fig.show()
fig.write_image(
    "../figures/stacked-bars/squamous_cell_full_patient_druggability.png",
    width=1920,
    height=1080,
    scale=1,
)
