import pandas
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio


pio.renderers.default = "browser"

"""
Read the data with duplicates removed and split into a csv for each histology - adenocarcinoma and squamous cell.
Filter out non-specified samples.
Count number of samples and create a pie chart of the distribution of histologies.
"""



# # # 1. Read filtered data
filtered_duplicates = pandas.read_csv('../00-database-csv/cosmic_data_duplicates_removed.csv', header=0, sep=',',
                                      low_memory=False)
print(filtered_duplicates)


# # # 2. create CSV that is filtered for adenocarcinoma (OAC)
oac_filtered_data = filtered_duplicates[filtered_duplicates[' HISTOLOGY_SUBTYPE_1'].str.contains('(?:adenocarcinoma)$')]
print(oac_filtered_data)
oac_filtered_data.to_csv('../adenocarcinoma/adenocarcinoma-csv/02_cosmic_adenocarcinoma_only.csv', index=False)


# # # 3. Create CSV that is filtered for squamous cell carcinoma (OSCC)
oscc_filtered_data = filtered_duplicates[
    filtered_duplicates[' HISTOLOGY_SUBTYPE_1'].str.contains('(?:squamous_cell_carcinoma)$')]
print(oscc_filtered_data)
oscc_filtered_data.to_csv('../squamous-cell-carcinoma/squamous-cell-carcinoma-csv/02_cosmic_squamous_cell_only.csv',
                          index=False)


# # # 4. filter for non-specified samples
ns_filtered_data = filtered_duplicates[filtered_duplicates[' HISTOLOGY_SUBTYPE_1'].str.contains('(?:NS)$')]
print(ns_filtered_data)


# # # 5. count sample numbers
print(oac_filtered_data[' ID_SAMPLE'].value_counts())  # 323
print(oscc_filtered_data[' ID_SAMPLE'].value_counts())  # 788
print(ns_filtered_data[' ID_SAMPLE'].value_counts())  # 192


# # # 6. pie chart of histologies
pie_cancer_type = ['Adenocarcinoma', 'Squamous cell carcinoma', 'Not specified']
pie_cancer_counts = [323, 788, 192]

colour = color_discrete_sequence = px.colors.qualitative.G10

fig = go.Figure(data=[go.Pie(labels=pie_cancer_type, values=pie_cancer_counts)])
fig.update_traces(textposition='inside', marker=dict(colors=colour))
fig.update_layout(uniformtext_minsize=40, uniformtext_mode='hide')
fig.update_layout(legend=dict(orientation='v', yanchor='middle', y=1.02, xanchor='right', x=1,), legend_font_size=35)
fig.update_layout(title_x=0.02, title_y=0.05, title_xanchor='left', title_yanchor='bottom', title_font_size=20)
fig.show()
fig.write_image('../figures/pie-charts/oc_histologies.png', width=1920, height=1080, scale=1)
