import pandas

"""
Compare approved drugs, clinical trial drugs and adenocarcinoma SSL drug targets to the patient data to find which
patients are treatable
"""

patients = pandas.read_csv('../adenocarcinoma/adenocarcinoma-csv/03_filtered_cosmic_adenocarcinoma.csv',
                           header=0, sep=',')

approved_drugs = pandas.read_csv('../00-database-csv/approved_drugs_full.csv', header=0, sep=',')
clinical_trial_drugs = pandas.read_csv('../00-database-csv/clinical_trial_drugs_full.csv', header=0, sep=',')
repositionable_drugs = pandas.read_csv('../00-database-csv/adenocarcinoma_SSL_druggable_targets.csv', header=0, sep=',')



# # # 1. compare approved drugs against adenocarcinoma patients
gene_name = 'Gene Name'
drugs = 'Drug Name'
approved = 'Drug Status'
approved_vs_patients = approved_drugs.loc[approved_drugs['Gene Name'].isin(patients['GENE_NAME'])][[gene_name, drugs,
                                                                                                    approved]]
# print(approved_vs_patients)

# # # strip the p from mutations
patients[' MUTATION_AA'] = patients[' MUTATION_AA'].str.lstrip('p.')
# print(patients[' MUTATION_AA'])

# # # get patient ids and genes
patient = ' ID_SAMPLE'
gene = 'GENE_NAME'
mutation_aa = ' MUTATION_AA'
mutation = ' MUTATION_DESCRIPTION'
patients_vs_approved = patients.loc[patients['GENE_NAME'].isin(approved_drugs['Gene Name'])][[patient, gene,
                                                                                               mutation_aa, mutation]]
# print(patients_vs_approved)


# # #  merge the columns to get lists of patient id, genes and drugs
merge_columns = approved_vs_patients.merge(patients_vs_approved, left_on='Gene Name', right_on='GENE_NAME')
# print(merge_columns)

# # # drop the extra gene name column
drop_column = merge_columns.drop(['GENE_NAME'], axis=1)
drop_column.rename(columns={' ID_SAMPLE': 'Patient ID', ' MUTATION_AA': 'Mutation Aa',
                            ' MUTATION_DESCRIPTION': 'Mutation Description'}, inplace=True)
# print(drop_column)
# print(type(drop_duplicate_column))

# # # sort ID_SAMPLE by ID
approved_final = drop_column.set_index('Patient ID')
approved_final.sort_values(by=['Patient ID'], inplace=True)
approved_final = approved_final.drop_duplicates()
print(approved_final)
approved_final.to_csv('../personalised-regimes/personalised-regimes-csv/13_adenocarcinoma_approved_final.csv')


# # # 2. compare clinical trial drugs against adenocarcinoma patients
gene_name = 'Gene Name'
drugs = 'Drug Name'
clinical_trial = 'Drug Status'
clinial_trial_vs_patients = clinical_trial_drugs.loc[clinical_trial_drugs['Gene Name'].isin(patients['GENE_NAME'])][[
    gene_name, drugs, clinical_trial]]
# print(approved_vs_patients)

# # # strip the p from mutations
patients[' MUTATION_AA'] = patients[' MUTATION_AA'].str.lstrip('p.')
# print(patients[' MUTATION_AA'])

# # # get patient ids and genes
patient = ' ID_SAMPLE'
gene = 'GENE_NAME'
mutation_aa = ' MUTATION_AA'
mutation = ' MUTATION_DESCRIPTION'
patients_vs_clinical_trial = patients.loc[patients['GENE_NAME'].isin(clinical_trial_drugs['Gene Name'])][[patient, gene,
                                                                                                          mutation_aa,
                                                                                                          mutation]]
# print(patients_vs_approved)


# # #  merge the columns to get lists of patient id, genes and drugs
merge_columns = clinial_trial_vs_patients.merge(patients_vs_clinical_trial, left_on='Gene Name', right_on='GENE_NAME')
# print(merge_columns)

# # # drop the extra gene name column
drop_column = merge_columns.drop(['GENE_NAME'], axis=1)
drop_column.rename(columns={' ID_SAMPLE': 'Patient ID', ' MUTATION_AA': 'Mutation Aa',
                            ' MUTATION_DESCRIPTION': 'Mutation Description'}, inplace=True)
# print(drop_column)
# print(type(drop_duplicate_column))

# # # sort ID_SAMPLE by ID
ct_final = drop_column.set_index('Patient ID')
ct_final.sort_values(by=['Patient ID'], inplace=True)
ct_final = ct_final.drop_duplicates()
print(ct_final)
ct_final.to_csv('../personalised-regimes/personalised-regimes-csv/13_adenocarcinoma_clinical_trial_final.csv')


# # # 3. compare repositionable drugs against adenocarcinoma patients
repositionable_drugs.rename(columns={'Drugs': 'Drug Name', 'Gene Name B': 'SSL Partner'}, inplace=True)

gene_name = 'Gene Name A'
SSL_target = 'SSL Partner'
drugs = 'Drug Name'
respositionable = 'Drug Status'
repositionable_vs_patients = repositionable_drugs.loc[repositionable_drugs['Gene Name A'].isin(patients['GENE_NAME'])][[
    gene_name, drugs, respositionable, SSL_target]]
# print(approved_vs_patients)

# # # strip the p from mutations
patients[' MUTATION_AA'] = patients[' MUTATION_AA'].str.lstrip('p.')
# print(patients[' MUTATION_AA'])

# # # get patient ids and genes
patient = ' ID_SAMPLE'
gene = 'GENE_NAME'
mutation_aa = ' MUTATION_AA'
mutation = ' MUTATION_DESCRIPTION'
patients_vs_repositionable = patients.loc[patients['GENE_NAME'].isin(repositionable_drugs['Gene Name A'])][[patient,
                                                                                                            gene,
                                                                                                            mutation_aa,
                                                                                                            mutation]]
# print(patients_vs_approved)

# # #  merge the columns to get lists of patient id, genes and drugs
merge_columns = repositionable_vs_patients.merge(patients_vs_repositionable,
                                                 left_on='Gene Name A',
                                                 right_on='GENE_NAME')
# print(merge_columns)

# # # drop the extra gene name column
drop_column = merge_columns.drop(['GENE_NAME'], axis=1)
drop_column.rename(columns={' ID_SAMPLE': 'Patient ID', ' MUTATION_AA': 'Mutation Aa',
                            ' MUTATION_DESCRIPTION': 'Mutation Description'}, inplace=True)
# print(drop_column)
# print(type(drop_duplicate_column))

# # # sort ID_SAMPLE by ID
repositionable_final = drop_column.set_index('Patient ID')
repositionable_final.reset_index()
repositionable_final.sort_values(by=['Patient ID'], inplace=True)
repositionable_final = repositionable_final.drop_duplicates()
# print(repositionable_final)
repositionable_final.rename(columns={'Gene Name A': 'Gene Name'}, inplace = True)
repositionable_final.to_csv('../personalised-regimes/personalised-regimes-csv/13_adenocarcinoma_repositionable_final'
                            '.csv')


# # # 4. concatenate the data into one file

approved_final = pandas.read_csv('../personalised-regimes/personalised-regimes-csv/13_adenocarcinoma_approved_final.csv'
                                 , header=0, sep=',')
ct_final = pandas.read_csv('../personalised-regimes/personalised-regimes-csv/13_adenocarcinoma_clinical_trial_final.csv'
                           , header=0, sep=',')
repositionable_final = pandas.read_csv('../personalised-regimes/personalised-regimes-csv/13_adenocarcinoma_'
                                       'repositionable_final.csv', header=0, sep=',')


data = [approved_final, ct_final, repositionable_final]

complete = pandas.concat(data)
complete.sort_values(by=['Patient ID'], inplace=True)
# print(complete)
columns = complete.columns.tolist()
# print(columns)
complete = complete[['Patient ID', 'Gene Name', 'Mutation Aa', 'Mutation Description', 'Drug Name', 'Drug Status',
                     'SSL Partner']]
# print(complete)
complete.replace("Repositionable cancer drug - synthetic lethality", "Repositionable - SSL target", inplace=True)
complete.replace("Repositionable cancer drug - GoF oncogene", "Repositionable - GoF oncogene target", inplace=True)
print(complete)

complete.to_csv('FINAL_adenocarcinoma_complete_patient_regimes.csv', index=False)


# # # 5. count how many patients are treated with each treatment type
patient_count = complete.value_counts('Patient ID')
print(patient_count)

x = approved_final['Patient ID'].value_counts()
print(x)

x = ct_final['Patient ID'].value_counts()
print(x)

x = repositionable_final['Patient ID'].value_counts()
print(x)
