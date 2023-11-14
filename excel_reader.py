import pandas as pd
df = pd.read_excel('microbiomeCorpusToDomain.xlsx')

df_airway = df[df['Airway'] == 1]
df_faecal = df[df['Faecal'] == 1]
df_skin = df[df['Skin'] == 1]
df_urinary = df[df['Urinary'] == 1]
df_vaginal = df[df['Vaginal'] == 1]

airway_ids = df_airway['PMCID'].tolist()
faecal_ids = df_faecal['PMCID'].tolist()
skin_ids = df_skin['PMCID'].tolist()
urinary_ids = df_urinary['PMCID'].tolist()
vaginal_ids = df_vaginal['PMCID'].tolist()

print(airway_ids)
'''
print(df_airway)
print(df_faecal)
print(df_skin)
print(df_urinary)
print(df_vaginal)
'''