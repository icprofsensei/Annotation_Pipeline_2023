import pandas as pd
import os
import shutil
from pathlib import Path
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

filelist_test = []
for filename in os.scandir('allfiles/test'):
      if filename.is_file():
            j = os.path.abspath(filename)
            print(j)
            filelist_test.append(j)



os.makedirs('Microbiomes_test/airway')
os.mkdir('Microbiomes_test/faecal')
os.mkdir('Microbiomes_test/skin')
os.mkdir('Microbiomes_test/urinary')
os.mkdir('Microbiomes_test/vaginal')
for i in filelist_test:
        if str(i[71:]) in airway_ids:
            print('air')
            shutil.copyfile(i, 'Microbiomes_test/airway/' + str(i[71:]))
        if str(i[71:]) in faecal_ids:
            print('fae')
            shutil.copyfile(i, 'Microbiomes_test/faecal/' + str(i[71:]))
        if str(i[71:]) in skin_ids:
            print('sk')
            shutil.copyfile(i, 'Microbiomes_test/skin/' + str(i[71:]))
        if str(i[71:]) in urinary_ids:
            print('ur')
            shutil.copyfile(i, 'Microbiomes_test/urinary/' + str(i[71:]))
        if str(i[71:]) in vaginal_ids:
            print('vag')
            shutil.copyfile(i, 'Microbiomes_test/vaginal/' + str(i[71:]))


filelist_train = []
for filename in os.scandir('allfiles/train'):
      if filename.is_file():
            j = os.path.abspath(filename)
            print(j)
            filelist_train.append(j)



os.makedirs('Microbiomes_train/airway')
os.mkdir('Microbiomes_train/faecal')
os.mkdir('Microbiomes_train/skin')
os.mkdir('Microbiomes_train/urinary')
os.mkdir('Microbiomes_train/vaginal')
for i in filelist_train:
        if str(i[72:]) in airway_ids:
            print('air')
            shutil.copyfile(i, 'Microbiomes_train/airway/' + str(i[72:]))
        if str(i[72:]) in faecal_ids:
            print('fae')
            shutil.copyfile(i, 'Microbiomes_train/faecal/' + str(i[72:]))
        if str(i[72:]) in skin_ids:
            print('sk')
            shutil.copyfile(i, 'Microbiomes_train/skin/' + str(i[72:]))
        if str(i[72:]) in urinary_ids:
            print('ur')
            shutil.copyfile(i, 'Microbiomes_train/urinary/' + str(i[72:]))
        if str(i[72:]) in vaginal_ids:
            print('vag')
            shutil.copyfile(i, 'Microbiomes_train/vaginal/' + str(i[72:]))
