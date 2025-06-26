#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import pyteomics.auxiliary as aux
import re, os

from datetime import date
DATE = date.today().strftime('%Y%m%d')

project_palette = {
    "canon":"orangered",
    "trembl":"yellowgreen",
    "openprot":"cornflowerblue"
}

def cohenD(X,Y):
    nx,mx,sdx = len(X),np.mean(X),np.std(X)
    ny,my,sdy = len(Y),np.mean(Y),np.std(Y)
    sdp = np.sqrt(((nx-1)*sdx**2 + (ny-1)*sdy**2)/(nx+ny-2))
    return round(abs(mx-my)/sdp, 3)

# Importing and filtering results
def import_pep_IDs(PATH, filtering=False, drop_contaminants=True):
    df = pd.read_csv(PATH, usecols=['spectrum_title','scan','spectrum_file','matched_peptide','database_peptide',
                                     'modifications','leadprot','database','precursor_mass',
                                     'isCanonical','isModified', 'psm_score',
                                     'q.value',
                                     'group_qval',
                                    'custom_q']
                    )
    if drop_contaminants:
        df = df[df.isCanonical!='Contam'].copy(deep=True)

    if filtering=='global':
        df = df[df['q.value']<0.01].copy(deep=True)
    elif filtering=='groupwalk':
        df = df[df.group_qval<0.01].copy(deep=True)
    elif filtering=='custom':
        df = df[df.custom_q<0.01].copy(deep=True)
    elif filtering=='hybrid':
        tmp = []
        for pippo,pluto in df.groupby('isCanonical').__iter__():
            if pippo=='Canonical':
                tmp.append(pluto[pluto['q.value']<0.01])
            elif pippo=='NonCanonical':
                tmp.append(pluto[pluto.custom_q<0.01])
        df = pd.concat(tmp, ignore_index=True)
        del tmp
    elif filtering: 
        # gives error if filtering is not False
        print(f'Error! Filtering = {filtering}')
        return filtering

    # fixes issue with some files being .RAW and other being .raw
    df.spectrum_file = df.spectrum_file.apply(lambda x: x.split('.')[0])
    # in some mgf files the 'spectrum title' includes the file name, making the spectrum title unique.
    # when the file name is NOT included, spectrum titles are NOT unique, and this can mess up some analysis.
    df.spectrum_title = df.spectrum_file + ':' + df.spectrum_title.apply(lambda x: x.split(':')[-1])
    
    df['modified_peptide'] = df.matched_peptide + '|' + df.modifications
    # to remove retention times in parentheses
    df['modifications_noRT'] = df.modifications.str.replace(re.compile(r'\(.+?\)'), '', regex=True)
    
    return df


# FDR recalculation
def process_proteins(protein_str):
    """
    Process a protein string:
    1. Split by '||'.
    2. If more than 100 proteins, consider only the first 100.
    3. For each protein, split using the regex pattern and extract the fourth element.
    """
    proteins = protein_str.split('||')
    proteins = proteins[:100] 
    return [re.split(r'\(\(|\)\)', p)[3] for p in proteins]

def classify_leadprot(x):
    if 'CONTAMINANT' in x.upper():
        return 'Contam'
    elif x.startswith('II_') or x.startswith('IP_'):
        return 'NonCanon'
        # Ensembl is canonical
    else:
        return 'Canon'

def is_peptide_canonical(x):
    '''x is the list of protein classes'''
    if np.array([_=='Contam' for _ in x]).any():
        return 'Contam'
    if np.array([_=='Canon' for _ in x]).any():
        return 'Canonical'
    return 'NonCanonical'

def classifiy_mods(row):
    if row.modifications=='Unmodified':
        return 'Unmodified'
    elif len(row.unexpected_modification)>1:
        return 'Unexpected'
    else:
        return 'Expected'

def custom_subgroup_filter(data_):
    filtered_subgroups = []
    for (c,m),df in data_.groupby(['isCanonical','isModified']).__iter__():
        tmp = aux.target_decoy.qvalues(df, key='psm_score', reverse=True, is_decoy=df.database=='D',
                                      formula=1, full_output=True, q_label='custom_q')
        filtered_subgroups.append(tmp)

    return pd.concat(filtered_subgroups, ignore_index=True)

# Protein inference and FDR
def get_max_psm_score(x, psms_):
    protein_group = x.split(',')
    tmp3 = psms_[psms_.proteins.isin(protein_group)]
    return np.max(tmp3.psm_score)

def group_is_decoy(x):
    protein_group = x.split(',')
    tmp = [_.startswith('decoy') for _ in protein_group]
    return set(tmp)

def classify_protein_group(x):
    protein_group = x.split(',')
    tmp = [classify_leadprot(_) for _ in protein_group]
    return is_peptide_canonical(tmp)



### SANKEY PLOT FUNCTIONS ###
import seaborn as sns
import plotly.graph_objects as go

def define_colors(x):
    mypalette = sns.color_palette("Paired").as_hex()
    mypalette
    if x.startswith('Canonical+Unmodified/Expected'):
        return mypalette[0]
    elif x.startswith('NonCanonical+Unmodified/Expected'):
        return mypalette[1]
    elif x.startswith('Canonical+Unexpected'):
        return mypalette[2]
    elif x.startswith('NonCanonical+Unexpected'):
        return mypalette[3]
    elif x.startswith('Unidentified'):
        return mypalette[-2]
    elif x.startswith('Decoy'):
        return mypalette[-4]
    
def get_sankey_label(row, suffixes_):
    tmp = []
    for x in suffixes_:
        if row['database'+x]=='D':
            tmp.append('Decoy')
        elif row['database'+x]=='T':
            tmp2 = 'Unexpected' if row['isModified'+x]=='Unexpected' else 'Unmodified/Expected'
            tmp.append(row['isCanonical'+x] + '+' + tmp2)
        else:
            tmp.append('Unidentified')
    return tmp

def make_sankey_plot_with_counts(data, suffixes=['_trembl','_open']):
    data0 = data.copy(deep=True)
    sankey_labels = ['sankey'+x for x in suffixes]
    data0[sankey_labels] = data0.apply(lambda row: get_sankey_label(row, suffixes), result_type='expand', axis=1)
    
    sankey_l_counts = pd.DataFrame(data0[sankey_labels[0]].value_counts()).reset_index()
    sankey_l_counts[sankey_labels[0]+'_2'] = sankey_l_counts[sankey_labels[0]] + ' (' + sankey_l_counts['count'].apply(str) + ')'
    sankey_l_counts = sankey_l_counts.set_index(sankey_labels[0]).to_dict()[sankey_labels[0]+'_2']
    
    sankey_r_counts = pd.DataFrame(data0[sankey_labels[1]].value_counts()).reset_index()
    sankey_r_counts[sankey_labels[1]+'_2'] = sankey_r_counts[sankey_labels[1]] + ' (' + sankey_r_counts['count'].apply(str) + ')'
    sankey_r_counts = sankey_r_counts.set_index(sankey_labels[1]).to_dict()[sankey_labels[1]+'_2']
    
    leftLabels  = list(sankey_l_counts.values())
    rightLabels = list(sankey_r_counts.values())
    labels = leftLabels + rightLabels
    
    labels_colors = {_:define_colors(_) for _ in labels}
    
    data2 = pd.DataFrame(data0[sankey_labels].value_counts()).reset_index()
    data2[sankey_labels[0]] = data2[sankey_labels[0]].apply(lambda x: sankey_l_counts[x])
    data2[sankey_labels[1]] = data2[sankey_labels[1]].apply(lambda x: sankey_r_counts[x])
    data2['source'] = data2[sankey_labels[0]].apply(labels.index)
    data2['target'] = data2[sankey_labels[1]].apply(labels.index)
    
    fig = go.Figure(data=[go.Sankey(
        node = dict(
          pad = 15,
          thickness = 20,
          line = dict(color = "black", width = 0.5),
          label = leftLabels+rightLabels,
          color = list(labels_colors.values())
        ),
        link = dict(
          source = data2['source'],
          target = data2['target'],
          value  = data2['count']
      ))])
    fig.update_layout(title_text=f'Closed search (left) vs Open PTM search (right)', font_size=10)
    fig.update_layout(height=800, width=800,)

    data3 = pd.DataFrame(data0[sankey_labels].value_counts()).reset_index()
    data3 = data3.pivot(index=sankey_labels[0], columns=sankey_labels[1] , values='count')
    data3 = data3.fillna(0).astype(int)
    return fig, data3

# Auto save & export notebook to html!!
from ipylab import JupyterFrontEnd
import subprocess
import time

def autosave(FLD='publication-data', extra_labels=''):
    app = JupyterFrontEnd()
    app.commands.execute('docmanager:save')
    time.sleep(5) 
    
    nbname = os.path.split(os.environ.get("JPY_SESSION_NAME"))[-1]
    
    command = ['jupyter', 'nbconvert', nbname, '--to', 'html', '--output', 
               os.path.join(FLD, f'{DATE}-{nbname.replace('.ipynb','')}{extra_labels}.html')]
    _ = subprocess.run(command)
    print(_.returncode==0)