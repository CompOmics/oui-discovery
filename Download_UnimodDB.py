#!/usr/bin/env python
# coding: utf-8

# > Unimod database schema:
# > https://www.unimod.org/pdf/unimod_schema.pdf

import urllib.error
import urllib.request
import numpy as np
import pandas as pd
import requests, xmltodict, time

def Download_UnimodDB():
    response = requests.get('https://www.unimod.org/xml/unimod_tables.xml')
    time.sleep(3)
    data = xmltodict.parse(response.content.decode())

    specificity = pd.DataFrame(data['unimod']['specificity']['specificity_row'])
    specificity['mod_to_specificity']   = specificity['@mod_key'].apply(int)
    specificity['class_to_specificity'] = specificity['@classifications_key'].apply(int)

    modifications = pd.DataFrame(data['unimod']['modifications']['modifications_row'])
    modifications['mod_to_specificity'] = modifications['@record_id'].apply(int)

    classifications = pd.DataFrame(data['unimod']['classifications']['classifications_row'])
    classifications['class_to_specificity'] = classifications['@record_id'].apply(int)


    final_cols = [
        "@record_id_x",
        "@code_name",
        "@full_name",
        "@avge_mass",
        "@mono_mass",
        "@composition",
        "@one_letter",
        "@classification",
        "misc_notes_x",
        "misc_notes_y"
    ]

    final_cols_2 = [
        "unimod_id",
        "code_name",
        "full_name",
        "avg_mass",
        "mono_mass",
        "composition",
        "residue",
        "classification",
        "misc_notes_x",
        "misc_notes_y"
    ]


    final_df = modifications.merge(specificity, on="mod_to_specificity", how="outer")
    final_df = final_df.merge(classifications, on='class_to_specificity', how="outer")
    final_df = final_df[final_cols]
    final_df.columns = final_cols_2
    final_df = final_df[final_df.unimod_id==final_df.unimod_id]
    
    final_df.mono_mass = final_df.mono_mass.apply(lambda x: round(float(x), 4))
    final_df.unimod_id = final_df.unimod_id.apply(int)
    
    return final_df

def Download_Unimod_Dict():
    unimod = Download_UnimodDB()
    condensed_unimod = {}
    for uni_id,df in unimod.groupby("unimod_id").__iter__():
        condensed_unimod[uni_id] = {}
        condensed_unimod[uni_id]['residues'] = "".join([_ for _ in df.residue if '-' not in _])
        condensed_unimod[uni_id]['mono_mass'] = [_ for _ in set(df.mono_mass)][0]
        
    return condensed_unimod    

def getPTMmass(x,unimod_dict_):
    try:
        return unimod_dict_[x]['mono_mass']
    except:
        return np.nan