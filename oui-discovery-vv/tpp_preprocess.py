from explore_import import  *

def combine_tpp(work_dir,work_subdir,database,file_sufix,verbose=False):
    ''' Combine TPP output from one directory into one data frame.
    file_sufix: _PeptideProphet.pep.xml | _ProteinProphet.pep.prot.xml | .pep.xml
    '''
    path = os.path.join(work_dir, work_subdir, f'{work_subdir[0:9]}-{database}')
    if verbose: print(path)
    toread=[os.path.join(lvl1[0],lvl2) for lvl1 in os.walk(path) for lvl2 in lvl1[2] if lvl2.endswith(file_sufix)]
    if verbose: print(toread)
    combined=pd.DataFrame()
    for file in toread:
        if "ProteinProphet" in file_sufix:
            data=pyteomics.protxml.DataFrame(file)
        elif "PeptideProphet" in file_sufix:
            data=pyteomics.pepxml.DataFrame(file)
        spectrum_file=re.sub(f"\{file_sufix}$", "", file.lstrip(path))+".mgf"
        data['spectrum_file']=spectrum_file
        combined=pd.concat([combined,data])
    return combined

def get_database_tpp(row):
    if any("decoy_" in p for p in row):
        return "D"
    else:
        return "T"

def get_qval_tpp(df,score="fval"):
    """ calculate FDR and q-values based on score"""
    df = df.sort_values(by=score, ascending=False).reset_index(drop=True)        
    # Calculate cumulative counts of targets and decoys
    df['cum_targets'] = (df['database'] == 'T').cumsum()
    df['cum_decoys'] = (df['database'] == 'D').cumsum()
    # Calculate FDR
    df['FDR'] = df['cum_decoys'] / df['cum_targets']
    # cumulative minimum from bottom to top
    df['q-value'] = df['FDR'][::-1].cummin()[::-1]  
    return df

def get_idrate_tpp(mgfFiles, df):    
    IDrate={}
    for file in mgfFiles:
        if "RAW" in file: sup=re.sub(f"\{'.RAW.mgf'}$", "", file)
        sup=re.sub(f"\{'.mgf'}$", "", file)
        id=np.sum([1 for sp in df.spectrum if sup in sp])
        IDrate[file]=id/mgfFiles[file]
    return IDrate

#def classify_peptide_tpp(row):
#    if any('decoy_' in i for i in row):
#        return 'decoy'        
#    row=[p.split("|")[1] for p in row]
#    if any('CONTAMINANT' in x or 'contaminant' in x for x in row):
#        return 'Contam'
#    elif all(x.startswith('II_') or x.startswith('IP_') for x in row):
#        return 'NonCanon_unique'
#    elif all(not x.startswith('II_') and not x.startswith('IP_') and not 'CONTAMINANT' in x and not 'contaminant' in x for x in row):
#        return 'Canon_unique'
#    else:
#        return 'Shared'

def classify_peptide_tpp(row):
    if any('decoy_' in i for i in row):
        return 'decoy'        
    ##row=[p.split("|")[1] for p in row]
    if any('CONTAMINANT' in x or 'contaminant' in x for x in row):
        return 'Contam'
    elif all(x.startswith('II_') or x.startswith('IP_') for x in row):
        if len(row)==1: 
            return "unique_to_Noncanon"
        else:
            return 'shared_in_Noncanon'
    elif all(not x.startswith('II_') and not x.startswith('IP_') and not 'CONTAMINANT' in x and not 'contaminant' in x for x in row):
        if len(row)==1: 
            return "unique_to_Canon"
        else:
            return 'shared_in_Canon'
    else:
        return 'shared_btw_can_noncan'

def classify_protein_tpp(row):
    if any('decoy_' in i for i in row):
        return 'decoy'        
    ##row=[p.split("|")[1] for p in row]
    if any('CONTAMINANT' in x or 'contaminant' in x for x in row):
        return 'Contam'
    elif all(x.startswith('II_') or x.startswith('IP_') for x in row):
        return 'NonCanon'
        # Ensembl is canonical
    else:
        return 'Canon'

def get_protein_groups_tpp(df):
    tmp=pd.DataFrame(columns=["protein_name","protein_class"])
    for group_number,group in df.groupby("group_number"):
        #if anu is canonical
        ##ifwcanon=any((not p.split("|")[1].startswith("IP_")) and (not p.split("|")[1].startswith("II_")) for p in group.protein_name.tolist())
        ifwcanon=any((not p.startswith("IP_")) and (not p.startswith("II_")) for p in group.protein_name.tolist())
        prot_class="Canon" if ifwcanon else "NonCanon"
        tmp=pd.concat([tmp,pd.DataFrame({"protein_name":["||".join(group.protein_name.tolist())],"protein_class":[prot_class]})])
    return tmp
    
def get_openprot_proteins_tpp(df_pep,df_prot,nc_peptide_classes):
    inffered_proteins={"Canon":0,"NonCanon":0}
    prot_grs=get_protein_groups_tpp(df_prot)
    inffered_proteins["Canon"]=set(prot_grs[prot_grs.protein_class=="Canon"].protein_name.tolist()) #set(a_f)
    a=df_pep[df_pep.peptide_class.isin(nc_peptide_classes)].protein.tolist()
    inffered_proteins["NonCanon"]=set(list(chain(*a)))
    return inffered_proteins

def get_protein_level_tpp(row):
    return "Singleton" if len(row.split("||"))==1 else "Group"

def get_proteins_byrun_tpp(dict,database):
    #singeltons - lower limit, indistinguishable proteins - upper limit
    a=[] ; b=[]
    for dataset in dict:
        for spectrum_file in dict[dataset][database]:
             for cl in dict[dataset][database][spectrum_file]:
                 prots=dict[dataset][database][spectrum_file][cl]
                 sing=[p for p in prots if "||" not in p]
                 a=a+sing
                 uplim=[i for p in prots if "||" in p for i in p.split("||")]
                 #take away those in singletone to avoid dulicates                 
                 b=b+list(set(uplim)-set(sing))
    return set(a),set(b)

def get_frompepdf_pepcounts_tpp(pep_df,peptide_class,proteins):
    d=dict(zip(set(proteins),[0]*len(set(proteins))))
    cl_pep_df=pep_df[(pep_df.peptide_class==peptide_class)&(pep_df.proteins_expand.isin(set(proteins)))]
    cl_pep_df.drop_duplicates(["peptide","proteins_expand"],keep="first",inplace=True)
    cl_pep_df_counts=Counter(cl_pep_df.proteins_expand.tolist())
    d.update(cl_pep_df_counts)
    return  d
    
def get_frompsmdf_psmcounts_tpp(pep_df,peptide_class,proteins):
    d=dict(zip(set(proteins),[0]*len(set(proteins))))
    cl_pep_df=pep_df[(pep_df.peptide_class==peptide_class)&(pep_df.proteins_expand.isin(set(proteins)))]
    cl_pep_df_counts=Counter(cl_pep_df.proteins_expand.tolist())
    d.update(cl_pep_df_counts)
    return  d

def get_reproduce_proteins_tpp(dict,database):
    #singeltons, upper limit
    a=[] ; b=[]
    for dataset in dict:
        for spectrum_file in dict[dataset][database]:
             for cl in dict[dataset][database][spectrum_file]:
                 prots=dict[dataset][database][spectrum_file][cl]
                 sing=[p for p in prots if "||" not in p]
                 a=a+list(set(sing)) #take 1 count on protein in file
                 uplim=[i for p in prots if "||" in p for i in p.split("||")]
                 #take away those in singletone to avoid dulicates                 
                 b=b+list(set(uplim)-set(sing))
    return Counter(a),Counter(b)

def get_class_op_proteins(dict,database,cl):
    #singeltons, upper limit
    a=[] ; b=[]
    for dataset in dict:
        for spectrum_file in dict[dataset][database]:
            prots=dict[dataset][database][spectrum_file][cl]
            sing=[p for p in prots if "||" not in p]
            a=a+sing
            uplim=[i for p in prots if "||" in p for i in p.split("||")]
            #take away those in singletone to avoid dulicates                 
            b=b+list(set(uplim)-set(sing))
    return set(a),set(b)