from explore_import import  *

def combine_ionbot(working_folder, dataset_name, db_sufix,file_name):
    #several folders, representing samples, containe ionbot results
    dir=os.path.join(working_folder, dataset_name, f'{dataset_name}{db_sufix}')
    toread=[os.path.join(lvl1[0],file_name) for lvl1 in os.walk(os.path.join(working_folder, dataset_name, f'{dataset_name}{db_sufix}')) if len(lvl1[0].split("/"))==6 and file_name in lvl1[2]]   
    combined=pd.DataFrame()
    #print(dir)
    #print(toread)
    for file in toread:
        #print(file)
        #print(file.lstrip(dir).rstrip(file_name))
        a=pd.read_csv(file)
        a['spectrum_file']=file.lstrip(dir).rstrip(file_name)
        combined=pd.concat([combined,a])
    return combined

def load_combined_first(working_folder,dataset_name,sufix="-combined-filtered-results.csv"):
    
    combined_first={"canon":0,"trembl":0,"openprot":0}
    #set path
    canonical_path = os.path.join(working_folder, dataset_name, f'{dataset_name}-canon', 
                              f'{dataset_name}-canon{sufix}')
    trembl_path    = os.path.join(working_folder, dataset_name, f'{dataset_name}-trembl', 
                                  f'{dataset_name}-trembl{sufix}')
    openprot_path  = os.path.join(working_folder, dataset_name, f'{dataset_name}-openprot', 
                                  f'{dataset_name}-openprot{sufix}')
    #load
    patt = re.compile(r'\(.+?\)') #for removing retention time in parentasies
    for path,db in zip([canonical_path,trembl_path,openprot_path],combined_first.keys()):
        df=pd.read_csv(path)
        df['spectrum_file'] = df.spectrum_file.str.replace('.RAW.mgf', '.mgf')
        df['modified_peptide'] = df.matched_peptide + '|' + df.modifications
        df['modifications_noRT'] = df.modifications.str.replace(patt, '', regex=True)
        combined_first[db]=df
    return combined_first

def classify_protein_io(row):
    if any('decoy' in i for i in row):
        return 'Decoy'
        
    row=[re.sub(r'\|\|.*','',x) for x in row]
    row=[re.sub(r'.*\(\(','',x) for x in row]
    row=[re.sub(r'\)\).*','',x) for x in row]
    if any('CONTAMINANT' in x or 'contaminant' in x for x in row):
        return 'Contam'
    elif all(x.startswith('II_') or x.startswith('IP_') for x in row):
        return 'NonCanon'
    else:
        return 'Canon'

def classify_protein_byprotname_io(x):
    if 'decoy' in x:
        return 'Decoy'

    if 'CONTAMINANT' in x or 'contaminant' in x:
        return 'Contam'
    elif x.startswith('II_') or x.startswith('IP_'):
        return 'NonCanon'
    else:
        return 'Canon'

def classify_protbypep_io(row):
    if row in ["Decoy","Contam"]:
        return row
    if "Noncanon" in row:
        return "NonCanon"
    elif "Canon" in row or "can_noncan" in row:
        return "Canon"
def extract_accession(x, mode):
    #the mode terms are confused here, but their application  is the following
    #open - is for canon and trembl database
    #closed - is for openprot
    #so to annotate all proteins correctly, use "closed" mode everywhere
    if mode=="open":
        match = re.match(r"^[^(]+", x)
        result = match.group(0) if match else None
    elif mode=="closed":
        match = re.search(r'\(\(([^()]+)\)\)$', x)
        result = match.group(1) if match else None
    return result
    
def classify_peptide_io(row,mode):
    #row=[re.sub(r'\(\(.*','',x) for x in row]
    row=[extract_accession(x,mode) for x in row]
    row=[x for x in row if not x is None]
    
    if any('decoy' in i for i in row):
        return 'Decoy'        
    #row=[re.sub(r'\|\|.*','',x) for x in row]
    #row=[re.sub(r'.*\(\(','',x) for x in row]
    #row=[re.sub(r'\)\).*','',x) for x in row]
    
    
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

def get_nonc_prots(pep_df,mode):
    return list(set([extract_accession(p, mode) for p in pep_df[pep_df.peptide_class.isin(["unique_to_Noncanon"])].proteins])) #'shared_in_Noncanon'

def get_frompepdf_pepcounts_io(pep_df,peptide_class,proteins):
    d=dict(zip(set(proteins),[0]*len(set(proteins))))
    cl_pep_df=pep_df[(pep_df.peptide_class==peptide_class)&(pep_df.proteins_expand.isin(set(proteins)))]
    cl_pep_df_counts=Counter(cl_pep_df.proteins_expand.tolist())
    d.update(cl_pep_df_counts)
    return  d
    
def get_frompsmdf_psmcounts_io(pep_df,psm_df,peptide_class,proteins):
    d=dict(zip(set(proteins),[0]*len(set(proteins))))
    cl_peps=pep_df[(pep_df.peptide_class==peptide_class)&(pep_df.proteins_expand.isin(set(proteins)))].database_peptide.tolist()
    cl_psm_df_counts=Counter(psm_df[psm_df.database_peptide.isin(cl_peps)].proteins_expand.tolist())
    del_keys=set(cl_psm_df_counts.keys())-set(proteins)
    for k in del_keys:
        del cl_psm_df_counts[k]
    d.update(cl_psm_df_counts)
    #d={k:v for k,v in d.items() if k in proteins} #exidantly sliped in proteins with same peps
    
    return d