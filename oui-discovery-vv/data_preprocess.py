from explore_import import  *

def list_files(startpath):
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * (level)
        print('{}{}/'.format(indent, os.path.basename(root)))
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            print('{}{}'.format(subindent, f))
            
def get_total_nspectra(mgf_path):
    mgfFiles={file.split("/")[-1]:np.nan for file in glob.glob(f"{mgf_path}*.mgf")}
    for mgf_file in mgfFiles.keys():
        n_spectrum=0
        with mgf.read(f"{mgf_path}{mgf_file}") as reader:
            for spectrum in reader:
                n_spectrum+=1
        mgfFiles[mgf_file]=n_spectrum
    return mgfFiles

def get_idrate(mgfFiles, database, df):    
    IDrate={}
    for file in mgfFiles:
        print(file,len(df[df["spectrum_file"]==file.split(".")+f"-{database}"]))
        #file=file.replace('.RAW.mgf', '.mgf')
        IDrate[file]=len(df[df["spectrum_file"]==file.split(".")+f"-{database}"])/mgfFiles[file]
    return IDrate

