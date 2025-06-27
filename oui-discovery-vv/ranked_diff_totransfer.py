#code transfer
import pickle
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#load data
with open('combo_optr_first_datasets.pickle', 'rb') as handle:
    combo_optr_first_datasets = pickle.load(handle)
with open('mgfspectraN_datasets.pickle', 'rb') as handle:
    mgfspectraN_datasets = pickle.load(handle)
#
def get_score_ranking(df,score_col,bin_size,bins_min,bins_max,prf_canon,prf_nonc):
    break_points=np.arange(bins_min, bins_max+bin_size, bin_size)[1:]#0,1
    a=(df[f"{score_col}_{prf_canon}"]-df[f"{score_col}_{prf_nonc}"])/df[f"{score_col}_{prf_canon}"]
    a=a.sort_values()
    #filter outliers
    data=np.array(a.values.tolist())
    filtered = data[(data > bins_min) & (data < bins_max)]
    if len(filtered)==0:
        return [np.nan]*len(break_points)
    bin_label_list=[]; missed=[]
    for idx,v in enumerate(filtered):
        appended=False
        for i,b in enumerate(break_points):
            if i!=0 and v<=b and v>break_points[i-1]:
                bin_label_list.append(i); appended=True
            elif i==0 and v<=b:
                bin_label_list.append(i); appended=True
            elif v==1:
                bin_label_list.append(len(break_points)); append=True
        if not appended:
            missed.append(v)
    count={k:v for k,v in sorted(Counter(bin_label_list).items())}
    #bin hight, sort by bin index
    bin_hights=[count[k] if k in list(count.keys()) else 0 for k in list(range(len(break_points)))]
    return bin_hights
#standartase bins spread for all samples
qval_bins_mins=[]; qval_bins_maxs=[]
psmscore_bins_mins=[]; psmscore_bins_maxs=[]
for dataset_name in combo_optr_first_datasets:
    mgfspectraN=mgfspectraN_datasets[dataset_name]
    for mgffile in mgfspectraN:
        mgffile=mgffile.replace('.RAW.mgf', '.mgf')
        df=combo_optr_first_datasets[dataset_name]
        dfalt=df.loc[(df["spectrum_file"]==mgffile) & (df["database_trembl"]=="T") & (df["database_open"]=="T") & ((df["only altprots"]==True) | (df["only novelisoprots"]==True)),]
        score_col="q-value"
        a=(dfalt[f"{score_col}_trembl"]-dfalt[f"{score_col}_open"])/dfalt[f"{score_col}_trembl"]
        a=a.sort_values()
        data=np.array(a.values.tolist())
        qval_bins_mins.append(np.min(data));qval_bins_maxs.append(np.max(data))
        score_col="psm_score"
        a=(dfalt[f"{score_col}_trembl"]-dfalt[f"{score_col}_open"])/dfalt[f"{score_col}_trembl"]
        a=a.sort_values()
        data=np.array(a.values.tolist())
        psmscore_bins_mins.append(np.min(data));psmscore_bins_maxs.append(np.max(data))
#calculate ranking for all samples in datasets, take median between samples of same datasets
q_bins_min,q_bins_max=int(np.mean(qval_bins_mins)-2*np.std(qval_bins_mins)),1
p_bins_min,p_bins_max=int(np.mean(psmscore_bins_mins)-2*np.std(psmscore_bins_mins)),1
qvalrank_tropalt_first_datasets={"PXD002057.v0.11.4":0,"PXD003594.v0.11.4":0,"PXD014258.v0.11.4":0}
psmscorerank_tropalt_first_datasets={"PXD002057.v0.11.4":0,"PXD003594.v0.11.4":0,"PXD014258.v0.11.4":0}
for dataset_name in combo_optr_first_datasets:
    qvalrank=[]; psmscorerank=[]
    mgfspectraN=mgfspectraN_datasets[dataset_name]
    for mgffile in mgfspectraN:
        mgffile=mgffile.replace('.RAW.mgf', '.mgf')
        df=combo_optr_first_datasets[dataset_name]
        dfalt=df.loc[(df["spectrum_file"]==mgffile) & (df["database_trembl"]=="T") & (df["database_open"]=="T") & ((df["only altprots"]==True) | (df["only novelisoprots"]==True)),]
        qvalrank.append(get_score_ranking(dfalt,"q-value",0.5,q_bins_min,q_bins_max,"trembl","open"))
        psmscorerank.append(get_score_ranking(dfalt,"psm_score",10,p_bins_min,p_bins_max,"trembl","open"))
    qvalrank_tropalt_first_datasets[dataset_name]=np.median(np.array(qvalrank),axis=0)
    psmscorerank_tropalt_first_datasets[dataset_name]=np.median(np.array(psmscorerank),axis=0)
#some elements for plot
q_bin_size=0.5
q_binsvals=np.arange(q_bins_min,q_bins_max+q_bin_size, q_bin_size)[1:]
qvalrank_med_error={i:0 for i in range(len(q_binsvals))}
qvalrank_med_avarage={i:0 for i in range(len(q_binsvals))}
qvalrank_med_swarmplot={i:0 for i in range(len(q_binsvals))}

for bins in range(len(q_binsvals)):
    qvalrank_medians=[l[bins] for l in qvalrank_tropalt_first_datasets.values()]
    qvalrank_med_error[bins]=stats.sem(qvalrank_medians, axis=None, ddof=0)
    qvalrank_med_avarage[bins]=np.mean(qvalrank_medians)
    qvalrank_med_swarmplot[bins]=qvalrank_medians
#same, psm score
p_bin_size=10
p_binsvals=np.arange(p_bins_min,p_bins_max+p_bin_size, p_bin_size)[1:]
psmrank_med_error={i:0 for i in range(len(p_binsvals))}
psmrank_med_avarage={i:0 for i in range(len(p_binsvals))}
psmrank_med_swarmplot={i:0 for i in range(len(p_binsvals))}
for bins in range(len(p_binsvals)):
    psmrank_medians=[l[bins] for l in psmscorerank_tropalt_first_datasets.values()]
    psmrank_med_error[bins]=stats.sem(psmrank_medians, axis=None, ddof=0)
    psmrank_med_avarage[bins]=np.mean(psmrank_medians)
    psmrank_med_swarmplot[bins]=psmrank_medians
#plot
dbs=np.array([round(i,3) for i in q_binsvals][::-1])*(-1)
dbs=[str(db) if i<=19 else f"{dbs[20]}-{dbs[len(dbs)-1]}" for i,db in enumerate(dbs[0:21])]
x=np.array(range(len(dbs)))
y=np.array(list(qvalrank_med_avarage.values()))[::-1]
y=[j if i<=19 else np.sum(y[20:]) for i,j in enumerate(y[0:21])]
yerr=np.array(list(qvalrank_med_error.values()))[::-1]
yerr=[j if i<=19 else np.sum(yerr[20:]) for i,j in enumerate(yerr[0:21])]
barcolor="mediumaquamarine"
fig, ax = plt.subplots(2,1,layout='constrained', figsize=(8,5))
fig.suptitle('Ranked difference between OpenProt and UniProt scores FDR 1%')
width=0.75
ax[0].bar(x, y, width=width, color=barcolor)
ax[0].errorbar(x, y, yerr=yerr, fmt="^",color="black",markersize=4)
ax[0].set_ylabel('Median count')
ax[0].set_xticks(x, dbs,fontsize=8,rotation=45)
ax[0].set_xlabel("binned q-value difference")
ax[0].set_ylim(0,150)

ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)

dbs=np.array([round(i,3) for i in p_binsvals][::-1])*(-1)
dbs=[str(db) if i<=19 else f"{dbs[20]}-{dbs[len(dbs)-1]}" for i,db in enumerate(dbs[0:21])]
x=np.array(range(len(dbs)))
y=np.array(list(psmrank_med_avarage.values()))[::-1]
y=[j if i<=19 else np.sum(y[20:]) for i,j in enumerate(y[0:21])]
yerr=np.array(list(psmrank_med_error.values()))[::-1]
yerr=[j if i<=19 else np.sum(yerr[20:]) for i,j in enumerate(yerr[0:21])]
barcolor="slateblue"

ax[1].bar(x, y, width=width, color=barcolor)
ax[1].errorbar(x, y, yerr=yerr, fmt="^",color="black",markersize=4)
a=[i for i in x for j in range(3)]

ax[1].set_ylabel('Median count')
ax[1].set_ylim(0,150)
ax[1].set_xticks(x, dbs,fontsize=8,rotation=45)
ax[1].set_xlabel("binned psm score difference")
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)

plt.savefig("F4BC_ranks_errorbar.svg")