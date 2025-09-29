import pandas as pd
import matplotlib.pyplot as plt
import math

from utils import get_q_vertical, plot_qval_thresholds, get_min_fdp, get_controled_q_vertical #get_q_fdp_firstids

MERGED_CSV=snakemake.input["merged_csv"] #contains reference file locations according to path_columns
INPUT_DIR=snakemake.input["input_dir"]
THRESHOLDS_CSV=snakemake.output["thresholds_csv"]
THRESHOLDS_SVG=snakemake.output["thresholds_svg"]
try:
    FILTER=snakemake.params["filter"] #dictionary column:value
except:
    FILTER=None
custom_order_search = snakemake.params["custom_order_search"]

#Convert to numeric
stats_df=pd.read_csv(MERGED_CSV)
stats_df["lower_bound"]=stats_df["Lower bound"].apply(lambda x: float(x.replace("%","")))
stats_df["upper_bound"]=stats_df["Paired method"].apply(lambda x: float(x.replace("%","")))

#Filter merged csv if needed; By default all conditions are logicaly AND (&)
if FILTER:
    for column,value in FILTER.items():
        if isinstance(value,list):
            stats_df = stats_df[stats_df[column].isin(value)]
        else:
            stats_df = stats_df[stats_df[column]==value]

#Create output df
path_columns=['sample', 'subset', 'dataset', 'database', 'search_type', 'group_type', 'approach']
threshold_df=pd.DataFrame(columns=path_columns)

for i, row in stats_df.iterrows():
    """
     If minimal FDP is < 1%, then select either q-value for FDP < 1% or a minimal FDP, which ever is available.
    """
    j=len(threshold_df)
    threshold_df.loc[j,path_columns]= row[path_columns]
    
    path = f"{INPUT_DIR}/{row.search_type}_{row.group_type}/{row.dataset}-{row.database}/{row.approach}/subsets/{row['sample']}.csv"
    print(path)
    df = pd.read_csv(path)

    df[["q_value", "paired_fdp", "lower_bound_fdp"]] *= 100

    #apply smoothing
    window_factor=0.005 #if len(df)>200 else 0.05
    ROLLING_WINDOW=int(math.ceil(window_factor*len(df)))
    ROLLING_WINDOW= 10 if ROLLING_WINDOW<1 else ROLLING_WINDOW
    df=df[["q_value", "paired_fdp", "lower_bound_fdp"]].rolling(window=ROLLING_WINDOW).mean()

    fdp_min, fdp_min_faild_5, fdp_min_q_max, fdp_1_q_max, q_sel , q_sel_fdp, q_sel_fail_5, fdp_sel_q_max_smooth , fdp_sel_q_min_controled,  q_sel_lower_bound_fdp = [None]*10

    #fdp_min = get_q_fdp_firstids(df)
    fdp_min = get_min_fdp(df)
    # Minimal FDP > 5% -> failed
    fdp_min_faild_5= fdp_min>=5
    # What is max q-value of smallest FDP
    indx_fdp_min=df.index[df["paired_fdp"] == fdp_min]
    fdp_min_q_max=df.loc[indx_fdp_min[-1], 'q_value']
    # max q-value for selected FDP - first guess - q-value coresponding to minomal FDP
    q_sel = fdp_min_q_max
    #FDP of selected q-value max - first guess - same as minimal FDP in sample
    q_sel_fdp = fdp_min
    #Next: if minimal FDP is under 1%, than give take 1%
    # What is max q-value for FDP 1%
    idx_fdp_1 = df.index[df["paired_fdp"] <= 1] #this can also select fdp=0 but not fdp>0, so..
    if len(idx_fdp_1) > 0:
        fdp_1_q_max = df.loc[idx_fdp_1[-1], 'q_value']
        # If max q-value for FDP 1% > max q-value of smallest FDP, select first
        q_sel= fdp_1_q_max #if fdp_1_q_max > fdp_min_q_max else fdp_min_q_max
        #Record q-value threshold, for which all identifications have FDP <= 1 %
        #q_sel_fdp = 1
        q_sel_fdp = df.loc[idx_fdp_1, 'paired_fdp'].max() #...consider it here
        fdp_sel_q_max_smooth = q_sel if q_sel <= 1 else get_q_vertical(df.copy(deep=True), q_sel_fdp)
        #Select lowest q-value (>0) for q_sel_fdp
        fdp_sel_q_min_controled = q_sel if q_sel <= 1 else get_controled_q_vertical(df.copy(deep=True),q_sel_fdp, 0, q_sel)
    else:
        fdp_1_q_max =  None
        #fdp_1_q_max_smooth = None
        fdp_sel_q_max_smooth = q_sel if q_sel <= 1 else get_q_vertical(df.copy(deep=True), q_sel_fdp)
        #Select lowest q-value (>0) for q_sel_fdp
        fdp_sel_q_min_controled = q_sel if q_sel <= 1 else get_controled_q_vertical(df.copy(deep=True),q_sel_fdp, 0, q_sel)
            
    #Does selected q-value excid 5%?
    q_sel_fail_5 = q_sel >= 5
    # FDP Lower bound for selected q-value threshold
    q_sel_lower_bound_fdp = df.loc[df.index[df["q_value"] == fdp_sel_q_min_controled][-1], 'lower_bound_fdp'] #q_sel
            
    threshold_df.loc[j, [
        "fdp_min", "fdp_min_faild_5",
        "fdp_min_q_max", "fdp_1_q_max",
        "q_sel" , "q_sel_fdp", "q_sel_fail_5", "fdp_sel_q_max_smooth", "fdp_sel_q_min_controled", "q_sel_lower_bound_fdp",
    ]] = [
        fdp_min, fdp_min_faild_5,
        fdp_min_q_max, fdp_1_q_max,
        q_sel , q_sel_fdp, q_sel_fail_5, fdp_sel_q_max_smooth, fdp_sel_q_min_controled, q_sel_lower_bound_fdp,
    ]

threshold_df.to_csv(THRESHOLDS_CSV,index=False)

#custom order on search_type for visualisation
custom_order_search_type = pd.CategoricalDtype(categories=custom_order_search, ordered=True)
threshold_df['search_type'] = threshold_df['search_type'].astype(custom_order_search_type)

#Plot
for (database, approach, subset),data in threshold_df.groupby(["database", "approach", "subset"]):
    suptitle=f"Relationship between selected q-value threshold and corresponding FDP\n in samples for {database} database, {approach} approach ({ 'all' if subset is None else subset} subset)"
    plot_qval_thresholds(data,suptitle)
    plt.savefig(THRESHOLDS_SVG)
    plt.close()
    