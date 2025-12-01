#!/usr/bin/env python3
import argparse
import sys
import logging
from pathlib import Path

import pandas as pd
import numpy as np
from pyteomics import pepxml
import pyteomics.auxiliary as aux
import ast
import os, re, subprocess

#python msfragger_prepare_to_fdrbench.py -f ./msfragger_pepXML/PXD002057/PXD002057-canon/canon-130327_o2_01_hu_c1_2hr.pep.xml \
# -out ./msfragger_pepXML/PXD002057/PXD002057-canon/ \
# -a /mnt/d/vasv1101/Documents_discC/pide_reanalysis/entrapment/entrapment-concat-fasta/entrapment-concat-fasta/Human_2023_01_canonical_ionbot_db/annotated_entrapment_fasta.csv \
# -v

class EntrapmentPreprocessFragPipe:
    def __init__(self, pepxml_path: Path, anndf: pd.DataFrame):
        self.pepxml_path = pepxml_path
        self.df: pd.DataFrame | None = None
        self.anndf: pd.DataFrame = anndf
        self._load()
        self.scores=['fval','peptideprophet_probability']
        self.select_columns=['peptide','global_q','custom_q','glob_cust_hybrid','fval','FDRGroup','remove']
        self.unmodified_threshold=0.01 #define as the highest mass difference which was not assigned any modification in non-entrapment run (module from Q3 of distribution, round up)
        #plt.boxplot(unmod_massdiff, showfliers=False);
        #Q1=np.percentile(unmod_massdiff, 25)
        #Q3=np.percentile(unmod_massdiff, 75)
        #IQR=Q3-Q1
        #print(Q1 - (1.5 * IQR),Q3+ (1.5 * IQR))

    def _load(self):
        try:
            self.df = pepxml.DataFrame(str(self.pepxml_path))
        except Exception as e:
            logging.error("Failed to load pepXML '%s': %s", self.pepxml_path, e)
            raise
    
    @staticmethod
    def load_annotation(ann_path: Path) -> pd.DataFrame:
        try:
            anndf = pd.read_csv(str(ann_path))
            anndf["isCanonical"] = anndf.peptide_class.apply(
                lambda x: x.split('_')[1] if '_' in x else x
            )
            anndf["isCanonical"].replace(
                {"Canon":"Canonical","NonCanon":"NonCanonical","Contam":"Contam"},
                inplace=True
            )
            anndf['proteins'] = anndf['proteins'].apply(
                lambda x: ast.literal_eval(x) if not pd.isna(x) else np.nan
            )
            #remove!
            #anndf["entrapment_proteins"] = anndf.apply(
            #    lambda x: ";".join(x.proteins)+x.peptide_type
            #    if isinstance(x.proteins,list) else np.nan, axis=1
            #)
            return anndf
        except Exception as e:
            logging.error("Failed to load peptides annotation file '%s': %s", ann_path, e)
            raise
    
    @staticmethod
    def classify_leadprot(x) -> str:
        x=x.replace("decoy_","")
        if 'CONTAMINANT' in x.upper():
            return 'Contam'
        elif x.startswith('II_') or x.startswith('IP_'):
            return 'NonCanon'
            # Ensembl is canonical
        else:
            return 'Canon'
    @staticmethod
    def is_peptide_canonical(x) -> str:
        '''x is the list of protein classes'''
        if np.array([_=='Contam' for _ in x]).any():
            return 'Contam'
        if np.array([_=='Canon' for _ in x]).any():
            return 'Canonical'
        return 'NonCanonical'
    @staticmethod
    def is_target(x) -> bool:
        return not 'decoy_' in x
    @staticmethod
    def mask_in_fdrbench(x) -> bool:
        '''mask decoys and contaminants'''
        return 'CONTAMINANT' in x.upper() or 'decoy_' in x
    
    def classify_mods(self, row) -> str:
        if len(row.modifications) == 0:
            if abs(row.massdiff) <= self.unmodified_threshold:
                return 'Unmodified'
            else:
                return 'Unexpected'
        else:
            return 'Expected'


    def custom_subgroup_filter(self,df,key):
        filtered_subgroups = []
        for (c,m),subdf in df.groupby(['isCanonical','isModified']).__iter__():
            tmp = aux.target_decoy.qvalues(subdf, key=key, reverse=True, is_decoy=subdf.isTarget==False,
                                        formula=1, full_output=True, q_label=f'custom_q_{key}')
            filtered_subgroups.append(tmp)

        return pd.concat(filtered_subgroups, ignore_index=True)
    
    def _validate_columns(self, df):
        required_cols = ("peptide", "protein", "modifications", "fval", "peptideprophet_probability")
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            logging.error("Missing required columns %s in '%s'; aborting.", missing, self.pepxml_path)
            return False
        if df.empty:
            logging.warning("Input DataFrame is empty for '%s'. Nothing to do.", self.pepxml_path)
            return False
        return True

    def _coerce_fval_column(self, df):
        if not np.issubdtype(df["fval"].dtype, np.number):
            logging.warning("Column 'fval' is not numeric; attempting coercion.")
            df["fval"] = pd.to_numeric(df["fval"], errors="coerce")
        nan_count = df["fval"].isna().sum()
        if nan_count:
            logging.info("Found %d NaNs in 'fval' after coercion.", nan_count)

    def _assign_target_decoy(self, df):
        try:
            #df["isTarget"] = df.protein.apply(lambda x: all(self.is_target(_) for _ in x))
            df["isTarget"] = df.peptide_header_class.apply(lambda x: x!='Decoy' if isinstance(x,str) else True)
            n_targets = df["isTarget"].sum()
            n_decoys = (~df["isTarget"]).sum()
            logging.info("Target/decoy breakdown: %d targets, %d decoys.", int(n_targets), int(n_decoys))
            if n_targets == 0 or n_decoys == 0:
                logging.warning("Only one class present; FDR estimation may be unreliable.")
            return True
        except Exception as e:
            logging.exception("Error during target/decoy assignment: %s", e)
            return False

    def _compute_global_qvalues(self, df, key):
        try:
            df_out = aux.target_decoy.qvalues(
                df, key=key, reverse=True,
                is_decoy=(df.isTarget == False),
                q_label=f"global_q_{key}", formula=1, full_output=True
            )
            return df_out
            #self.df = df_out  # update if needed later
            #return "global_q" in df_out.columns
        except Exception as e:
            logging.exception("Failed to compute global q-values: %s", e)
            return False

    def _classify_protein_and_mods(self, df, anndf):
        try:
            #df["protein_classes"] = df.protein.apply(lambda x: np.unique([self.classify_leadprot(_) for _ in x]))
            #df["isCanonical"] = df.protein_classes.apply(self.is_peptide_canonical)
            logging.info("Merge annotated fasta with pepxml by peptide sequence to identify peptides class.")
            df=df.merge(anndf[["sequence_x","isCanonical","peptide_header_class"]],left_on="peptide",right_on="sequence_x",how="left") #,"entrapment_proteins"
            #optionally drop 'sequence_x' column after merge
            df.drop(columns=["sequence_x"], inplace=True)
            #some peptides are duplicated due to I->L substitution. If the ntrie is matched more than once to same peptide class - deduplicate
            before_dedup = len(df)
            df.drop_duplicates(["peptide","isCanonical","spectrumNativeID"], inplace=True)
            after_dedup = len(df)
            dropped = before_dedup - after_dedup
            logging.info("Dropped %d duplicated rows based on peptide, isCanonical, and spectrumNativeID.", dropped)

            n_canon = len(df[df["isCanonical"]=="Canonical"])
            n_noncanon = len(df[df["isCanonical"]=="NonCanonical"])
            n_cantom = len(df[df["isCanonical"]=="Contam"])
            logging.info("Classes breakdown: %d Canon, %d NonCanon, %d Contam.", int(n_canon), int(n_noncanon), int(n_cantom))

            df["isModified"] = df.apply(self.classify_mods, axis=1)
            invalid = df[~df["isCanonical"].isin(["Canonical", "NonCanonical", "Contam"])].shape[0]
            if invalid:
                logging.warning("%d rows have invalid 'isCanonical' values.", invalid)
            return df
        except Exception as e:
            logging.exception("Error classifying proteins/modifications: %s", e)
            return False

    def _filter_custom_subgroups(self, df, key):
        try:
            return self.custom_subgroup_filter(df, key)
        except Exception as e:
            logging.exception("custom_subgroup_filter failed: %s", e)
            return None

    def _add_hybrid_scores(self, df2, key):
        try:
            if f"custom_q_{key}" not in df2.columns or f"global_q_{key}" not in df2.columns:
                logging.error(f"Missing 'custom_q_{key}' or 'global_q_{key}' for hybrid scoring.")
                return False
            df2[f"glob_cust_hybrid_{key}"] = df2.apply(
                lambda x: x[f"custom_q_{key}"] if x.isCanonical == "NonCanonical" else x[f"global_q_{key}"], axis=1
            )
            return True
        except Exception as e:
            logging.exception("Error computing hybrid scores: %s", e)
            return False

    def _prepare_group_walk(self, df2):
        try:
            df2["FDRGroup"] = df2.isCanonical.astype(str) + "_" + df2.isModified.astype(str)
            return True
        except Exception as e:
            logging.exception("Error preparing Group-walk grouping: %s", e)
            return False

    def _mark_psms_to_remove(self, df2):
        try:
            #df2["remove"] = df2.protein.apply(lambda x: any(self.mask_in_fdrbench(_) for _ in x))
            df2["remove"] = df2.peptide_header_class.apply(lambda x: isinstance(x,str))
            return True
        except Exception as e:
            logging.exception("Error applying FDRBench mask: %s", e)
            return False

    def _log_summary(self, df2):
        total = len(df2)
        removed = df2["remove"].sum()
        logging.info("Prepared %d PSMs for FDRBench; %d marked for removal (%.1f%%).", total, int(removed), 100 * removed / total if total else 0)
        for key in self.scores:
            nulls = df2[[colbase+f"_{key}" for colbase in ["global_q", "custom_q", "glob_cust_hybrid"]] ].isna().sum()
            logging.debug("Null counts in q-value columns: %s", nulls.to_dict())

    def run_groupwalk(self, r_path, script_path, input_path):
        logging.debug("GroupWalk run start...")
        dataset_dir=input_path.parent
        working_dir=script_path.parent
        file_name=input_path.name
        _ = subprocess.run([r_path, script_path, dataset_dir, working_dir, 
                            file_name])


    def prepare_to_fdrbench(self):
        if self.df is None:
            raise RuntimeError("DataFrame not loaded")

        df = self.df.copy()
        anndf=self.anndf.copy()

        if not self._validate_columns(df):
            return False

        self._coerce_fval_column(df)

        df = self._classify_protein_and_mods(df, anndf)
        if df is False:
            return False

        if not self._assign_target_decoy(df):
            return False
        
        for key in self.scores:

            df = self._compute_global_qvalues(df,key)

            if df is None:
                return False

            #if not self._classify_protein_and_mods(df, anndf):
            #    return False

            df = self._filter_custom_subgroups(df,key)
            #logging.debug("Columns in df2 after custom_subgroup_filter: %s", df2.columns.tolist())

            if df is None:
                return False

            if not self._add_hybrid_scores(df,key):
                return False

        if not self._prepare_group_walk(df):
            return False

        if not self._mark_psms_to_remove(df):
            return False

        self._log_summary(df)

        return df


def find_pepxml_files(directory: Path) -> list[Path]:
    return sorted(directory.glob("*.pep.xml"))

def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare PeptideProphet pepXML output to FDRBench run. Calculate FDR and mask decoys/contaminants."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--pepxml", "-f",
        type=Path,
        help="Path to a single PeptideProphet pepXML file to validate."
    )
    group.add_argument(
        "--working-dir", "-d",
        type=Path,
        help="Directory containing one or more .pep.xml files to validate."
    )
    parser.add_argument(
        "--out-dir", "-out", required=True,
        type=Path,
        help="Directory where to save preprocessed files."
    )
    parser.add_argument(
        "--ann-path", "-a", required=True,
        type=Path,
        help="Database annotation file with peptide sequence and peptide_class information."
    )
    parser.add_argument(
        "--r-path", "-r", required=True,
        type=Path,
        help="Path to r executable."
    )
    parser.add_argument(
        "--script-path", "-s", required=True,
        type=Path,
        help="Path to Run_group_walk.R. Note, Group_walk.R file has to be in the same directory."
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging."
    )

    return parser.parse_args()

def setup_logging(verbose: bool, log_file: Path | None = None):
    handlers = []
    fmt = "%(asctime)s %(levelname)s: %(message)s"
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(logging.Formatter(fmt))
        handlers.append(fh)
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(logging.Formatter(fmt))
    handlers.append(sh)
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, handlers=handlers, force=True)

def main():
    args = parse_args()

    # determine where to write the log: next to the pepXML(s)
    if args.pepxml:
        if not args.pepxml.exists():
            logging.error("Specified pepXML file does not exist: %s", args.pepxml)
            sys.exit(1)
        log_dir = args.pepxml.parent
        pepxml_paths = [args.pepxml]
    else:
        if not args.working_dir.is_dir():
            logging.error("Working directory is not a directory or does not exist: %s", args.working_dir)
            sys.exit(1)
        log_dir = args.working_dir
        pepxml_paths = find_pepxml_files(args.working_dir)
        if not pepxml_paths:
            logging.error("No .pep.xml files found in directory: %s", args.working_dir)
            sys.exit(1)

    log_path = log_dir / "entrapment_preprocess.log"
    setup_logging(args.verbose, log_file=log_path)

    logging.debug("Parsed arguments: %s", args)

    if not args.out_dir.exists():
        logging.info("Output directory does not exist. Creating: %s", args.out_dir)
        args.out_dir.mkdir(parents=True, exist_ok=True)
    elif not args.out_dir.is_dir():
        logging.error("Specified output path exists but is not a directory: %s", args.out_dir)
        sys.exit(1)  

    # Load annotation once
    if not args.ann_path.exists():
        logging.error("The annotation file does not exist %s", args.ann_path)
        sys.exit(1)
    anndf = EntrapmentPreprocessFragPipe.load_annotation(args.ann_path)           

    any_failed = False
    for pepxml_file in pepxml_paths:
        try:
            preprocessor = EntrapmentPreprocessFragPipe(pepxml_file, anndf)
            df_prepped = preprocessor.prepare_to_fdrbench()
            if df_prepped is False:
                logging.error("Preprocessing returned no data for %s; skipping.", pepxml_file)
                any_failed = True
                continue
            output_path = args.out_dir / (str(pepxml_file.stem).split(".")[0] + "_prep.csv")
            df_prepped.to_csv(output_path, index=False)
            logging.info("Wrote preprocessed output to %s", output_path)
            try:
                preprocessor.run_groupwalk(args.r_path, args.script_path, output_path)
                logging.info("Wrote group-walk output to %s", "groupwalk_output_"+output_path.name)
            except Exception:
                logging.exception("Exception while runing group-walk for '%s'", output_path)
                any_failed = True
        except Exception:
            logging.exception("Exception while preprocessing '%s'", pepxml_file)
            any_failed = True

    if any_failed:
        logging.error("One or more preprocessing files failed.")
        sys.exit(2)



if __name__ == "__main__":
    main()


        