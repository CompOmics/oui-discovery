#!/usr/bin/env python3
import argparse
import sys
import logging
from pathlib import Path

#example
#python msfragger_generate_commands.py \
#    --working-dir /path/to/mzXML/PXD014258_mzXML \
#    --search-tag openprot \
#    --msfragger /path/to/MSFragger-4.2/MSFragger-4.2.jar \
#    --msfr-params /path/to/msfragger_params/open_fragger_OUIdiscovery_entrapment_openprot.params \
#    --java-xmx 100G \
#    --slurm-account X \
#    --slurm-time 06:00:00 \
#    --slurm-cpu 10 \
#    --slurm-mem 10G \
#    --out-dir /path/to/msfragger_commands \
#    --verbose


class MSFraggerGenSh:
    def __init__(self, mzxml_path: Path, out_dir: Path, search_tag: str,
                 msfragger_path: Path, msfr_params: Path, java_xmx: str,
                 slurm_account: str, slurm_time: str, slurm_cpu: str, slurm_mem: str):
        self.mzxml_path = mzxml_path
        self.out_dir = out_dir
        self.search_tag = search_tag
        self.msfragger_path = msfragger_path
        self.msfr_params = msfr_params
        self.java_xmx = java_xmx
        self.slurm_account = slurm_account
        self.slurm_time = slurm_time
        self.slurm_cpu = slurm_cpu
        self.slurm_mem = slurm_mem

    def write_sh(self):
        path = Path(self.mzxml_path)
        pxd = path.parent.name[:9]
        file_stem = path.stem
        out_sh = Path(self.out_dir) / f"{pxd}_{file_stem}_{self.search_tag}.sh"

        lines = [
            "#!/bin/bash\n",
            f"#SBATCH --account={self.slurm_account}\n",
            f"#SBATCH --time={self.slurm_time}\n",
            f"#SBATCH --cpus-per-task={self.slurm_cpu}\n",
            f"#SBATCH --mem-per-cpu={self.slurm_mem}\n",
            "module load java\n",
            f"java -Xmx{self.java_xmx} -jar {self.msfragger_path} {self.msfr_params} {self.mzxml_path}\n"
        ]

        try:
            with open(out_sh, "w") as f:
                f.writelines(lines)
            logging.info(f"Successfully wrote SLURM script to: {out_sh}")
            return True
        except Exception as e:
            logging.error(f"Failed to write SLURM script to: {out_sh}. Error: {e}")
            return False

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate .sh files to run MSFragger with slurm input."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--mzxml", "-f",
        type=Path,
        help="Path to a single mzXML or other MSFragger-suited file to search."
    )
    group.add_argument(
        "--working-dir", "-d",
        type=Path,
        help="Directory containing one or more mzXML or other MSFragger-suited files to search."
    )
    parser.add_argument(
        "--search-tag", "-tag", required=True,
        type=str,
        help="If the search is done several times on same files - put a tag for this search. It will appear in dire"
    )
    parser.add_argument(
        "--msfragger", "-msf", required=True,
        type=Path,
        help="Path to MSFragger on server."
    )
    parser.add_argument(
        "--msfr-params", "-p", required=True,
        type=Path,
        help="Path to MSFragger parameters file on server."
    )
    parser.add_argument(
        "--java-xmx", "-xmx", required=True,
        type=str,
        help="Xmx parameter for java with weight units (e.g. 100G)."
    )
    parser.add_argument(
        "--slurm-account", "-a", required=True,
        type=str,
        help="Slurm parameter --account (e.g. def-pi1)."
    )
    parser.add_argument(
        "--slurm-time", "-t", required=True,
        type=str,
        help="Slurm parameter --time (e.g. 1:00:00)."
    )
    parser.add_argument(
        "--slurm-cpu", "-cpu", required=True,
        type=str,
        help="Slurm parameter --cpus-per-task (e.g. 20)."
    )
    parser.add_argument(
        "--slurm-mem", "-mem", required=True,
        type=str,
        help="Slurm parameter --mem-per-cpu (e.g. 10G)."
    )
    parser.add_argument(
        "--out-dir", "-out", required=True,
        type=Path,
        help="Directory where to save .sh files."
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging."
    )

    return parser.parse_args()

def find_files(directory: Path) -> list[Path]:
    return sorted(directory.glob("*"))

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

    if args.mzxml:
        if not args.mzxml.exists():
            logging.error("Specified MSMS file does not exist: %s", args.mzxml)
            sys.exit(1)
        log_dir = args.out_dir.parent
        mzxml_paths = [args.mzxml]
    else:
        if not args.working_dir.is_dir():
            logging.error("Working directory is not a directory or does not exist: %s", args.working_dir)
            sys.exit(1)
        log_dir = args.out_dir
        mzxml_paths = find_files(args.working_dir)
        if not mzxml_paths:
            logging.error("No files found in directory: %s", args.working_dir)
            sys.exit(1)

    log_path = log_dir / "gen_msfr_sh.log"
    setup_logging(args.verbose, log_file=log_path)

    logging.debug("Parsed arguments: %s", args)

    if not args.out_dir.exists():
        logging.info("Output directory does not exist. Creating: %s", args.out_dir)
        args.out_dir.mkdir(parents=True, exist_ok=True)
    elif not args.out_dir.is_dir():
        logging.error("Specified output path exists but is not a directory: %s", args.out_dir)
        sys.exit(1)  

    any_failed = False
    for mzxml_file in mzxml_paths:
        try:
            gen = MSFraggerGenSh(
                mzxml_path=mzxml_file,
                out_dir=args.out_dir,
                search_tag=args.search_tag,
                msfragger_path=args.msfragger,
                msfr_params=args.msfr_params,
                java_xmx=args.java_xmx,
                slurm_account=args.slurm_account,
                slurm_time=args.slurm_time,
                slurm_cpu=args.slurm_cpu,
                slurm_mem=args.slurm_mem,
            )
            gen_res=gen.write_sh()
            if gen_res is False:
                logging.error("Saving .sh file failed for %s; skipping.", mzxml_file)
                any_failed = True
                continue
        except Exception:
            logging.exception("Exception while writing .sh for '%s'", mzxml_file)
            any_failed = True

    if any_failed:
        logging.error("One or more preprocessing files failed.")
        sys.exit(2)

if __name__ == "__main__":
    main()
